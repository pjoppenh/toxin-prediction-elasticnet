library(foreach)
library(doParallel)
library(glmnet)
library(tidyverse)

# Set the number of cores
num_cores <- 20

# Load csv file
path <- ifelse(Sys.getenv("USERNAME") == "qdread", 'data/oppenheimer/ssim_jan2024', '/mnt/lili/pjo/rawdata')

# Toxins contain concentrations of all toxins of interest in each sample
toxins <- read.csv(file = file.path(path, "all_toxin_italy_2.csv"))
# Contains counts of each ASV of interest in each sample.
ASV_abundance_3 <- read.csv(file = file.path(path, "culled_asvs_nochim_2.csv"))
# Contains qPCR estimates of biomass of 8 different Fusarium species in each sample
qPCR <- read.csv(file = "/mnt/lili/pjo/rawdata/Italy_qPCRs.csv")

relative_abundance<-read.csv(file = "/mnt/lili/pjo/rawdata/Italy_relative_abundance.csv")

ASV_abundance_3$Toxin <- NA
ASV_abundance_3 <- ASV_abundance_3[, c(1, which(names(ASV_abundance_3) == "Toxin"), 2:(ncol(ASV_abundance_3) - 1))]
relative_abundance$Toxin <- NA
relative_abundance <- relative_abundance[, c(1, which(names(relative_abundance) == "Toxin"), 2:(ncol(relative_abundance) - 1))]

ASV_abundance_scaled <- scale(as.matrix(ASV_abundance_3[, -(1:2)]), center = TRUE, scale = TRUE)
qPCR_scaled <- scale(as.matrix(qPCR[, -(1)]), center = TRUE, scale = TRUE)
relative_abundance_scaled <- scale(as.matrix(relative_abundance[, -(1:2)]), center = TRUE, scale = TRUE)

n_runs <- 1000
n_folds <- 12

set.seed(130)

tox_names <- colnames(toxins)

# Create a cluster for parallel processing
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define the foreach loop
results <- foreach(i = 1:length(tox_names), .packages = c("glmnet")) %dopar% {
  ASV_abundance_3$Toxin <- toxins[, which(tox_names[i] == colnames(toxins))]

  lambda_mins_ASV <- rep(NA, n_runs)
  lambda_mins_qPCR <- rep(NA, n_runs)
  rmse_values_ASV <- rep(NA, n_runs)
  rmse_values_qPCR <- rep(NA, n_runs)

  for (run in 1:n_runs) {

    fold_id <- rep(sample(rep(c(1:n_folds), 2)), each = 3)
    alpha0.5fit_ASV <- cv.glmnet(
      x = ASV_abundance_scaled,
      y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
      type.measure = "mse",
      foldid = fold_id,  # Use the same fold_id for each run
      alpha = 0.5,
      family = "gaussian"
    )
    lambda_min_ASV <- alpha0.5fit_ASV$lambda.min
    lambda_mins_ASV[run] <- lambda_min_ASV

    predictions_ASV <- predict(alpha0.5fit_ASV, s = lambda_min_ASV, newx = ASV_abundance_scaled)
    rmse_values_ASV[run] <- sqrt(mean((log(ASV_abundance_3$Toxin + 1) - predictions_ASV)^2))

    alpha0.5fit_qPCR <- cv.glmnet(
      x = qPCR_scaled,
      y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
      type.measure = "mse",
      foldid = fold_id,  # Use the same fold_id for each run
      alpha = 0.5,
      family = "gaussian"
    )
    lambda_min_qPCR <- alpha0.5fit_qPCR$lambda.min
    lambda_mins_qPCR[run] <- lambda_min_qPCR
    predictions_qPCR <- predict(alpha0.5fit_qPCR, s = lambda_min_qPCR, newx = qPCR_scaled)
    rmse_values_qPCR[run] <- sqrt(mean((log(ASV_abundance_3$Toxin + 1) - predictions_qPCR)^2))
  }

  # Calculate the median of lambda_mins

  final_lambda_ASV <- median(lambda_mins_ASV)
  final_lambda_qPCR <- median(lambda_mins_qPCR)
  lambda2.5_ASV <- quantile(lambda_mins_ASV, 0.025)
  lambda97.5_ASV <- quantile(lambda_mins_ASV, 0.975)
  lambda2.5_qPCR <- quantile(lambda_mins_qPCR, 0.025)
  lambda97.5_qPCR <- quantile(lambda_mins_qPCR, 0.975)
  model0.5fit_ASV <- glmnet(
    x = ASV_abundance_scaled,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = final_lambda_ASV
  )
  model0.5fit_qPCR <- glmnet(
    x = qPCR_scaled,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = final_lambda_qPCR
  )
  model0.025fit_ASV <- glmnet(
    x = ASV_abundance_scaled,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = lambda2.5_ASV
  )
  model0.025fit_qPCR <- glmnet(
    x = qPCR_scaled,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = lambda2.5_qPCR
  )
  model0.975fit_ASV <- glmnet(
    x = ASV_abundance_scaled,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = lambda97.5_ASV
  )
  model0.975fit_qPCR <- glmnet(
    x = qPCR_scaled,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = lambda97.5_qPCR
  )
  # Perform t-test to compare RMSE values
  if (var(rmse_values_ASV) > 0 && var(rmse_values_qPCR) > 0) {
    t_test_result_qPCR <- t.test(rmse_values_ASV, rmse_values_qPCR)
  } else {
    t_test_result_qPCR <- list(
      statistic = NA,
      p.value = NA,
      conf.int = NA,
      estimate = NA,
      null.value = NA,
      alternative = NA,
      method = NA,
      data.name = NA
    )
  }

  # Store the results for this toxin
  list(
    rsquared_values_ASV = model0.5fit_ASV$dev.ratio,
    rsquared_values_qPCR = model0.5fit_qPCR$dev.ratio,
    rsquared_values2.5_ASV =  model0.025fit_ASV$dev.ratio,
    rsquared_values2.5_qPCR =  model0.025fit_qPCR$dev.ratio,
    rsquared_values97.5_ASV =  model0.975fit_ASV$dev.ratio,
    rsquared_values97.5_qPCR =  model0.975fit_qPCR$dev.ratio,
    variable_coefs_ASV = as.data.frame(coef(model0.5fit_ASV)[, 1]),
    variable_coefs_qPCR = as.data.frame(coef(model0.5fit_qPCR)[, 1]),
    rmse_values_ASV = rmse_values_ASV,
    rmse_values_qPCR = rmse_values_qPCR,
    t_test_result_qPCR = t_test_result_qPCR
  )
}

# Stop the parallel backend
stopCluster(cl)

cl <- makeCluster(num_cores)
registerDoParallel(cl)


results_2 <- foreach(i = 1:length(tox_names), .packages = c("glmnet")) %dopar% {
  ASV_abundance_3$Toxin <- toxins[, which(tox_names[i] == colnames(toxins))]
  relative_abundance$Toxin <- toxins[, which(tox_names[i] == colnames(toxins))]

  lambda_mins_ASV <- rep(NA, n_runs)
  lambda_mins_RA <- rep(NA, n_runs)
  rmse_values_ASV <- rep(NA, n_runs)
  rmse_values_RA <- rep(NA, n_runs)

  for (run in 1:n_runs) {

    fold_id <- rep(sample(rep(c(1:n_folds), 2)), each = 3)
    alpha0.5fit_ASV <- cv.glmnet(
      x = ASV_abundance_scaled,
      y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
      type.measure = "mse",
      foldid = fold_id,  # Use the same fold_id for each run
      alpha = 0.5,
      family = "gaussian"
    )
    lambda_min_ASV <- alpha0.5fit_ASV$lambda.min
    lambda_mins_ASV[run] <- lambda_min_ASV
    predictions_ASV <- predict(alpha0.5fit_ASV, s = lambda_min_ASV, newx = ASV_abundance_scaled)
    rmse_values_ASV[run] <- sqrt(mean((log(ASV_abundance_3$Toxin + 1) - predictions_ASV)^2))

    alpha0.5fit_RA <- cv.glmnet(
      x = relative_abundance_scaled ,
      y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
      type.measure = "mse",
      foldid = fold_id,  # Use the same fold_id for each run
      alpha = 0.5,
      family = "gaussian"
    )
    lambda_min_RA <- alpha0.5fit_RA$lambda.min
    lambda_mins_RA[run] <- lambda_min_RA
    predictions_RA <- predict(alpha0.5fit_RA, s = lambda_min_RA, newx = relative_abundance_scaled )
    rmse_values_RA[run] <- sqrt(mean((log(relative_abundance$Toxin + 1) - predictions_RA)^2))
  }

  # Calculate the median of lambda_mins
  final_lambda_ASV <- median(lambda_mins_ASV)
  final_lambda_RA <- median(lambda_mins_RA)
  lambda2.5_ASV <- quantile(lambda_mins_ASV, 0.025)
  lambda97.5_ASV <- quantile(lambda_mins_ASV, 0.975)
  lambda2.5_RA <- quantile(lambda_mins_RA, 0.025)
  lambda97.5_RA <- quantile(lambda_mins_RA, 0.975)
  model0.5fit_ASV <- glmnet(
    x = ASV_abundance_scaled,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = final_lambda_ASV
  )
  model0.5fit_RA <- glmnet(
    x = relative_abundance_scaled ,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = final_lambda_RA
  )
  model0.025fit_ASV <- glmnet(
    x = ASV_abundance_scaled,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = lambda2.5_ASV
  )
  model0.025fit_RA <- glmnet(
    x = relative_abundance_scaled ,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = lambda2.5_RA
  )
  model0.975fit_ASV <- glmnet(
    x = ASV_abundance_scaled,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = lambda97.5_ASV
  )
  model0.975fit_RA <- glmnet(
    x = relative_abundance_scaled ,
    y = scale(as.matrix(log(ASV_abundance_3$Toxin + 1)), center = TRUE, scale = TRUE),
    alpha = 0.5,
    family = "gaussian",
    lambda = lambda97.5_RA
  )

  # Perform t-test to compare RMSE values
  if (var(rmse_values_ASV) > 0 && var(rmse_values_RA) > 0) {
    t_test_result_RA <- t.test(rmse_values_ASV, rmse_values_RA)
  } else {
    t_test_result_RA <- list(
      statistic = NA,
      p.value = NA,
      conf.int = NA,
      estimate = NA,
      null.value = NA,
      alternative = NA,
      method = NA,
      data.name = NA
    )
  }

  # Store the results for this toxin
  list(
    rsquared_values_RA = model0.5fit_RA$dev.ratio,
    rsquared_values2.5_ASV =  model0.025fit_ASV$dev.ratio,
    rsquared_values2.5_RA =  model0.025fit_RA$dev.ratio,
    rsquared_values97.5_ASV =  model0.975fit_ASV$dev.ratio,
    rsquared_values97.5_RA =  model0.975fit_RA$dev.ratio,
    variable_coefs_RA = as.data.frame(coef(model0.5fit_RA)[, 1]),
    rmse_values_RA = rmse_values_RA,
    t_test_result_RA = t_test_result_RA
  )
}

# Stop the parallel backend
stopCluster(cl)



# Aggregate results
ASV_coefs <- list()
ASV_rsquared_values_toxins <- list()
ASV_rsquared_values_toxins_2.5 <- list()
ASV_rsquared_values_toxins_97.5 <- list()
for (i in 1:length(results)) {
  ASV_coefs[[i]] <- do.call(cbind, results[[i]]$variable_coefs_ASV)
  ASV_rsquared_values_toxins[[i]] <- results[[i]]$rsquared_values_ASV
  ASV_rsquared_values_toxins_2.5[[i]] <- results[[i]]$rsquared_values2.5_ASV
  ASV_rsquared_values_toxins_97.5[[i]] <- results[[i]]$rsquared_values97.5_ASV
}

# Combine the list of data frames into a single data frame
combined_ASV_coefs <- do.call(cbind, ASV_coefs)

# Optionally, set column names based on the names of the original list
colnames(combined_ASV_coefs) <- tox_names
combined_ASV_coefs<-combined_ASV_coefs[-1,]
rownames(combined_ASV_coefs) <- colnames(ASV_abundance_3[,3:52])
# View the structure of the combined data frame
str(combined_ASV_coefs)

str(ASV_rsquared_values_toxins)

# Extract the R^2 values from the list
rsquared_values_ASV <- unlist(lapply(ASV_rsquared_values_toxins, function(x) unname(unlist(x))))
rsquared_values_ASV_2.5 <- unlist(lapply(ASV_rsquared_values_toxins_2.5, function(x) unname(unlist(x))))
rsquared_values_ASV_97.5 <- unlist(lapply(ASV_rsquared_values_toxins_97.5, function(x) unname(unlist(x))))
# Create a data frame with one column for R^2 values
rsquared_df_ASV <- data.frame(RSquared_ASV = rsquared_values_ASV, row.names = tox_names, RSquared_ASV_2.5 = rsquared_values_ASV_2.5, RSquared_ASV_97.5 = rsquared_values_ASV_97.5)

# View the structure of the data frame
str(rsquared_df_ASV)



# Aggregate results
qPCR_coefs <- list()
qPCR_rsquared_values_toxins <- list()
qPCR_rsquared_values_toxins_2.5 <- list()
qPCR_rsquared_values_toxins_97.5 <- list()
for (i in 1:length(results)) {
  qPCR_coefs[[i]] <- do.call(cbind, results[[i]]$variable_coefs_qPCR)
  qPCR_rsquared_values_toxins[[i]] <- results[[i]]$rsquared_values_qPCR
  qPCR_rsquared_values_toxins_2.5[[i]] <- results[[i]]$rsquared_values2.5_qPCR
  qPCR_rsquared_values_toxins_97.5[[i]] <- results[[i]]$rsquared_values97.5_qPCR
}

# Combine the list of data frames into a single data frame
combined_qPCR_coefs <- do.call(cbind, qPCR_coefs)

# Optionally, set column names based on the names of the original list
colnames(combined_qPCR_coefs) <- tox_names
combined_qPCR_coefs<-combined_qPCR_coefs[-1,]
rownames(combined_qPCR_coefs) <- colnames(qPCR[,2:9])
# View the structure of the combined data frame
str(combined_qPCR_coefs)

str(qPCR_rsquared_values_toxins)

# Extract the R^2 values from the list
rsquared_values_qPCR <- unlist(lapply(qPCR_rsquared_values_toxins, function(x) unname(unlist(x))))
rsquared_values_qPCR_2.5 <- unlist(lapply(qPCR_rsquared_values_toxins_2.5, function(x) unname(unlist(x))))
rsquared_values_qPCR_97.5 <- unlist(lapply(qPCR_rsquared_values_toxins_97.5, function(x) unname(unlist(x))))
# Create a data frame with one column for R^2 values
rsquared_df_ASV <- data.frame(RSquared_ASV = rsquared_values_ASV, row.names = tox_names, RSquared_ASV_2.5 = rsquared_values_ASV_2.5, RSquared_ASV_97.5 = rsquared_values_ASV_97.5)

# Create a data frame with one column for R^2 values
rsquared_df_qPCR <- data.frame(RSquared_qPCR = rsquared_values_qPCR, row.names = tox_names, RSquared_qPCR_2.5 = rsquared_values_qPCR_2.5, RSquared_qPCR_97.5 = rsquared_values_qPCR_97.5)

# View the structure of the data frame
str(rsquared_df_qPCR)


# Aggregate results_2
RA_coefs <- list()
RA_rsquared_values_toxins <- list()
RA_rsquared_values_toxins_2.5 <- list()
RA_rsquared_values_toxins_97.5 <- list()
for (i in 1:length(results_2)) {
  RA_coefs[[i]] <- do.call(cbind, results_2[[i]]$variable_coefs_RA)
  RA_rsquared_values_toxins[[i]] <- results_2[[i]]$rsquared_values_RA
  RA_rsquared_values_toxins_2.5[[i]] <- results_2[[i]]$rsquared_values2.5_RA
  RA_rsquared_values_toxins_97.5[[i]] <- results_2[[i]]$rsquared_values97.5_RA
}

# Combine the list of data frames into a single data frame
combined_RA_coefs <- do.call(cbind, RA_coefs)

# Optionally, set column names based on the names of the original list
colnames(combined_RA_coefs) <- tox_names
combined_RA_coefs<-combined_RA_coefs[-1,]
rownames(combined_RA_coefs) <- colnames(relative_abundance[1,3:52])
# View the structure of the combined data frame
str(combined_RA_coefs)

str(RA_rsquared_values_toxins)

# Extract the R^2 values from the list
rsquared_values_RA <- unlist(lapply(RA_rsquared_values_toxins, function(x) unname(unlist(x))))
rsquared_values_RA_2.5 <- unlist(lapply(RA_rsquared_values_toxins_2.5, function(x) unname(unlist(x))))
rsquared_values_RA_97.5 <- unlist(lapply(RA_rsquared_values_toxins_97.5, function(x) unname(unlist(x))))
# Create a data frame with one column for R^2 values
rsquared_df_RA <- data.frame(RSquared_RA = rsquared_values_RA, row.names = tox_names, RSquared_RA_2.5 = rsquared_values_RA_2.5, RSquared_RA_97.5 = rsquared_values_RA_97.5)

# View the structure of the data frame
str(rsquared_df_RA)


write.csv(combined_ASV_coefs, file = '/mnt/lili/pjo/rawdata/combined_ASV_coefs_fin.csv')
write.csv(rsquared_df_ASV, file = "/mnt/lili/pjo/rawdata/rsquared_df_ASV_fin.csv")
write.csv(combined_qPCR_coefs, file = '/mnt/lili/pjo/rawdata/combined_qPCR_coefs_fin.csv')
write.csv(rsquared_df_qPCR, file = "/mnt/lili/pjo/rawdata/rsquared_df_qPCR_fin.csv")
write.csv(combined_RA_coefs, file = '/mnt/lili/pjo/rawdata/combined_RA_coefs_fin.csv')
write.csv(rsquared_df_RA, file = "/mnt/lili/pjo/rawdata/rsquared_df_RA_fin.csv")


saveRDS(ASV_coefs, file = "/mnt/lili/pjo/rawdata/ASV_coefs.rds")
saveRDS(ASV_rsquared_values_toxins, file = "/mnt/lili/pjo/rawdata/ASV_rsquared_values_toxins.rds")
saveRDS(qPCR_coefs, file = "/mnt/lili/pjo/rawdata/qPCR_coefs.rds")
saveRDS(qPCR_rsquared_values_toxins, file = "/mnt/lili/pjo/rawdata/qPCR_rsquared_values_toxins.rds")
saveRDS(RA_coefs, file = "/mnt/lili/pjo/rawdata/RA_coefs.rds")
saveRDS(RA_rsquared_values_toxins, file = "/mnt/lili/pjo/rawdata/RA_rsquared_values_toxins.rds")


# Initialize a list to store the t-test results
t_test_results_all_qPCR <- list()

# Initialize an empty datafqPCRme to store the final results
final_results_df_qPCR <- data.frame(
  Toxin = character(),
  p_value = numeric(),
  lower_conf_int = numeric(),
  upper_conf_int = numeric(),
  stringsAsFactors = FALSE
)

# Process results for each toxin
for (i in 1:length(tox_names)) {
  t_test_results_all_qPCR[[tox_names[i]]] <- results[[i]]$t_test_result_qPCR

  # ExtqPCRct the p-value and confidence intervals
  p_value <- t_test_results_all_qPCR[[tox_names[i]]]$p.value
  lower_conf_int <- t_test_results_all_qPCR[[tox_names[i]]]$conf.int[1]
  upper_conf_int <- t_test_results_all_qPCR[[tox_names[i]]]$conf.int[2]

  # Append the results to the datafqPCRme
  final_results_df_qPCR <- rbind(final_results_df_qPCR, data.frame(
    Toxin = tox_names[i],
    p_value = p_value,
    lower_conf_int = lower_conf_int,
    upper_conf_int = upper_conf_int
  ))
}

# Set the row names to be the Toxin names
rownames(final_results_df_qPCR) <- final_results_df_qPCR$Toxin
final_results_df_qPCR$Toxin <- NULL

# Print the final results dataframe
print(final_results_df_qPCR)



# Initialize a list to store the t-test results
t_test_results_all_RA <- list()

# Initialize an empty dataframe to store the final results
final_results_df_RA <- data.frame(
  Toxin = character(),
  p_value = numeric(),
  lower_conf_int = numeric(),
  upper_conf_int = numeric(),
  stringsAsFactors = FALSE
)

# Process results for each toxin
for (i in 1:length(tox_names)) {
  t_test_results_all_RA[[tox_names[i]]] <- results_2[[i]]$t_test_result_RA

  # Extract the p-value and confidence intervals
  p_value <- t_test_results_all_RA[[tox_names[i]]]$p.value
  lower_conf_int <- t_test_results_all_RA[[tox_names[i]]]$conf.int[1]
  upper_conf_int <- t_test_results_all_RA[[tox_names[i]]]$conf.int[2]

  # Append the results to the dataframe
  final_results_df_RA <- rbind(final_results_df_RA, data.frame(
    Toxin = tox_names[i],
    p_value = p_value,
    lower_conf_int = lower_conf_int,
    upper_conf_int = upper_conf_int
  ))
}

# Set the row names to be the Toxin names
rownames(final_results_df_RA) <- final_results_df_RA$Toxin
final_results_df$Toxin <- NULL

# Print the final results dataframe
print(final_results_df_RA)

write.csv(final_results_df_qPCR, file = '/mnt/lili/pjo/rawdata/final_results_df_qPCR.csv')
write.csv(final_results_df_RA, file = '/mnt/lili/pjo/rawdata/final_results_df_RA.csv')
saveRDS(final_results_df_qPCR, file = "/mnt/lili/pjo/rawdata/final_results_df_qPCR.rds")
saveRDS(final_results_df_RA, file = "/mnt/lili/pjo/rawdata/final_results_df_RA.rds")

