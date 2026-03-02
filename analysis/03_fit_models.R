# analysis/03_fit_models.R
library(glmnet)
library(foreach)
library(doParallel)

stopifnot(
  exists("toxins"),
  exists("tox_names"),
  exists("ASV_abundance_scaled"),
  exists("qPCR_scaled"),
  exists("relative_abundance_scaled"),
  "Species" %in% names(toxins)
)

# ---------------------------
# PARAMETERS
# ---------------------------
n_runs    <- 25
n_folds   <- 12
alpha     <- 0.5
num_cores <- 4   # adjust as needed

set.seed(130)

# ---------------------------
# Helper: grouped folds (keeps _1/_2/_3 together)
# ---------------------------
make_grouped_folds <- function(sample_id_vector, n_folds) {
  bio_ids <- sub("_(\\d+)$", "", sample_id_vector)
  unique_ids <- unique(bio_ids)

  # assign folds at biological-sample level
  fold_assignments <- sample(rep(seq_len(n_folds), length.out = length(unique_ids)))
  fold_map <- setNames(fold_assignments, unique_ids)

  fold_id <- as.integer(fold_map[bio_ids])

  # defensive
  if (anyNA(fold_id)) stop("Fold assignment produced NA values. Check sample IDs.")
  fold_id
}

sample_ids <- toxins$Species

model_inputs <- list(
  ASV  = ASV_abundance_scaled,
  qPCR = qPCR_scaled,
  RA   = relative_abundance_scaled
)

# ---------------------------
# Parallel cluster
# ---------------------------
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results <- foreach(i = seq_along(tox_names), .packages = "glmnet") %dopar% {

  toxin_name <- tox_names[i]
  toxin_vec  <- toxins[[toxin_name]]

  # ---- guard against NA/negative values ----
  if (anyNA(toxin_vec)) {
    stop(paste0("Toxin column '", toxin_name, "' contains NA values."))
  }
  if (any(toxin_vec < 0)) {
    stop(paste0("Toxin column '", toxin_name, "' has negative values; log(toxin+1) will produce NaNs."))
  }

  y_log <- log(toxin_vec + 1)
  y_scaled <- as.numeric(scale(y_log, center = TRUE, scale = TRUE))

  lambda_store <- list()
  rmse_store   <- list()
  dev_store    <- list()

  for (model_name in names(model_inputs)) {

    X <- model_inputs[[model_name]]

    lambda_mins <- numeric(n_runs)
    rmse_vals   <- numeric(n_runs)

    for (run in seq_len(n_runs)) {

      fold_id <- make_grouped_folds(sample_ids, n_folds)

      cvfit <- cv.glmnet(
        x = X,
        y = y_scaled,
        alpha = alpha,
        family = "gaussian",
        foldid = fold_id,
        type.measure = "mse"
      )

      lambda_mins[run] <- cvfit$lambda.min

      preds <- predict(cvfit, s = cvfit$lambda.min, newx = X)
      preds <- as.numeric(preds)

      # RMSE on log scale (unscaled y)
      rmse_vals[run] <- sqrt(mean((y_log - preds)^2))
    }

    final_lambda <- median(lambda_mins)

    final_model <- glmnet(
      x = X,
      y = y_scaled,
      alpha = alpha,
      family = "gaussian",
      lambda = final_lambda
    )

    lambda_store[[model_name]] <- lambda_mins
    rmse_store[[model_name]]   <- rmse_vals
    dev_store[[model_name]]    <- as.numeric(final_model$dev.ratio)
  }

  list(
    toxin_name = toxin_name,
    lambda     = lambda_store,
    rmse       = rmse_store,
    dev_ratio  = dev_store
  )
}

stopCluster(cl)
message("Model fitting complete.")
