# analysis/04_summarize_results.R
stopifnot(exists("results"))

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()
rmse_rows <- list()
rmse_long_rows <- list()  # optional: for boxplots

for (res in results) {

  toxin_name <- res$toxin_name

  for (model_type in names(res$dev_ratio)) {

    # ---- Deviance ratio ----
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Toxin = toxin_name,
      Model = model_type,
      Dev_Ratio = as.numeric(res$dev_ratio[[model_type]])
    )

    # ---- RMSE summaries ----
    rmse_vals <- as.numeric(res$rmse[[model_type]])
    rmse_vals <- rmse_vals[is.finite(rmse_vals)]

    rmse_rows[[length(rmse_rows) + 1]] <- data.frame(
      Toxin = toxin_name,
      Model = model_type,
      N_Runs = length(rmse_vals),
      RMSE_Mean = mean(rmse_vals),
      RMSE_SD = sd(rmse_vals),
      RMSE_Median = median(rmse_vals),
      RMSE_2.5 = as.numeric(quantile(rmse_vals, 0.025)),
      RMSE_97.5 = as.numeric(quantile(rmse_vals, 0.975))
    )

    # ---- RMSE long (optional) ----
    rmse_long_rows[[length(rmse_long_rows) + 1]] <- data.frame(
      Toxin = toxin_name,
      Model = model_type,
      RMSE = rmse_vals
    )
  }
}

dev_ratio_df <- do.call(rbind, summary_rows)
rmse_summary_df <- do.call(rbind, rmse_rows)
rmse_long_df <- do.call(rbind, rmse_long_rows)

write.csv(dev_ratio_df, "outputs/tables/dev_ratio_summary.csv", row.names = FALSE)
write.csv(rmse_summary_df, "outputs/tables/rmse_summary.csv", row.names = FALSE)
write.csv(rmse_long_df, "outputs/tables/rmse_long.csv", row.names = FALSE)

message("Summary tables written to outputs/tables/")
