# analysis/02_preprocess.R
stopifnot(
  exists("toxins"),
  exists("ASV_abundance_3"),
  exists("qPCR"),
  exists("relative_abundance")
)

# ---- toxins: multiple columns ----
stopifnot("Species" %in% names(toxins))
tox_names <- setdiff(names(toxins), "Species")
if (length(tox_names) == 0) stop("No toxin columns found in `toxins` besides Species.")

# ---- alignment checks (strict; assumes same row order) ----
stopifnot("Species" %in% names(ASV_abundance_3))
stopifnot("Species" %in% names(relative_abundance))
stopifnot("Species" %in% names(qPCR))

stopifnot(identical(toxins$Species, ASV_abundance_3$Species))
stopifnot(identical(toxins$Species, relative_abundance$Species))
stopifnot(identical(toxins$Species, qPCR$Species))

# ---- helper: numeric feature matrix excluding ID col ----
numeric_matrix_excluding <- function(df, drop_cols) {
  x_df <- df[, setdiff(names(df), drop_cols), drop = FALSE]
  is_num <- vapply(x_df, is.numeric, logical(1))
  x_df <- x_df[, is_num, drop = FALSE]
  if (ncol(x_df) == 0) stop("No numeric feature columns found after dropping: ", paste(drop_cols, collapse = ", "))
  as.matrix(x_df)
}

# Build & scale predictors (exclude ID columns only)
ASV_abundance_scaled        <- scale(numeric_matrix_excluding(ASV_abundance_3, drop_cols = c("Species")))
relative_abundance_scaled   <- scale(numeric_matrix_excluding(relative_abundance, drop_cols = c("Species")))
qPCR_scaled                 <- scale(numeric_matrix_excluding(qPCR, drop_cols = c("Species")))

# ---- messages ----
message("Preprocess complete:")
message("  tox_names: ", length(tox_names), " (", paste(tox_names, collapse = ", "), ")")
message("  ASV features: ", ncol(ASV_abundance_scaled))
message("  RA features:  ", ncol(relative_abundance_scaled))
message("  qPCR features:", ncol(qPCR_scaled))
message("  samples (rows): ", nrow(ASV_abundance_scaled))
