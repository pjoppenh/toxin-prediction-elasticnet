#' Resolve the data directory.
#'
get_data_dir <- function(default = "data") {
  data_dir <- Sys.getenv("TOXIN_DATA_DIR", unset = NA_character_)

  if (is.na(data_dir) || !nzchar(data_dir)) {
    data_dir <- default
  }

  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = FALSE)

  if (!dir.exists(data_dir)) {
    stop(
      "Data directory does not exist: ", data_dir, "\n",
      "Set it with: Sys.setenv(TOXIN_DATA_DIR = '/path/to/data')"
    )
  }

  return(data_dir)
}

# Helper to build file paths safely
path_in_data <- function(..., data_dir = get_data_dir()) {
  file.path(data_dir, ...)
}
