# analysis/01_load_data.R
source("R/config.R")

data_dir <- get_data_dir()

toxins <- read.csv(path_in_data("all_toxin_italy_2.csv", data_dir = data_dir), stringsAsFactors = FALSE)
ASV_abundance_3 <- read.csv(path_in_data("culled_asvs_nochim_2.csv", data_dir = data_dir), stringsAsFactors = FALSE)
qPCR <- read.csv(path_in_data("Italy_qPCRs.csv", data_dir = data_dir), stringsAsFactors = FALSE)
relative_abundance <- read.csv(path_in_data("Italy_relative_abundance.csv", data_dir = data_dir), stringsAsFactors = FALSE)

#sanity checks
stopifnot(nrow(toxins) > 0, nrow(ASV_abundance_3) > 0, nrow(qPCR) > 0, nrow(relative_abundance) > 0)

message("Loaded files from: ", data_dir)
message("toxins: ", nrow(toxins), " x ", ncol(toxins))
message("ASV_abundance_3: ", nrow(ASV_abundance_3), " x ", ncol(ASV_abundance_3))
message("qPCR: ", nrow(qPCR), " x ", ncol(qPCR))
message("relative_abundance: ", nrow(relative_abundance), " x ", ncol(relative_abundance))
