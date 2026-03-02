# data-raw/01_generate_synthetic_data.R

set.seed(130)

dir.create("data", showWarnings = FALSE)

n_bio <- 60
n_rep <- 3
n <- n_bio * n_rep

Species <- paste0(rep(seq_len(n_bio), each = n_rep), "_", rep(seq_len(n_rep), times = n_bio))

p_asv <- 60
p_ra  <- 60
p_qpcr <- 8

# Feature matrices (example)
ASV <- matrix(rnorm(n * p_asv), nrow = n, ncol = p_asv)
colnames(ASV) <- paste0("ASV_", seq_len(p_asv))

RA <- matrix(rnorm(n * p_ra), nrow = n, ncol = p_ra)
colnames(RA) <- paste0("RA_", seq_len(p_ra))

qPCR <- matrix(rnorm(n * p_qpcr), nrow = n, ncol = p_qpcr)
colnames(qPCR) <- paste0("qPCR_", seq_len(p_qpcr))

# Build multiple toxins with different "true" drivers (so models differ)
signal_asv  <-  1.2 * ASV[, 1] - 0.8 * ASV[, 2] + 0.5 * ASV[, 3]
signal_qpcr <-  1.5 * qPCR[, 1] - 1.0 * qPCR[, 2]
signal_ra   <-  1.0 * RA[, 1] + 0.6 * RA[, 2] - 0.4 * RA[, 3]

noise <- function(sd = 0.6) rnorm(n, 0, sd)

toxins <- data.frame(
  Species = Species,
  Toxin_A = pmax(0, exp(0.8 * signal_asv  + noise()) - 1),
  Toxin_B = pmax(0, exp(0.8 * signal_qpcr + noise()) - 1),
  Toxin_C = pmax(0, exp(0.6 * signal_ra   + noise()) - 1),
  Toxin_D = pmax(0, exp(0.4 * signal_asv + 0.4 * signal_qpcr + noise()) - 1)
)

# Write CSVs in the same filenames your pipeline expects
write.csv(toxins, "data/all_toxin_italy_2.csv", row.names = FALSE)

write.csv(data.frame(Species = Species, ASV), "data/culled_asvs_nochim_2.csv", row.names = FALSE)
write.csv(data.frame(Species = Species, qPCR), "data/Italy_qPCRs.csv", row.names = FALSE)
write.csv(data.frame(Species = Species, RA),  "data/Italy_relative_abundance.csv", row.names = FALSE)

message("Wrote synthetic inputs to ./data/")
