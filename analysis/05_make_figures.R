# analysis/05_make_figures.R

library(ggplot2)
library(dplyr)
library(tidyr)

stopifnot(file.exists("outputs/tables/dev_ratio_summary.csv"))
stopifnot(file.exists("outputs/tables/rmse_long.csv"))

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/figures", recursive = TRUE, showWarnings = FALSE)  # for README previews

# ---- Load data ----
dev_ratio_df <- read.csv("outputs/tables/dev_ratio_summary.csv")
rmse_long_df <- read.csv("outputs/tables/rmse_long.csv")

# ---- Consistent ordering ----
tox_order <- unique(dev_ratio_df$Toxin)
model_order <- c("ASV", "qPCR", "RA")

dev_ratio_df <- dev_ratio_df %>%
  mutate(
    Toxin = factor(Toxin, levels = tox_order),
    Model = factor(Model, levels = model_order)
  )

rmse_long_df <- rmse_long_df %>%
  mutate(
    Toxin = factor(Toxin, levels = tox_order),
    Model = factor(Model, levels = model_order)
  )

# ---- Clean theme ----
theme_clean <- theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 15, 10, 15)
  )

plot_width <- max(10, 0.9 * length(tox_order) + 6)

# =============================================================================
# 1) Deviance Ratio Plot (ensure zeros show + ensure all Toxin x Model combos)
# =============================================================================

dev_ratio_plot_df <- dev_ratio_df %>%
  tidyr::complete(Toxin, Model, fill = list(Dev_Ratio = 0))

# Tiny floor purely for visibility of true zeros
eps <- 0.005
dev_ratio_plot_df <- dev_ratio_plot_df %>%
  mutate(Dev_Ratio_plot = ifelse(is.na(Dev_Ratio) | Dev_Ratio == 0, eps, Dev_Ratio))

p_dev <- ggplot(dev_ratio_plot_df, aes(x = Toxin, y = Dev_Ratio_plot, fill = Model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = sprintf("%.2f", Dev_Ratio)),
    position = position_dodge(width = 0.75),
    vjust = -0.35,
    size = 3
  ) +
  coord_cartesian(ylim = c(0, max(dev_ratio_plot_df$Dev_Ratio_plot, na.rm = TRUE) * 1.15)) +
  theme_clean +
  labs(
    title = "Elastic Net Model Performance (Deviance Ratio)",
    subtitle = "Deviance explained (R²-equivalent) across diagnostic data types\n(Zero values shown as small ticks for visibility)",
    y = "Deviance Ratio",
    x = "Toxin",
    fill = "Model"
  )

ggsave("outputs/figures/dev_ratio_comparison.png", plot = p_dev, width = plot_width, height = 6, dpi = 300)
ggsave("docs/figures/dev_ratio_comparison.png",    plot = p_dev, width = plot_width, height = 6, dpi = 300)

# =============================================================================
# 2) RMSE Plot (boxplots across repeated runs)
# =============================================================================

p_rmse <- ggplot(rmse_long_df, aes(x = Toxin, y = RMSE, fill = Model)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    outlier.alpha = 0.25
  ) +
  theme_clean +
  labs(
    title = "Elastic Net RMSE Comparison",
    subtitle = "Distribution of RMSE across repeated grouped CV runs",
    y = "RMSE (log-scale response)",
    x = "Toxin",
    fill = "Model"
  )

ggsave("outputs/figures/rmse_comparison.png", plot = p_rmse, width = plot_width, height = 6, dpi = 300)
ggsave("docs/figures/rmse_comparison.png",    plot = p_rmse, width = plot_width, height = 6, dpi = 300)

message("Figures saved to outputs/figures/ and docs/figures/")
