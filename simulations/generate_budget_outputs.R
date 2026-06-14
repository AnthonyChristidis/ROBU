# -------------------------------------------------------------------------
# Supplementary Study: ROBU Blockwise Resampling Budget (Visualization)
# -------------------------------------------------------------------------

library(ggplot2)

# Ensure figures directory exists
dir.create("simulations/figures", recursive = TRUE, showWarnings = FALSE)

# Load the raw results
raw_results_file <- "simulations/results/robu_budget_supplement_raw.rds"

if (!file.exists(raw_results_file)) {
  stop("Results file not found. Please run simulations/run_budget_study.R first.")
}

results <- readRDS(raw_results_file)

cat("Generating plots...\n")

# Format epsilon for nicer x-axis labels
results$Eps_Factor <- factor(results$Eps, levels = c(0.20, 0.30), 
                             labels = c("0.20", "0.30"))

# Set a colorblind-friendly palette
cb_palette <- c("ROBU (Full Budget)" = "#0072B2", 
                "ROBU (Distributed Budget)" = "#D55E00")

# __________________________
# Plot 1: MSE Comparison
# __________________________

p_mse <- ggplot(results, aes(x = Eps_Factor, y = MSE, fill = Method)) +
  geom_boxplot(alpha = 0.85, outlier.size = 1.2, width = 0.45, 
               position = position_dodge(width = 0.65), color = "black") +
  scale_y_log10() +
  scale_fill_manual(values = cb_palette) +
  labs(x = expression("Contamination Level (" * epsilon * ")"),
       y = "Mean Squared Error (Log Scale)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    axis.title = element_text(face = "bold", margin = margin(t = 10, r = 10)),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

# __________________________
# Plot 2: Time Comparison
# __________________________

p_time <- ggplot(results, aes(x = Eps_Factor, y = Time, fill = Method)) +
  geom_boxplot(alpha = 0.85, outlier.size = 1.2, width = 0.45, 
               position = position_dodge(width = 0.65), color = "black") +
  scale_fill_manual(values = cb_palette) +
  labs(x = expression("Contamination Level (" * epsilon * ")"),
       y = "Computation Time (Seconds)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    axis.title = element_text(face = "bold", margin = margin(t = 10, r = 10)),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

# __________________________
# Save to PDF
# __________________________

ggsave("simulations/figures/robu_budget_supplement_mse.pdf", plot = p_mse, width = 7, height = 5.5, dpi = 300)
ggsave("simulations/figures/robu_budget_supplement_time.pdf", plot = p_time, width = 7, height = 5.5, dpi = 300)

cat("Done! Plots saved to 'simulations/figures/'.\n")