# -----------------------------------------------------------
# ROBU: Supplementary S4 - Severe Contamination Plots (35%)
# -----------------------------------------------------------

library(ggplot2)
library(dplyr)
library(gridExtra)
library(scales)

cat("\n--- Generating Supplementary S4 Severe Contamination Plots ---\n")

# Ensure the figures directory exists
if (!dir.exists("simulations/figures")) {
  dir.create("simulations/figures", recursive = TRUE)
}

# 1. Load the data
results_file <- "simulations/results/severe_contamination_eps0.35.rds"
if(!file.exists(results_file)) {
  stop("Result file not found! Please run 'run_severe_contamination.R' first.")
}
results <- readRDS(results_file)

# Lock in factor levels so plotting order is MM -> ROBU (OLS removed)
results$Method <- factor(results$Method, levels = c("Standard MM", "ROBU"))

# _______________________________________
# 2. Generate High-End Publication Plots 
# _______________________________________

# Professional Theme
pub_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = "black", face="bold"), 
    axis.text.y = element_text(size = 11, color = "black"),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "none" # Hide legend since x-axis labels are clear
  )

# Define custom colors for the methods (OLS removed)
method_colors <- c("Standard MM"  = "#D55E00",  # Vermilion/Red (Failure)
                   "ROBU"         = "#313695")  # Deep Blue (Success)

method_fills <- c("Standard MM"  = "#FDDBC7", 
                  "ROBU"         = "#E0F3F8")

# A. Plot Time vs Method (Left Panel)
p_time <- ggplot(results, aes(x = Method, y = Time, color = Method, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5, shape = 16) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_fills) +
  labs(
    x = "",
    y = "Computation Time (Seconds)",
    title = "Computational Effort"
  ) +
  pub_theme

# B. Plot Log-MSE vs Method (Right Panel)
p_mse <- ggplot(results, aes(x = Method, y = MSE, color = Method, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5, shape = 16) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_fills) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1, big.mark = ","),
    breaks = c(0.1, 1, 10, 100, 1000)
  ) +
  labs(
    x = "",
    y = "Mean Squared Error (Log Scale)",
    title = "Empirical Breakdown Resistance"
  ) +
  pub_theme

# Combine the two plots side-by-side
combined_plot <- grid.arrange(p_time, p_mse, ncol = 2)

# Save to PDF matching the LaTeX inclusion name
plot_file <- "simulations/figures/Supplementary_S4_Extreme_Contamination.pdf"

ggsave(filename = plot_file, plot = combined_plot, width = 10, height = 5.5, dpi = 300)
cat(sprintf("Saved Supplementary S4 figure to: %s\n", plot_file))