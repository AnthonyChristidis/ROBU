# -----------------------------------------------------------
# ROBU: Generate Plots and Tables for k-Sensitivity Analysis 
# -----------------------------------------------------------

library(ggplot2)
library(dplyr)
library(gridExtra)
library(scales)

cat("\n--- Generating k-Sensitivity Analysis Results ---\n")

# 1. Load the data
result_files <- list.files("simulations/results", pattern = "k_sensitivity_k.*\\.rds$", full.names = TRUE)
if(length(result_files) == 0) {
  result_files <- "simulations/results/k_sensitivity_results.rds"
}
results <- bind_rows(lapply(result_files, readRDS)) |> distinct() 

results$K_factor <- as.factor(results$K)

# __________________________
# 2. Generate Summary Table 
# __________________________

summary_table <- results |>
  group_by(K, BlockSize) |>
  summarize(
    Conv_Rate = mean(S_Conv, na.rm = TRUE) * 100, 
    .groups = "drop"
  ) |> arrange(K)

# Custom x-axis labels 
x_labels <- as.character(summary_table$K)

# _______________________________________
# 3. Generate High-End Publication Plots 
# _______________________________________

cat("\n--- Generating Professional PDF Plots ---\n")

# A. Professional Theme (Minimalist, academic, publication-ready)
pub_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = "black"), 
    axis.text.y = element_text(size = 11, color = "black"),
    # Subtle horizontal grid lines to make tracking values easier
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

# B. Plot Time vs K (Left Panel - Professional Blue)
p_time <- ggplot(results, aes(x = K_factor, y = Time)) +
  # Boxplot with outliers hidden (since jitter will show them)
  geom_boxplot(fill = "#E0F3F8", color = "#313695", outlier.shape = NA, width = 0.6, linewidth = 0.6) +
  # Jittered raw points on top
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5, color = "#313695", shape = 16) +
  scale_x_discrete(labels = x_labels) +
  labs(
    x = "Number of Blocks (k)",
    y = "Computation Time (Seconds)",
    title = "Computational Scaling"
  ) +
  pub_theme

# C. Plot MSE vs K (Right Panel - Okabe-Ito Vermilion/Rust)
p_mse <- ggplot(results, aes(x = K_factor, y = MSE)) +
  # Boxplot with very soft, muted peach/rust fill
  geom_boxplot(fill = "#FDDBC7", color = "#D55E00", outlier.shape = NA, width = 0.6, linewidth = 0.6) +
  # Jittered raw points in Deep Vermilion
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5, color = "#D55E00", shape = 16) +
  # Log scale with clean formatting
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1, big.mark = ","),
    breaks = c(0.1, 1, 10, 100, 1000)
  ) +
  scale_x_discrete(labels = x_labels) +
  labs(
    x = "Number of Blocks (k)",
    y = "Mean Squared Error (Log Scale)",
    title = "Empirical Breakdown Resistance"
  ) +
  pub_theme

# Combine the two plots side-by-side
combined_plot <- grid.arrange(p_time, p_mse, ncol = 2)

# Save to PDF
plot_file <- "simulations/figures/k_sensitivity_plot.pdf"
ggsave(filename = plot_file, plot = combined_plot, width = 12, height = 5.5, dpi = 300)

cat(sprintf("Sensitivity plots successfully saved to: %s\n", plot_file))