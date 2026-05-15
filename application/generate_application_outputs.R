# -----------------------------------------------------------
# ROBU: Generate Plots and Tables for k-Sensitivity Analysis 
# -----------------------------------------------------------

library(ggplot2)
library(dplyr)
library(patchwork) 
library(scales)

cat("\n--- Generating k-Sensitivity Analysis Results ---\n")

if (!dir.exists("simulations/figures")) {
  dir.create("simulations/figures", recursive = TRUE)
}

# 1. Load the data 
result_files <- list.files("simulations/results", pattern = "k_sensitivity_k.*\\.rds$", full.names = TRUE)
if(length(result_files) == 0) stop("No sensitivity simulation files found! Check the directory.")

results <- bind_rows(lapply(result_files, readRDS)) |> distinct() 
results$K_factor <- as.factor(results$K)

# Set epsilon target 
eps <- 0.20
plot_data <- results |> filter(Contamination == eps)

# __________________________
# 2. Print Summary Table 
# __________________________

cat(sprintf("\n--- Sensitivity Table Data (eps = %.2f) ---\n", eps))

summary_table <- plot_data |>
  group_by(K, BlockSize) |>
  summarize(
    MSE = mean(MSE, na.rm = TRUE),
    Time = mean(Time, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  arrange(K) |>
  mutate(
    MSE = ifelse(MSE > 1000, 
                 formatC(round(MSE, 1), format = "f", big.mark = ",", digits = 1), 
                 sprintf("%.1f", MSE)),
    Time = sprintf("%.1f", Time)
  )

print(as.data.frame(summary_table), row.names = FALSE)

# _______________________________________
# 3. Generate High-End Publication Plots 
# _______________________________________

cat("\n--- Generating PDF Plots ---\n")

x_labels <- as.character(summary_table$K)

pub_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_blank(), 
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = "black"), 
    axis.text.y = element_text(size = 11, color = "black"),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5, linetype = "dashed"),
    axis.line = element_line(color = "black", linewidth = 0.6)
  )

# A. Plot Time vs K (Left Panel - Log Scale for exponential drop)
p_time <- ggplot(plot_data, aes(x = K_factor, y = Time)) +
  geom_boxplot(fill = "#E0F3F8", color = "#313695", outlier.shape = NA, width = 0.6, linewidth = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5, color = "#313695", shape = 16) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1, big.mark = ","),
    breaks = scales::trans_breaks("log10", function(x) 10^x, n = 5)
  ) +
  scale_x_discrete(labels = x_labels) +
  labs(
    x = "Number of Blocks (k)",
    y = "Computation Time (Seconds, Log Scale)"
  ) +
  pub_theme

# B. Plot MSE vs K (Right Panel)
p_mse <- ggplot(plot_data, aes(x = K_factor, y = MSE)) +
  geom_boxplot(fill = "#FDDBC7", color = "#D55E00", outlier.shape = NA, width = 0.6, linewidth = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5, color = "#D55E00", shape = 16) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1, big.mark = ","),
    breaks = scales::trans_breaks("log10", function(x) 10^x, n = 5)
  ) +
  scale_x_discrete(labels = x_labels) +
  labs(
    x = "Number of Blocks (k)",
    y = "Mean Squared Error (Log Scale)"
  ) +
  pub_theme

# Combine the two plots side-by-side using patchwork
combined_plot <- p_time + p_mse + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold", size = 16))

plot_file <- "simulations/figures/k_sensitivity_plot_20.pdf"
ggsave(filename = plot_file, plot = combined_plot, width = 11, height = 5.5, dpi = 300)

cat(sprintf("Saved: %s\n", plot_file))