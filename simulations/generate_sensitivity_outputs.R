# -----------------------------------------------------------
# ROBU: Generate Plots and Tables for k-Sensitivity Analysis 
# -----------------------------------------------------------

library(ggplot2)
library(dplyr)
library(patchwork) 
library(scales)

cat("\n--- Generating k-Sensitivity Analysis Results ---\n")

# Ensure the figures directory exists
if (!dir.exists("simulations/figures")) {
  dir.create("simulations/figures", recursive = TRUE)
}

# 1. Load the data 
result_files <- list.files("simulations/results", pattern = "k_sensitivity_k.*\\.rds$", full.names = TRUE)
if(length(result_files) == 0) stop("No sensitivity simulation files found! Check the directory.")

results <- bind_rows(lapply(result_files, readRDS)) |> distinct() 

results$K_factor <- as.factor(results$K)

# Custom x-axis labels based on unique K values
x_labels <- as.character(sort(unique(results$K)))

# __________________________
# 2. Generate Summary Table 
# __________________________

# Set epsilon target
eps <- 0.20
plot_data <- results |> filter(Contamination == eps)

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
    MSE_fmt = ifelse(MSE > 1000, 
                 formatC(round(MSE, 1), format = "f", big.mark = ",", digits = 1), 
                 sprintf("%.1f", MSE)),
    Time_fmt = sprintf("%.1f", Time)
  )

print(as.data.frame(summary_table |> select(K, BlockSize, MSE=MSE_fmt, Time=Time_fmt)), row.names = FALSE)

# _______________________________________
# 3. Generate High-End Publication Plots 
# _______________________________________

cat("\n--- Generating PDF Plots ---\n")

pub_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_blank(), 
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = "black"), 
    axis.text.y = element_text(size = 11, color = "black"),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5, linetype = "dashed"),
    axis.line = element_line(color = "black", linewidth = 0.6)
  )

# The official ROBU color from Figure 1
robu_color <- "#0072B2"

# Data for Line Plot (Averages)
plot_data_line <- plot_data |>
  group_by(K_factor) |>
  summarize(Mean_Time = mean(Time, na.rm = TRUE), .groups = "drop")

# A. Panel 1 (Left): MSE vs K (Boxplot)
p_mse <- ggplot(plot_data, aes(x = K_factor, y = MSE)) +
  geom_boxplot(fill = robu_color, color = "black", alpha = 0.5, outlier.shape = NA, width = 0.6, linewidth = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5, color = robu_color, shape = 16) +
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

# B. Panel 2 (Right): Time vs K (Line Plot - Log Scale for exponential drop)
p_time <- ggplot(plot_data_line, aes(x = K_factor, y = Mean_Time, group = 1)) +
  geom_line(color = robu_color, linewidth = 1) +
  geom_point(color = robu_color, size = 3.5) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1, big.mark = ","),
    breaks = scales::trans_breaks("log10", function(x) 10^x, n = 5)
  ) +
  scale_x_discrete(labels = x_labels) +
  labs(
    x = "Number of Blocks (k)",
    y = "Average Computation Time (Seconds, Log Scale)"
  ) +
  pub_theme

# Combine the two plots side-by-side using patchwork (MSE on Left, Time on Right)
combined_plot <- p_mse + p_time 

plot_file <- "simulations/figures/k_sensitivity_plot_20.pdf"
ggsave(filename = plot_file, plot = combined_plot, width = 11, height = 5.5, dpi = 300)

cat(sprintf("Saved: %s\n", plot_file))