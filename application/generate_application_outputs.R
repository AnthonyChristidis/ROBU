# ----------------------------------------------------------------
# ROBU: Generate Plots and Tables for Proteogenomics Application
# ----------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr) 
library(scales) 

cat("\n--- Generating Proteogenomics Application Results ---\n")

# Ensure the figures directory exists
if (!dir.exists("application/figures")) dir.create("application/figures", recursive = TRUE)


# 1. Load the data
table_data <- readRDS("application/results/tcga_performance_table.rds")
plot_data_full <- readRDS("application/results/tcga_residuals_plotdata.rds")

# ________________________
# 2. Print Summary Table 
# ________________________

cat("\n--- Application Table Data ---\n")

# Reorder factor levels
table_data$Method <- factor(table_data$Method, levels = c("OLS", "Standard MM", "Deterministic MM", "ROBU"))
table_data <- table_data[order(table_data$Method), ]

# Calculate population sizes to compute rates
n_total <- 882
cont_prop <- 0.15
n_outliers <- floor(n_total * cont_prop)  # 132
n_clean <- n_total - n_outliers           # 750

# Format numbers and convert TP/FP to percentages
table_data <- table_data |>
  mutate(
    Time = sprintf("%.1f", Time_Seconds),
    MSE = sprintf("%.1f", MSE_vs_Baseline),
    `TP Rate` = sprintf("%.1f%%", (TP / n_outliers) * 100),
    `FP Rate` = sprintf("%.1f%%", (FP / n_clean) * 100)
  )

# Fix NAs for OLS
table_data$`TP Rate`[table_data$Method == "OLS"] <- "--"
table_data$`FP Rate`[table_data$Method == "OLS"] <- "--"

# Order columns exactly to match the manuscript: Method, MSE, TP Rate, FP Rate, Time
table_data <- table_data |> select(Method, MSE, `TP Rate`, `FP Rate`, Time)

# Print cleanly to console
print(as.data.frame(table_data), row.names = FALSE)

# ____________________________________________
# 3. Generate Two-Panel Residuals Scatterplot
# ____________________________________________

cat("\n--- Generating PDF Plot ---\n")

# Extract only the 50th replication to keep the plot clean (882 points)
plot_data <- plot_data_full |> filter(Rep == 50)

# A. Format the data for plotting
plot_data$Observation <- ifelse(plot_data$Is_Contaminated, "Adversarial Outlier", "Clean Observation")
plot_data$Observation <- factor(plot_data$Observation, levels = c("Clean Observation", "Adversarial Outlier"))

# Reshape the data using RAW Residuals!
plot_long <- pivot_longer(
  data = plot_data,
  cols = c(Res_MM_Contam, Res_ROBU_Contam), 
  names_to = "Method",
  values_to = "Contaminated_Residual"
)

# Clean up the method names for the plot titles
plot_long$Method <- ifelse(plot_long$Method == "Res_MM_Contam", 
                           "Standard MM-Estimator", 
                           "ROBU Algorithm")

plot_long$Method <- factor(plot_long$Method, levels = c("Standard MM-Estimator", "ROBU Algorithm"))

# B. Build the Y=X plot WITHOUT squish!
p <- ggplot(plot_long, aes(x = Res_CleanBaseline, y = Contaminated_Residual, 
                           color = Observation, shape = Observation)) +
  
  # Identity line (Perfect recovery)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = c("Clean Observation" = "#2C7BB6", "Adversarial Outlier" = "#D7191C")) +
  scale_shape_manual(values = c("Clean Observation" = 16, "Adversarial Outlier" = 4)) +
  facet_wrap(~ Method, scales = "fixed") + 
  
  # Zoom the plot to [-5, 5]. By removing "oob = scales::squish", 
  # the massive outliers will naturally fall off the canvas!
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(-5, 5)) +
  
  labs(
    x = "Baseline Raw Residuals (Clean Data)",
    y = "Estimated Raw Residuals (Contaminated Data)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "#f0f0f0"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  )

# Save the plot (suppressing warnings about removed outlier points off the canvas)
plot_file <- "application/figures/tcga_residuals_plot.pdf"
suppressWarnings(ggsave(filename = plot_file, plot = p, width = 10, height = 5.5, dpi = 300))

cat(sprintf("Plot successfully saved to: %s\n", plot_file))