# ----------------------------------------------------------------
# ROBU: Generate Plots and Tables for Proteogenomics Application
# ----------------------------------------------------------------

library(ggplot2)
library(xtable)
library(tidyr) 
library(scales) 

cat("\n--- Generating Proteogenomics Application Results ---\n")

# Ensure the figures directory exists
if (!dir.exists("application/figures")) dir.create("application/figures", recursive = TRUE)

# 1. Load the data
table_data <- readRDS("application/results/tcga_performance_table.rds")
plot_data <- readRDS("application/results/tcga_residuals_plotdata.rds")

# ________________________
# 2. Generate LaTeX Table 
# ________________________

cat("\n--- LaTeX Table Code ---\n")

# Reorder factor levels so Baseline OLS is first, then Standard MM, Det MM, then ROBU
table_data$Method <- factor(table_data$Method, levels = c("OLS", "Standard MM", "Deterministic MM", "ROBU"))
table_data <- table_data[order(table_data$Method), ]

# Format numbers
table_data$Time <- sprintf("%.2f", table_data$Time)
table_data$MSE <- sprintf("%.2f", table_data$MSE) 
table_data$TP <- sprintf("%.1f", table_data$TP)
table_data$FP <- sprintf("%.1f", table_data$FP)

# Fix NAs for OLS
table_data$TP[table_data$Method == "OLS"] <- "--"
table_data$FP[table_data$Method == "OLS"] <- "--"

colnames(table_data) <- c("Method", "Time (s)", "MSE (vs. Baseline)", "TP", "FP")

latex_table <- xtable(
  table_data, 
  caption = "Performance of estimators on the TCGA breast cancer proteogenomics dataset. Results are averaged over 50 replications with 15\\% artificial adversarial contamination.",
  label = "tab:real_data",
  align = c("l", "l", "c", "c", "c", "c")
)

print(latex_table, include.rownames = FALSE, caption.placement = "top", booktabs = TRUE)

# ____________________________________________
# 3. Generate Two-Panel Residuals Scatterplot
# ____________________________________________

cat("\n--- Generating PDF Plot ---\n")

# A. Format the data for plotting
plot_data$Observation <- ifelse(plot_data$Is_Contaminated, "Adversarial Outlier", "Clean Observation")
plot_data$Observation <- factor(plot_data$Observation, levels = c("Clean Observation", "Adversarial Outlier"))

# Reshape the data from "wide" to "long" format so we can use facet_wrap
plot_long <- pivot_longer(
  data = plot_data,
  cols = c(StdRes_MM_Contam, StdRes_ROBU_Contam), 
  names_to = "Method",
  values_to = "Contaminated_Residual"
)

# Clean up the method names for the plot titles
plot_long$Method <- ifelse(plot_long$Method == "StdRes_MM_Contam", 
                           "Standard MM-Estimator", 
                           "ROBU Algorithm")

# Set the factor levels so Standard MM is on the left, ROBU is on the right
plot_long$Method <- factor(plot_long$Method, 
                           levels = c("Standard MM-Estimator", "ROBU Algorithm"))

# Tukey Bisquare rejection boundary (approx 4.685 for 95% efficiency)
tukey_c <- 4.685

# B. Build the plot
p <- ggplot(plot_long, aes(x = StdRes_CleanBaseline, y = Contaminated_Residual, 
                           color = Observation, shape = Observation)) +
  
  # Identity line (Perfect recovery)
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", linewidth = 0.5) +
  
  # Tukey rejection boundaries
  geom_hline(yintercept = tukey_c, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = -tukey_c, linetype = "dashed", color = "red", alpha = 0.5) +
  
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = c("Clean Observation" = "#2C7BB6", "Adversarial Outlier" = "#D7191C")) +
  scale_shape_manual(values = c("Clean Observation" = 16, "Adversarial Outlier" = 4)) +
  facet_wrap(~ Method, scales = "fixed") + 
  
  # Zoom the plot to [-6, 6] to show the Tukey boundary clearly
  scale_x_continuous(limits = c(-6, 6), oob = scales::squish) +
  scale_y_continuous(limits = c(-6, 6), oob = scales::squish) +
  
  labs(
    x = "Baseline Standardized Residuals (Clean Data)",
    y = "Estimated Standardized Residuals (Contaminated Data)",
    title = NULL
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

# Save the plot 
plot_file <- "application/figures/tcga_residuals_plot.pdf"
ggsave(filename = plot_file, plot = p, width = 10, height = 5.5, dpi = 300)

cat(sprintf("Plot successfully saved to: %s\n", plot_file))