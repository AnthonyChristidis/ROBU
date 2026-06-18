# -------------------------------------------------------------------------
# Supplementary Study: ROBU Blockwise Resampling Budget (Visualization)
# -------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

cat("\n--- Generating Budget Supplement Results ---\n")

# Ensure figures directory exists
dir.create("simulations/figures", recursive = TRUE, showWarnings = FALSE)

# Load the raw results
raw_results_file <- "simulations/results/robu_budget_supplement_raw.rds"

if (!file.exists(raw_results_file)) {
  stop("Results file not found. Please run simulations/run_budget_study.R first.")
}

results <- readRDS(raw_results_file)

# Format epsilon for nicer x-axis labels in the plot
results$Eps_Factor <- factor(results$Eps, levels = c(0.20, 0.30), 
                             labels = c("0.20", "0.30"))

# Set a colorblind-friendly palette with TWO SHADES OF BLUE for ROBU
cb_palette <- c("ROBU (Full Budget)" = "#0072B2",         # Dark/Standard Blue
                "ROBU (Distributed Budget)" = "#56B4E9")  # Light Blue

# Universal Theme (matching your main simulation script)
pub_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_blank(), 
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 11),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.5, linetype = "dashed"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    legend.margin = margin(t = 0)
  )

# __________________________
# Plot 1: MSE Boxplot Only
# __________________________

cat("Generating MSE Boxplot...\n")

p_mse <- ggplot(results, aes(x = Eps_Factor, y = MSE, fill = Method)) +
  geom_boxplot(color = "black", outlier.size = 1.5, outlier.alpha = 0.6, 
               width = 0.45, position = position_dodge(width = 0.65), 
               linewidth = 0.5, alpha = 0.8) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1, big.mark = ","),
    breaks = scales::trans_breaks("log10", function(x) 10^x, n = 5)
  ) +
  scale_fill_manual(values = cb_palette) +
  labs(x = expression("Contamination Level (" * epsilon * ")"),
       y = "Mean Squared Error (Log Scale)") +
  pub_theme

# Save to PDF
plot_file <- "simulations/figures/robu_budget_supplement_mse.pdf"
ggsave(plot_file, plot = p_mse, width = 6.5, height = 5.5, dpi = 300)
cat(sprintf("Saved plot to: %s\n", plot_file))

# ____________________________________________________
# Table: Computation Time and MSE (For LaTeX writing)
# ____________________________________________________

cat("\n--- Table S2 Data (Budget Supplement: MSE & Time) ---\n")

table_supp <- results |>
  group_by(Method, Eps) |>
  summarize(
    MSE = mean(MSE, na.rm = TRUE),
    Time = mean(Time, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    MSE_fmt = ifelse(MSE > 1000, 
                     formatC(round(MSE, 1), format = "f", big.mark = ",", digits = 1), 
                     sprintf("%.1f", MSE)),
    Time_fmt = sprintf("%.1f", Time)
  ) |>
  select(Method, Eps, MSE = MSE_fmt, Time = Time_fmt) |>
  pivot_wider(
    names_from = Eps, 
    values_from = c(MSE, Time),
    names_glue = "{.value}_{Eps}"
  ) |>
  select(Method, 
         MSE_0.2, Time_0.2,
         MSE_0.3, Time_0.3)

# Print cleanly to console so you can copy-paste into LaTeX
print(as.data.frame(table_supp), row.names = FALSE)
cat("\nDone!\n")