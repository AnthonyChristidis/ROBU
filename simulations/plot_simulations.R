# ------------------------------------------------------------------
# ROBU: Generate Plots and Tables for Main Simulations (Section 4)
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)
library(gridExtra)
library(scales)

cat("\n--- Generating Main Simulation Results ---\n")

# 1. Load all simulation .rds files
result_files <- list.files("simulations/results", pattern = "sim_n1000_p.*\\.rds$", full.names = TRUE)
if (length(result_files) == 0) stop("No simulation files found! Check the directory.")

results <- bind_rows(lapply(result_files, readRDS)) |> distinct()

# Rename OLS to Standard OLS to match the manuscript text
results$Method <- ifelse(results$Method == "OLS", "Standard OLS", as.character(results$Method))

# Ensure factors are ordered logically
results$Method <- factor(results$Method, levels = c("Standard OLS", "Standard MM", "Deterministic MM", "ROBU"))
results$Scenario <- factor(results$Scenario, levels = c("Clean", "Vertical Outliers", "Leverage Points"))

# _______________________________________________________
# 2. Generate Master LaTeX Table (Leverage Points Only)
# _______________________________________________________

cat("\n--- Generating LaTeX Table ---\n")

# Filter to only Leverage Points for the main text table
table_data <- results |> 
  filter(Scenario == "Leverage Points") |>
  group_by(P, Method, Contamination) |>
  summarize(
    MSE = mean(MSE, na.rm = TRUE),
    Time = mean(Time, na.rm = TRUE),
    TP = mean(TP, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    .groups = "drop"
  )

# Format the metrics
table_data <- table_data |>
  mutate(
    MSE = ifelse(MSE > 1000, formatC(round(MSE), format = "f", big.mark = ",", digits = 0), sprintf("%.2f", MSE)),
    Time = sprintf("%.2f", Time),
    TP = ifelse(is.na(TP), "--", sprintf("%.1f", TP)),
    FP = ifelse(is.na(FP), "--", sprintf("%.1f", FP))
  )

# Pivot wider so that Contamination = 0.10, 0.20, 0.30 become columns
table_wide <- table_data |>
  pivot_wider(
    names_from = Contamination, 
    values_from = c(MSE, TP, FP, Time),
    names_glue = "{.value}_{Contamination}"
  ) |>
  arrange(P, Method)

# Reorder columns to match: P, Method, (MSE, TP, FP, Time) for 0.10, 0.20, 0.30
table_wide <- table_wide |>
  select(P, Method, 
         MSE_0.1, TP_0.1, FP_0.1, Time_0.1,
         MSE_0.2, TP_0.2, FP_0.2, Time_0.2,
         MSE_0.3, TP_0.3, FP_0.3, Time_0.3)

# Rename columns for xtable
colnames(table_wide) <- c("$\\boldsymbol{p}$", "\\textbf{Method}", 
                          "\\textbf{MSE}", "\\textbf{TP}", "\\textbf{FP}", "\\textbf{Time}",
                          "\\textbf{MSE}", "\\textbf{TP}", "\\textbf{FP}", "\\textbf{Time}",
                          "\\textbf{MSE}", "\\textbf{TP}", "\\textbf{FP}", "\\textbf{Time}")

latex_xtable <- xtable(
  table_wide, 
  caption = "Average performance metrics under Scenario 3 (Adversarial Leverage Points) across 50 replications. True Positives (TP) and False Positives (FP) indicate the average number of observations trimmed by the robust estimators.",
  label = "tab:sim_leverage",
  align = c("l", "l", "l",  "c","c","c","c",  "c","c","c","c",  "c","c","c","c")
)

print(latex_xtable, include.rownames = FALSE, sanitize.text.function = identity, 
      sanitize.colnames.function = identity, caption.placement = "top", booktabs = TRUE)


# _________________________________________________
# 3. Generate Two-Panel Figure (Boxplots)
# _________________________________________________

cat("\n--- Generating PDF Plots ---\n")

# Filter raw data for boxplots (Fixing epsilon at 0.20 to show scaling behavior)
plot_data <- results |>
  filter(Scenario == "Leverage Points", Contamination == 0.20) |>
  mutate(P_factor = factor(P)) # Convert P to factor so boxplots group correctly by x-axis

# Professional Academic Theme
pub_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 11),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

# Custom color palette for boxplot fills
method_fills <- c("Standard OLS" = "#f0f0f0", "Standard MM" = "#Fcae91", 
                  "Deterministic MM" = "#bae4b3", "ROBU" = "#bdd7e7")
method_colors <- c("Standard OLS" = "#969696", "Standard MM" = "#D7191C", 
                   "Deterministic MM" = "#31a354", "ROBU" = "#2C7BB6")

# A. Left Panel: Computation Time vs P (Boxplots)
p_time <- ggplot(plot_data |> filter(Method != "Standard OLS"), 
                 aes(x = P_factor, y = Time, fill = Method, color = Method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5, position = position_dodge(0.8)) +
  scale_fill_manual(values = method_fills) +
  scale_color_manual(values = method_colors) +
  labs(
    x = "Number of Predictors (p)",
    y = "Computation Time (Seconds)",
    title = "Computational Scaling"
  ) +
  pub_theme

# B. Right Panel: MSE vs P (Boxplots)
p_mse <- ggplot(plot_data, aes(x = P_factor, y = MSE, fill = Method, color = Method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5, position = position_dodge(0.8)) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1, big.mark = ","),
    breaks = scales::trans_breaks("log10", function(x) 10^x, n = 5)
  ) +
  scale_fill_manual(values = method_fills) +
  scale_color_manual(values = method_colors) +
  labs(
    x = "Number of Predictors (p)",
    y = "Mean Squared Error (Log Scale)",
    title = "Algorithmic Failure Resistance"
  ) +
  pub_theme

# Combine the two plots side-by-side
combined_plot <- grid.arrange(p_time, p_mse, ncol = 2)

# Save to PDF
plot_file <- "simulations/figures/main_simulations_plot.pdf"
if (!dir.exists("simulations/figures")) dir.create("simulations/figures", recursive = TRUE)
ggsave(filename = plot_file, plot = combined_plot, width = 11, height = 5.5, dpi = 300)

cat(sprintf("Simulation plots successfully saved to: %s\n", plot_file))