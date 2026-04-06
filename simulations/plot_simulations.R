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
results$Method <- factor(results$Method, levels = c("Standard OLS", "Standard MM", "ROBU"))
results$Scenario <- factor(results$Scenario, levels = c("Clean", "Vertical Outliers", "Leverage Points"))

# Filter to ONLY 0% and 15% contamination for the main manuscript
results_main <- results |> filter(Contamination %in% c(0, 0.15))

# _______________________________________________________
# 2. Generate Master LaTeX Table (MSE Only, 15% Contam)
# _______________________________________________________

cat("\n--- Generating LaTeX Table ---\n")

# Aggregate results (Mean MSE)
summary_table <- results_main |>
  group_by(Scenario, Method, P) |>
  summarize(
    Mean_MSE = mean(MSE, na.rm = TRUE),
    .groups = "drop"
  )

# Format the metrics: Drop decimals and add commas if > 100, else keep 2 decimals
summary_table <- summary_table |>
  mutate(
    Formatted = ifelse(Mean_MSE > 100, 
                       formatC(round(Mean_MSE), format = "f", big.mark = ",", digits = 0), 
                       sprintf("%.2f", Mean_MSE))
  ) |>
  select(Scenario, Method, P, Formatted)

# Pivot wider so that P = 100, 200, 400 become columns
table_wide <- summary_table |>
  pivot_wider(names_from = P, values_from = Formatted, names_prefix = "p_") |>
  arrange(Scenario, Method)

# Clean up Scenario labels for the LaTeX table
table_wide$Scenario <- as.character(table_wide$Scenario)
table_wide$Scenario[table_wide$Scenario == "Clean"] <- "\\textbf{Clean Data} (0\\% Contam.)"
table_wide$Scenario[table_wide$Scenario == "Vertical Outliers"] <- "\\textbf{Vertical Outliers} (15\\% Contam.)"
table_wide$Scenario[table_wide$Scenario == "Leverage Points"] <- "\\textbf{Leverage Points} (15\\% Contam.)"

# Rename columns for xtable
colnames(table_wide) <- c("\\textbf{Scenario}", "\\textbf{Method}", 
                          "\\textbf{\\boldmath $p=100$}", "\\textbf{\\boldmath $p=200$}", "\\textbf{\\boldmath $p=400$}")

latex_xtable <- xtable(
  table_wide, 
  caption = "Average Mean Squared Error (MSE) across 100 replications. Results highlight the empirical breakdown of the Standard MM-estimator at $p=400$, and the successful signal recovery by ROBU under 15\\% contamination.",
  label = "tab:sim_mse",
  align = c("l", "l", "l", "r", "r", "r") # Left align text, Right align numbers
)

print(latex_xtable, include.rownames = FALSE, sanitize.text.function = identity, 
      sanitize.colnames.function = identity, caption.placement = "top", booktabs = TRUE)


# _________________________________________________
# 3. Generate Two-Panel Figure (Time and MSE vs p)
# _________________________________________________

cat("\n--- Generating PDF Plots ---\n")

# We will focus the plot on the most critical scenario: 15% Bad Leverage Points
# This perfectly illustrates the empirical breakdown transition.
plot_data <- results_main |>
  filter(Scenario == "Leverage Points") |>
  group_by(Method, P) |>
  summarize(
    Mean_Time = mean(Time, na.rm = TRUE),
    Mean_MSE = mean(MSE, na.rm = TRUE),
    .groups = "drop"
  )

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

# Custom color palette
method_colors <- c("Standard OLS" = "#969696", "Standard MM" = "#D7191C", "ROBU" = "#2C7BB6")
method_shapes <- c("Standard OLS" = 15, "Standard MM" = 17, "ROBU" = 16)

# A. Left Panel: Computation Time vs P
p_time <- ggplot(plot_data |> filter(Method != "Standard OLS"), # OLS is too fast/irrelevant for time comparison
                 aes(x = P, y = Mean_Time, color = Method, shape = Method, group = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3.5) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_x_continuous(breaks = c(100, 200, 400)) +
  labs(
    x = "Number of Predictors (p)",
    y = "Average Computation Time (Seconds)",
    title = "Computational Scaling"
  ) +
  pub_theme

# B. Right Panel: MSE vs P
p_mse <- ggplot(plot_data, aes(x = P, y = Mean_MSE, color = Method, shape = Method, group = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3.5) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1, big.mark = ","),
    breaks = c(0.1, 1, 10, 100, 1000)
  ) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_x_continuous(breaks = c(100, 200, 400)) +
  labs(
    x = "Number of Predictors (p)",
    y = "Mean Squared Error (Log Scale)",
    title = "Empirical Breakdown Resistance"
  ) +
  pub_theme

# Combine the two plots side-by-side
combined_plot <- grid.arrange(p_time, p_mse, ncol = 2)

# Save to PDF
plot_file <- "simulations/figures/main_simulations_plot.pdf"
if (!dir.exists("simulations/figures")) dir.create("simulations/figures", recursive = TRUE)
ggsave(filename = plot_file, plot = combined_plot, width = 11, height = 5.5, dpi = 300)

cat(sprintf("Simulation plots successfully saved to: %s\n", plot_file))