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

# Rename OLS to Standard OLS
results$Method <- ifelse(results$Method == "OLS", "Standard OLS", as.character(results$Method))

# Ensure factors are ordered logically
results$Method <- factor(results$Method, levels = c("Standard OLS", "Standard MM", "Deterministic MM", "ROBU"))
results$Scenario <- factor(results$Scenario, levels = c("Clean", "Vertical Outliers", "Leverage Points"))

# _______________________________________________________
# 2. Helper Function to Build the Wide LaTeX Tables
# _______________________________________________________

build_latex_table <- function(data, scenario_name, metric_type = "Statistical", table_caption, table_label) {
  
  # Filter data to the specific scenario
  scenario_data <- data |> filter(Scenario == scenario_name)
  
  # Aggregate means
  agg_data <- scenario_data |>
    group_by(P, Method, Contamination) |>
    summarize(
      MSE = mean(MSE, na.rm = TRUE),
      Time = mean(Time, na.rm = TRUE),
      TP = mean(TP, na.rm = TRUE),
      FP = mean(FP, na.rm = TRUE),
      S_Conv = mean(S_Conv, na.rm = TRUE) * 100,
      M_Conv = mean(M_Conv, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  # Format metrics
  agg_data <- agg_data |>
    mutate(
      MSE = ifelse(MSE > 1000, formatC(round(MSE), format = "f", big.mark = ",", digits = 0), sprintf("%.2f", MSE)),
      Time = sprintf("%.2f", Time),
      TP = ifelse(is.na(TP), "--", sprintf("%.1f", TP)),
      FP = ifelse(is.na(FP), "--", sprintf("%.1f", FP)),
      S_Conv = ifelse(Method == "Standard OLS", "--", sprintf("%.0f\\%%", S_Conv)),
      M_Conv = ifelse(Method == "Standard OLS", "--", sprintf("%.0f\\%%", M_Conv))
    )
  
  # Pivot wider based on requested metric type
  if (metric_type == "Statistical") {
    # Keep MSE, TP, FP
    cols_to_pivot <- c("MSE", "TP", "FP")
    align_str <- c("l", "l", "c",  "c","c","c",  "c","c","c",  "c","c","c")
    col_names <- c("$\\boldsymbol{p}$", "\\textbf{Method}", 
                   "\\textbf{MSE}", "\\textbf{MSE}", "\\textbf{TP}", "\\textbf{FP}",
                   "\\textbf{MSE}", "\\textbf{TP}", "\\textbf{FP}",
                   "\\textbf{MSE}", "\\textbf{TP}", "\\textbf{FP}")
  } else {
    # Keep Time, S_Conv, M_Conv
    cols_to_pivot <- c("Time", "S_Conv", "M_Conv")
    align_str <- c("l", "l", "c",  "c","c","c",  "c","c","c",  "c","c","c")
    col_names <- c("$\\boldsymbol{p}$", "\\textbf{Method}", 
                   "\\textbf{Time}", "\\textbf{Time}", "\\textbf{S-Conv}", "\\textbf{M-Conv}",
                   "\\textbf{Time}", "\\textbf{S-Conv}", "\\textbf{M-Conv}",
                   "\\textbf{Time}", "\\textbf{S-Conv}", "\\textbf{M-Conv}")
  }
  
  # Handle Clean Data (which only has Contamination == 0) vs others
  if (scenario_name == "Clean") {
    table_wide <- agg_data |>
      filter(Contamination == 0) |>
      select(P, Method, all_of(cols_to_pivot[1])) |> # Only keep MSE or Time for Clean
      arrange(P, Method)
    
    col_names <- c("$\\boldsymbol{p}$", "\\textbf{Method}", ifelse(metric_type == "Statistical", "\\textbf{MSE}", "\\textbf{Time}"))
    align_str <- c("l", "l", "c")
  } else {
    table_wide <- agg_data |>
      pivot_wider(
        names_from = Contamination, 
        values_from = all_of(cols_to_pivot),
        names_glue = "{.value}_{Contamination}"
      ) |>
      arrange(P, Method)
    
    # Reorder columns explicitly
    if (metric_type == "Statistical") {
      table_wide <- table_wide |> select(P, Method, MSE_0.1, TP_0.1, FP_0.1, MSE_0.2, TP_0.2, FP_0.2, MSE_0.3, TP_0.3, FP_0.3)
    } else {
      table_wide <- table_wide |> select(P, Method, Time_0.1, S_Conv_0.1, M_Conv_0.1, Time_0.2, S_Conv_0.2, M_Conv_0.2, Time_0.3, S_Conv_0.3, M_Conv_0.3)
    }
  }
  
  colnames(table_wide) <- col_names
  
  latex_xtable <- xtable(
    table_wide, 
    caption = table_caption,
    label = table_label,
    align = align_str
  )
  
  cat(sprintf("\n--- %s ---\n", table_label))
  print(latex_xtable, include.rownames = FALSE, sanitize.text.function = identity, 
        sanitize.colnames.function = identity, caption.placement = "top", booktabs = TRUE)
}

# _______________________________________________________
# 3. Print all 4 Tables for Copy/Pasting into LaTeX
# _______________________________________________________

# Note: In the final paper, you will manually combine the "Clean" column next to the "Leverage Points" columns.
# This script prints them perfectly formatted so you can easily copy the rows.

build_latex_table(results, "Clean", "Statistical", "Clean Data - Statistical", "tab:clean_stat")
build_latex_table(results, "Leverage Points", "Statistical", "Leverage Points - Statistical", "tab:lev_stat")

build_latex_table(results, "Clean", "Computational", "Clean Data - Computational", "tab:clean_comp")
build_latex_table(results, "Leverage Points", "Computational", "Leverage Points - Computational", "tab:lev_comp")

build_latex_table(results, "Vertical Outliers", "Statistical", "Supp: Vertical Outliers - Statistical", "tab:supp_vert_stat")
build_latex_table(results, "Vertical Outliers", "Computational", "Supp: Vertical Outliers - Computational", "tab:supp_vert_comp")


# _________________________________________________
# 4. Generate Two-Panel Figure (Boxplots)
# _________________________________________________

cat("\n--- Generating PDF Plots ---\n")

plot_data <- results |>
  filter(Scenario == "Leverage Points", Contamination == 0.20) |>
  mutate(P_factor = factor(P))

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

method_fills <- c("Standard OLS" = "#f0f0f0", "Standard MM" = "#Fcae91", 
                  "Deterministic MM" = "#bae4b3", "ROBU" = "#bdd7e7")
method_colors <- c("Standard OLS" = "#969696", "Standard MM" = "#D7191C", 
                   "Deterministic MM" = "#31a354", "ROBU" = "#2C7BB6")

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

combined_plot <- grid.arrange(p_time, p_mse, ncol = 2)

plot_file <- "simulations/figures/main_simulations_plot.pdf"
if (!dir.exists("simulations/figures")) dir.create("simulations/figures", recursive = TRUE)
ggsave(filename = plot_file, plot = combined_plot, width = 11, height = 5.5, dpi = 300)

cat(sprintf("Simulation plots successfully saved to: %s\n", plot_file))