# ------------------------------------------------------------------
# ROBU: Generate Plots and Tables for Main Simulations (Section 4)
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # Replaces gridExtra
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
# 2. Aggregate and Format the Data
# _______________________________________________________

# Aggregate means and format strings globally (No TP/FP/Conv needed for tables)
agg_data <- results |>
  group_by(Scenario, P, Method, Contamination) |>
  summarize(
    MSE = mean(MSE, na.rm = TRUE),
    Time = mean(Time, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    MSE_fmt = ifelse(MSE > 1000, formatC(round(MSE), format = "f", big.mark = ",", digits = 0), sprintf("%.2f", MSE)),
    Time_fmt = sprintf("%.2f", Time)
  )

# _______________________________________________________
# 3. Build Table 1 (Main Text: Clean + Leverage Points)
# _______________________________________________________

cat("\n--- Table 1 Data (Clean Data + Leverage Points) ---\n")

# A. Extract Clean Data (eps = 0)
clean_df <- agg_data |> 
  filter(Scenario == "Clean", Contamination == 0) |>
  select(P, Method, MSE_0 = MSE_fmt, Time_0 = Time_fmt)

# B. Extract Leverage Points and pivot wider
lev_df <- agg_data |> 
  filter(Scenario == "Leverage Points") |>
  select(P, Method, Contamination, MSE = MSE_fmt, Time = Time_fmt) |>
  pivot_wider(
    names_from = Contamination, 
    values_from = c(MSE, Time),
    names_glue = "{.value}_{Contamination}"
  ) |>
  select(P, Method, 
         MSE_0.1, Time_0.1,
         MSE_0.2, Time_0.2,
         MSE_0.3, Time_0.3)

# C. Join them together
table1_main <- left_join(clean_df, lev_df, by = c("P", "Method")) |> arrange(P, Method)

# Print cleanly to console
print(as.data.frame(table1_main), row.names = FALSE)


# _______________________________________________________
# 4. Build Table S1 (Supplement: Vertical Outliers)
# _______________________________________________________

cat("\n--- Table S1 Data (Vertical Outliers) ---\n")

table_supp <- agg_data |> 
  filter(Scenario == "Vertical Outliers") |>
  select(P, Method, Contamination, MSE = MSE_fmt, Time = Time_fmt) |>
  pivot_wider(
    names_from = Contamination, 
    values_from = c(MSE, Time),
    names_glue = "{.value}_{Contamination}"
  ) |>
  select(P, Method, 
         MSE_0.1, Time_0.1,
         MSE_0.2, Time_0.2,
         MSE_0.3, Time_0.3) |>
  arrange(P, Method)

# Print cleanly to console
print(as.data.frame(table_supp), row.names = FALSE)


# _________________________________________________
# 5. Generate Two-Panel Figure (Line + Boxplots)
# _________________________________________________

cat("\n--- Generating PDF Plots ---\n")

# A. Data Prep
plot_data_box <- results |>
  filter(Scenario == "Leverage Points", Contamination == 0.20) |>
  mutate(P_factor = factor(P))

plot_data_line <- plot_data_box |>
  filter(Method != "Standard OLS") |>
  group_by(Method, P) |>
  summarize(Mean_Time = mean(Time, na.rm = TRUE), .groups = "drop")

# Colors and Shapes (Okabe-Ito Palette)
method_colors <- c("Standard OLS" = "#999999", 
                   "Standard MM" = "#D55E00", 
                   "Deterministic MM" = "#009E73", 
                   "ROBU" = "#0072B2")
method_shapes <- c("Standard OLS" = 15, "Standard MM" = 17, "Deterministic MM" = 18, "ROBU" = 16)

# Universal Theme
pub_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_blank(), # Removed titles for academic formatting
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 11),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.5, linetype = "dashed"),
    axis.line = element_line(color = "black", linewidth = 0.6)
  )

# B. Left Panel: Computation Time vs P (Line Plot)
p_time <- ggplot(plot_data_line, aes(x = P, y = Mean_Time, color = Method, shape = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3.5) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_shape_manual(values = method_shapes, drop = FALSE) +
  scale_x_continuous(breaks = c(100, 200, 400)) +
  scale_y_continuous(labels = scales::comma) + 
  labs(
    x = "Number of Predictors (p)",
    y = "Average Computation Time (Seconds)"
  ) +
  pub_theme

# C. Right Panel: MSE vs P (Boxplots)
p_mse <- ggplot(plot_data_box, aes(x = P_factor, y = MSE, fill = Method)) +
  geom_boxplot(color = "black", outlier.size = 1.5, outlier.alpha = 0.6, 
               position = position_dodge(0.8), linewidth = 0.5, alpha = 0.6) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.1, big.mark = ","),
    breaks = scales::trans_breaks("log10", function(x) 10^x, n = 5)
  ) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(
    x = "Number of Predictors (p)",
    y = "Mean Squared Error (Log Scale)"
  ) +
  pub_theme

# Combine and force ONE perfectly centered legend (No A/B tags)
combined_plot <- p_time + p_mse + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.margin = margin(t = 0)) &
  guides(
    fill = guide_legend(
      nrow = 1, 
      override.aes = list(
        shape = method_shapes,
        color = method_colors,
        fill = method_colors,
        linetype = c("blank", "solid", "solid", "solid"),
        alpha = 0.6
      )
    ),
    color = "none",
    shape = "none"
  )

plot_file <- "simulations/figures/main_simulations_plot.pdf"
if (!dir.exists("simulations/figures")) dir.create("simulations/figures", recursive = TRUE)
ggsave(filename = plot_file, plot = combined_plot, width = 11, height = 5.5, dpi = 300)

cat(sprintf("Simulation plots successfully saved to: %s\n", plot_file))