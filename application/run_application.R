# --------------------------------------------------------
# ROBU: Proteogenomics Real Data Application (TCGA BRCA)
# --------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load required libraries
library(curatedTCGAData)
library(TCGAutils)
library(robustbase)
library(RobStatTM)

# Source ROBU algorithm scripts
source("R/irwls.R")
source("R/robu.R")

# Set seed for reproducibility
set.seed(0)

# ____________________________________
# 1. Data Downloading and Processing 
# ____________________________________

cat("\n--- Preparing TCGA BRCA Dataset ---\n")

if (!dir.exists("application/data")) dir.create("application/data", recursive = TRUE)
if (!dir.exists("application/results")) dir.create("application/results", recursive = TRUE)

data_file <- "application/data/TCGA_BRCA_Matched.rds"

if(!file.exists(data_file)) {
  cat("Downloading TCGA BRCA data (may take a minute)...\n")
  
  brca_data <- curatedTCGAData::curatedTCGAData("BRCA", c("RNASeq2GeneNorm", "RPPAArray"), "2.1.1", FALSE)
  matched_data <- MultiAssayExperiment::intersectColumns(brca_data)
  
  assay_names <- names(experiments(matched_data))
  rna_name <- assay_names[grep("RNASeq2GeneNorm", assay_names)]
  prot_name <- assay_names[grep("RPPAArray", assay_names)]
  
  rna_mat <- t(assay(matched_data[[rna_name]]))
  prot_mat <- t(assay(matched_data[[prot_name]]))
  
  rownames(rna_mat) <- substr(rownames(rna_mat), 1, 12)
  rownames(prot_mat) <- substr(rownames(prot_mat), 1, 12)
  
  rna_mat <- rna_mat[!duplicated(rownames(rna_mat)), ]
  prot_mat <- prot_mat[!duplicated(rownames(prot_mat)), ]
  
  common_patients <- intersect(rownames(rna_mat), rownames(prot_mat))
  rna_mat <- rna_mat[common_patients, ]
  prot_mat <- prot_mat[common_patients, ]
  
  cat(sprintf("Matched %d patients with both RNA and Protein data.\n", length(common_patients)))
  saveRDS(list(rna = rna_mat, prot = prot_mat), data_file)
  
} else {
  cat("Loading TCGA data from local file...\n")
  mats <- readRDS(data_file)
  rna_mat <- mats$rna
  prot_mat <- mats$prot
}

# ________________________________________________________________
# 2. Define the Prediction Task (ER-alpha) & Dimension Reduction
# ________________________________________________________________

target_protein <- "ER-alpha"
cat(sprintf("\n--- Filtering Data for %s ---\n", target_protein))

y_full <- prot_mat[, target_protein]

valid_y_idx <- which(!is.na(y_full))
y_clean <- y_full[valid_y_idx]
X_full <- rna_mat[valid_y_idx, ]

gene_vars <- apply(X_full, 2, var)
top_genes <- order(gene_vars, decreasing = TRUE)[1:300]
X_clean <- scale(X_full[, top_genes])
y_clean <- scale(y_clean)

n <- nrow(X_clean)
p <- ncol(X_clean)
cat(sprintf("Final Dimensions: n = %d, p = %d\n", n, p))

# __________________________________________________
# 3. Establish the Clean Baseline ("Ground Truth")
# __________________________________________________

cat("\n--- Establishing Clean Baseline ---\n")

# Use the massively expanded "fair" budget
my.control <- lmrob.control(
  method = "MM", 
  fast.s.large.n = Inf, 
  k.max = 5000,           
  nResample = 5000,
  k.m_s = 2000,
  max.it = 1000,
  refine.tol = 1e-5      
)

cat("Fitting Standard MM to clean data (this may take a few minutes)...\n")
fit_clean <- robustbase::lmrob(y_clean ~ X_clean - 1, control = my.control)
beta_clean <- coef(fit_clean)
scale_clean <- fit_clean$scale

# ______________________________________________________________
# 4. Repeated Contamination & Benchmarking (50 Replications)
# ______________________________________________________________

# Calculate optimal block size (10 variables per block)
k.blocks <- max(1, floor(p / 10))
n_reps <- 50
cont_prop <- 0.15
n_cont <- floor(cont_prop * n)

results_list <- list()
plot_data_list <- list() # List to store all 50 reps of residuals

cat(sprintf("\n--- Running %d Replications of Adversarial Contamination ---\n", n_reps))

# Helper to calculate MSE relative to clean baseline (Sum of Squared Errors)
calc_mse <- function(beta_est) {
  sum((beta_est - beta_clean)^2)
}

# Helper to calculate TP/FP
calc_tp_fp <- function(fit, is_outlier) {
  wt <- fit$rweights
  if (is.null(wt)) wt <- fit$weights
  if (is.null(wt)) return(c(NA, NA))
  
  # Weights strictly 0 indicate trimmed/rejected outliers
  is_detected <- wt == 0
  tp <- sum(is_detected & is_outlier)
  fp <- sum(is_detected & !is_outlier)
  return(c(tp, fp))
}

# Helper to silently catch warnings to prevent console spam
evaluate_method <- function(method_name, expr, is_outlier) {
  s_conv <- TRUE
  m_conv <- TRUE
  start_time <- proc.time()["elapsed"]
  
  fit <- tryCatch({
    withCallingHandlers(
      expr,
      warning = function(w) {
        if (grepl("S refinements did not converge", w$message) || grepl("S-step", w$message)) {
          s_conv <<- FALSE; invokeRestart("muffleWarning") 
        }
        if (grepl("not converged", w$message) && !grepl("S refinements", w$message)) {
          m_conv <<- FALSE; invokeRestart("muffleWarning") 
        }
      }
    )
  }, error = function(e) { s_conv <<- FALSE; m_conv <<- FALSE; return(NULL) })
  
  time_elapsed <- proc.time()["elapsed"] - start_time
  
  if (is.null(fit) || any(is.na(fit$coefficients))) {
    return(list(Method = method_name, MSE = NA, Time = time_elapsed, TP = NA, FP = NA, 
                S_Conv = s_conv, M_Conv = m_conv, fit = NULL))
  }
  
  mse <- calc_mse(coef(fit))
  
  if (method_name == "OLS") {
    tp <- NA; fp <- NA
  } else {
    stats <- calc_tp_fp(fit, is_outlier)
    tp <- stats[1]; fp <- stats[2]
  }
  
  return(list(Method = method_name, MSE = mse, Time = unname(time_elapsed), 
              TP = tp, FP = fp, S_Conv = s_conv, M_Conv = m_conv, fit = fit))
}

for (rep in 1:n_reps) {
  cat(sprintf("Running Rep %d / %d...\n", rep, n_reps))
  
  # A. Inject Sneaky Adversarial Leverage Points N(10, 1)
  X_cont <- X_clean
  y_cont <- y_clean
  
  bad_patients <- sample(1:n, n_cont)
  bad_genes <- sample(1:p, 30)
  
  X_cont[bad_patients, bad_genes] <- X_cont[bad_patients, bad_genes] + 
    matrix(rnorm(n_cont * 30, mean = 10, sd = 1), nrow = n_cont, ncol = 30)
  
  beta_bad <- rnorm(p, mean = -5, sd = 2)
  y_cont[bad_patients] <- X_cont[bad_patients, ] %*% beta_bad + rnorm(n_cont, 0, 1)
  
  is_outlier <- rep(FALSE, n)
  is_outlier[bad_patients] <- TRUE
  
  # B. Evaluate Methods
  res_ols  <- evaluate_method("OLS", lm(y_cont ~ X_cont - 1), is_outlier)
  res_mm   <- evaluate_method("Standard MM", robustbase::lmrob(y_cont ~ X_cont - 1, control = my.control), is_outlier)
  res_det  <- evaluate_method("Deterministic MM", RobStatTM::lmrobdetMM(y_cont ~ X_cont - 1), is_outlier)
  res_robu <- evaluate_method("ROBU", robu(x = X_cont, y = y_cont, k = k.blocks, robu.control = my.control, m.control = my.control), is_outlier)
  
  # C. Store Metrics
  results_list[[rep]] <- data.frame(
    Rep = rep,
    Method = c("OLS", "Standard MM", "Deterministic MM", "ROBU"),
    Time_Seconds = c(res_ols$Time, res_mm$Time, res_det$Time, res_robu$Time),
    MSE_vs_Baseline = c(res_ols$MSE, res_mm$MSE, res_det$MSE, res_robu$MSE),
    TP = c(res_ols$TP, res_mm$TP, res_det$TP, res_robu$TP),
    FP = c(res_ols$FP, res_mm$FP, res_det$FP, res_robu$FP)
  )
  
  # D. Store Raw AND Standardized Residuals for all 50 Reps
  res_baseline <- drop(y_clean - (X_clean %*% beta_clean))
  res_mm_raw   <- if(is.null(res_mm$fit)) rep(0, n) else drop(y_cont - (X_cont %*% coef(res_mm$fit)))
  res_robu_raw <- if(is.null(res_robu$fit)) rep(0, n) else drop(y_cont - (X_cont %*% coef(res_robu$fit)))
  
  std_baseline <- res_baseline / scale_clean
  std_mm_std   <- if(is.null(res_mm$fit)) rep(0, n) else res_mm_raw / res_mm$fit$scale
  std_robu_std <- if(is.null(res_robu$fit)) rep(0, n) else res_robu_raw / res_robu$fit$scale
  
  plot_data_list[[rep]] <- data.frame(
    Rep = rep,
    Patient_ID = rownames(X_cont),
    Is_Contaminated = is_outlier,
    Res_CleanBaseline = res_baseline,
    Res_MM_Contam = res_mm_raw,
    Res_ROBU_Contam = res_robu_raw,
    StdRes_CleanBaseline = std_baseline,
    StdRes_MM_Contam = std_mm_std,
    StdRes_ROBU_Contam = std_robu_std
  )
}

# ______________________________
# 5. Summarize and Save Results
# ______________________________

all_results <- do.call(rbind, results_list)
saveRDS(all_results, "application/results/tcga_full_results.rds")

summary_table <- aggregate(cbind(Time_Seconds, MSE_vs_Baseline, TP, FP) ~ Method, data = all_results, FUN = mean, na.rm = TRUE)

cat("\n======================================================\n")
cat("   PROTEOGENOMICS APPLICATION: AVERAGE OVER 50 REPS   \n")
cat("======================================================\n")
print(summary_table, row.names = FALSE)
cat("======================================================\n")

saveRDS(summary_table, "application/results/tcga_performance_table.rds")

# Bind all 50 reps of residuals and save
plot_data <- do.call(rbind, plot_data_list)
saveRDS(plot_data, "application/results/tcga_residuals_plotdata.rds")

cat("\nAnalysis complete. Results and plot data saved to application/results/\n")