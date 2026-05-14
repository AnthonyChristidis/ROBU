# --------------------------------------------------------
# ROBU: Proteogenomics Real Data Application (TCGA BRCA)
# --------------------------------------------------------

rm(list = ls())

library(curatedTCGAData)
library(TCGAutils)
library(robustbase)
library(RobStatTM)

source("R/irwls.R")
source("R/robu.R")

set.seed(0)

# ____________________________________
# 1. Data Downloading and Processing 
# ____________________________________
cat("\n--- Preparing TCGA BRCA Dataset ---\n")

if (!dir.exists("application/data")) dir.create("application/data", recursive = TRUE)
if (!dir.exists("application/results")) dir.create("application/results", recursive = TRUE)

data_file <- "application/data/TCGA_BRCA_Matched.rds"

if(!file.exists(data_file)) {
  cat("Downloading TCGA BRCA data...\n")
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
  
  saveRDS(list(rna = rna_mat, prot = prot_mat), data_file)
} else {
  mats <- readRDS(data_file)
  rna_mat <- mats$rna
  prot_mat <- mats$prot
}

# ________________________________________________________________
# 2. Define the Prediction Task (ER-alpha) & Dimension Reduction
# ________________________________________________________________
target_protein <- "ER-alpha"
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

# __________________________________________________
# 3. Establish the Clean Baseline ("Ground Truth")
# __________________________________________________
my.control <- lmrob.control(
  method = "MM", fast.s.large.n = Inf, 
  k.max = 5000, nResample = 5000, k.m_s = 2000, max.it = 1000, refine.tol = 1e-5      
)

fit_clean <- robustbase::lmrob(y_clean ~ X_clean - 1, control = my.control)
beta_clean <- coef(fit_clean)
scale_clean <- fit_clean$scale

# ______________________________________________________________
# 4. Repeated Contamination & Benchmarking (50 Replications)
# ______________________________________________________________
k.blocks <- max(1, floor(p / 20))
n_reps <- 50
cont_prop <- 0.15
n_cont <- floor(cont_prop * n)

results_list <- list()
plot_data <- NULL 

# Helper function (Includes S_Conv and M_Conv tracking!)
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
  
  mse <- sum((coef(fit) - beta_clean)^2)
  
  if (method_name == "OLS") {
    tp <- NA; fp <- NA
  } else {
    wt <- fit$rweights
    if (is.null(wt)) wt <- fit$weights
    tp <- sum(wt < 0.5 & is_outlier)
    fp <- sum(wt < 0.5 & !is_outlier)
  }
  
  return(list(Method = method_name, MSE = mse, Time = unname(time_elapsed), 
              TP = tp, FP = fp, S_Conv = s_conv, M_Conv = m_conv, fit = fit))
}

for (rep in 1:n_reps) {
  cat(sprintf("Running Rep %d / %d...\n", rep, n_reps))
  
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
  
  res_ols <- evaluate_method("OLS", lm(y_cont ~ X_cont - 1), is_outlier)
  res_mm  <- evaluate_method("Standard MM", robustbase::lmrob(y_cont ~ X_cont - 1, control = my.control), is_outlier)
  res_det <- evaluate_method("Deterministic MM", RobStatTM::lmrobdetMM(y_cont ~ X_cont - 1), is_outlier)
  res_robu<- evaluate_method("ROBU", robu(x = X_cont, y = y_cont, k = k.blocks, robu.control = my.control, m.control = my.control), is_outlier)
  
  results_list[[rep]] <- data.frame(
    Rep = rep, Method = c("OLS", "Standard MM", "Deterministic MM", "ROBU"),
    Time = c(res_ols$Time, res_mm$Time, res_det$Time, res_robu$Time),
    MSE = c(res_ols$MSE, res_mm$MSE, res_det$MSE, res_robu$MSE),
    TP = c(res_ols$TP, res_mm$TP, res_det$TP, res_robu$TP),
    FP = c(res_ols$FP, res_mm$FP, res_det$FP, res_robu$FP),
    S_Conv = c(res_ols$S_Conv, res_mm$S_Conv, res_det$S_Conv, res_robu$S_Conv),
    M_Conv = c(res_ols$M_Conv, res_mm$M_Conv, res_det$M_Conv, res_robu$M_Conv)
  )
  
  if (rep == n_reps) {
    # Extract standardized residuals (r_i / scale)
    std_res_baseline <- drop(y_clean - (X_clean %*% beta_clean)) / scale_clean
    
    std_res_mm <- if(is.null(res_mm$fit)) rep(0, n) else drop(y_cont - (X_cont %*% coef(res_mm$fit))) / res_mm$fit$scale
    std_res_robu <- if(is.null(res_robu$fit)) rep(0, n) else drop(y_cont - (X_cont %*% coef(res_robu$fit))) / res_robu$fit$scale
    
    plot_data <- data.frame(
      Patient_ID = rownames(X_cont),
      Is_Contaminated = is_outlier,
      StdRes_CleanBaseline = std_res_baseline,
      StdRes_MM_Contam = std_res_mm,
      StdRes_ROBU_Contam = std_res_robu
    )
  }
}

all_results <- do.call(rbind, results_list)
saveRDS(all_results, "application/results/tcga_full_results.rds")

summary_table <- aggregate(cbind(Time, MSE, TP, FP) ~ Method, data = all_results, FUN = mean, na.rm = TRUE)
saveRDS(summary_table, "application/results/tcga_performance_table.rds")
saveRDS(plot_data, "application/results/tcga_residuals_plotdata.rds")