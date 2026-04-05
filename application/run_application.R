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
  
  # Download the specific assays
  brca_data <- curatedTCGAData::curatedTCGAData(
    diseaseCode = "BRCA",
    assays = c("RNASeq2GeneNorm", "RPPAArray"),
    version = "2.1.1",
    dry.run = FALSE
  )
  
  # Match samples that have BOTH RNA and Protein data
  matched_data <- MultiAssayExperiment::intersectColumns(brca_data)
  
  assay_names <- names(experiments(matched_data))
  rna_name <- assay_names[grep("RNASeq2GeneNorm", assay_names)]
  prot_name <- assay_names[grep("RPPAArray", assay_names)]
  
  rna_mat <- t(assay(matched_data[[rna_name]]))
  prot_mat <- t(assay(matched_data[[prot_name]]))
  
  # Clean up patient IDs
  rownames(rna_mat) <- substr(rownames(rna_mat), 1, 12)
  rownames(prot_mat) <- substr(rownames(prot_mat), 1, 12)
  
  # Handle duplicates
  rna_mat <- rna_mat[!duplicated(rownames(rna_mat)), ]
  prot_mat <- prot_mat[!duplicated(rownames(prot_mat)), ]
  
  # Find exact intersection
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

# Remove NAs in the target protein
valid_y_idx <- which(!is.na(y_full))
y_clean <- y_full[valid_y_idx]
X_full <- rna_mat[valid_y_idx, ]

# Unsupervised Variance Filter: Keep top 300 most variable genes
gene_vars <- apply(X_full, 2, var)
top_genes <- order(gene_vars, decreasing = TRUE)[1:300]
X_clean <- X_full[, top_genes]

# Standardize predictors and response for numerical stability
X_clean <- scale(X_clean)
y_clean <- scale(y_clean)

n <- nrow(X_clean)
p <- ncol(X_clean)
cat(sprintf("Final Dimensions: n = %d, p = %d\n", n, p))

# __________________________________________________
# 3. Establish the Clean Baseline ("Ground Truth")
# __________________________________________________

cat("\n--- Establishing Clean Baseline ---\n")

my.control <- lmrob.control(
  method = "MM", 
  fast.s.large.n = Inf, 
  k.max = 200,           
  refine.tol = 1e-5      
)

cat("Fitting Standard MM to clean data (this may take 1-2 minutes)...\n")
fit_clean <- robustbase::lmrob(y_clean ~ X_clean - 1, control = my.control)
beta_clean <- coef(fit_clean)

# ______________________________________________________________
# 4. Repeated Contamination & Benchmarking (50 Replications)
# ______________________________________________________________

k.blocks <- max(1, floor(p / 20))
n_reps <- 50
cont_prop <- 0.15
n_cont <- floor(cont_prop * n)

results_list <- list()
plot_data <- NULL # Will store the residuals from the final rep for plotting

cat(sprintf("\n--- Running %d Replications of Adversarial Contamination ---\n", n_reps))

# Helper to calculate MSE relative to clean baseline
calc_mse <- function(beta_est) {
  mean((beta_est - beta_clean)^2)
}

# Helper to silently catch warnings to prevent console spam
run_silently <- function(expr, p_vars) {
  tryCatch({
    withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
  }, error = function(e) list(coefficients = rep(NA, p_vars), residuals = rep(NA, n)))
}

for (rep in 1:n_reps) {
  cat(sprintf("Running Rep %d / %d...\n", rep, n_reps))
  
  # A. Inject Adversarial Leverage Points
  X_cont <- X_clean
  y_cont <- y_clean
  
  bad_patients <- sample(1:n, n_cont)
  bad_genes <- sample(1:p, 30)
  
  X_cont[bad_patients, bad_genes] <- X_cont[bad_patients, bad_genes] + 
    matrix(rnorm(n_cont * 30, mean = 100, sd = 10), nrow = n_cont, ncol = 30)
  
  beta_bad <- rnorm(p, mean = -5, sd = 2)
  y_cont[bad_patients] <- X_cont[bad_patients, ] %*% beta_bad + rnorm(n_cont, 0, 1)
  
  is_outlier <- rep(FALSE, n)
  is_outlier[bad_patients] <- TRUE
  
  # B. Run OLS
  t_ols <- system.time({ fit_ols <- lm(y_cont ~ X_cont - 1) })["elapsed"]
  mse_ols <- calc_mse(coef(fit_ols))
  
  # C. Run Standard MM
  t_mm <- system.time({
    fit_mm <- run_silently(robustbase::lmrob(y_cont ~ X_cont - 1, control = my.control), p)
  })["elapsed"]
  mse_mm <- calc_mse(fit_mm$coefficients)
  
  # D. Run ROBU
  t_robu <- system.time({
    fit_robu <- run_silently(robu(x = X_cont, y = y_cont, k = k.blocks, robu.control = my.control, m.control = my.control), p)
  })["elapsed"]
  mse_robu <- calc_mse(fit_robu$coefficients)
  
  # E. Store Rep Results
  results_list[[rep]] <- data.frame(
    Rep = rep,
    Method = c("OLS", "Standard MM", "ROBU"),
    Time_Seconds = unname(c(t_ols, t_mm, t_robu)),
    MSE_vs_Baseline = c(mse_ols, mse_mm, mse_robu)
  )
  
  # F. Store residual data from the final replication for the 2-panel plot
  if (rep == n_reps) {
    
    # Safely extract MM residuals (fill with 0s if it completely crashed to avoid ggplot errors)
    res_mm <- if (any(is.na(fit_mm$coefficients))) {
      rep(0, n)
    } else {
      drop(y_cont - (X_cont %*% fit_mm$coefficients))
    }
    
    # Safely extract ROBU residuals
    res_robu <- if (any(is.na(fit_robu$coefficients))) {
      rep(0, n)
    } else {
      drop(y_cont - (X_cont %*% fit_robu$coefficients))
    }
    
    plot_data <- data.frame(
      Patient_ID = rownames(X_cont),
      Is_Contaminated = is_outlier,
      Residuals_CleanBaseline = drop(y_clean - (X_clean %*% beta_clean)),
      Residuals_MM_Contam = res_mm,
      Residuals_ROBU_Contam = res_robu
    )
  }
}

# ______________________________
# 5. Summarize and Save Results
# ______________________________

all_results <- do.call(rbind, results_list)
saveRDS(all_results, "application/results/tcga_full_results.rds")

# Calculate Averages
summary_table <- aggregate(cbind(Time_Seconds, MSE_vs_Baseline) ~ Method, data = all_results, FUN = mean, na.rm = TRUE)

cat("\n======================================================\n")
cat("   PROTEOGENOMICS APPLICATION: AVERAGE OVER 50 REPS   \n")
cat("======================================================\n")
print(summary_table, row.names = FALSE)
cat("======================================================\n")

# Save the summary table
saveRDS(summary_table, "application/results/tcga_performance_table.rds")

# Save residuals for plotting (Standard MM & ROBU vs Baseline)
saveRDS(plot_data, "application/results/tcga_residuals_plotdata.rds")

cat("\nAnalysis complete. Results and plot data saved to application/results/\n")