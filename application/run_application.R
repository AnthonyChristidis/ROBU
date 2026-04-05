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
set.seed(2024)

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

# ________________________________________
# 4. Introduce Adversarial Contamination
# ________________________________________

cat("\n--- Injecting Adversarial Leverage Points ---\n")

X_cont <- X_clean
y_cont <- y_clean

cont_prop <- 0.15
n_cont <- floor(cont_prop * n)

# Randomly select 15% of patients and 30 genes to corrupt
bad_patients <- sample(1:n, n_cont)
bad_genes <- sample(1:p, 30)

# Corrupt X: Add massive noise N(100, 10^2)
X_cont[bad_patients, bad_genes] <- X_cont[bad_patients, bad_genes] + 
  matrix(rnorm(n_cont * 30, mean = 100, sd = 10), nrow = n_cont, ncol = 30)

# Corrupt y: Generate adversarial response
beta_bad <- rnorm(p, mean = -5, sd = 2)
y_cont[bad_patients] <- X_cont[bad_patients, ] %*% beta_bad + rnorm(n_cont, 0, 1)

# Create a boolean vector of true outliers for our plots later
is_outlier <- rep(FALSE, n)
is_outlier[bad_patients] <- TRUE

# ______________________________________________
# 5. Benchmarking Methods on Contaminated Data
# ______________________________________________

cat("\n--- Benchmarking Algorithms on Contaminated Data ---\n")

# Helper to calculate MSE relative to clean baseline
calc_mse <- function(beta_est) {
  mean((beta_est - beta_clean)^2)
}

# A. Standard OLS
cat("Running OLS...\n")
t_ols <- system.time({
  fit_ols <- lm(y_cont ~ X_cont - 1)
})["elapsed"]
mse_ols <- calc_mse(coef(fit_ols))

# B. Standard MM
cat("Running Standard MM...\n")
t_mm <- system.time({
  fit_mm <- tryCatch({
    robustbase::lmrob(y_cont ~ X_cont - 1, control = my.control)
  }, error = function(e) list(coefficients = rep(NA, p), residuals = rep(NA, n)))
})["elapsed"]
mse_mm <- calc_mse(coef(fit_mm))

# C. ROBU
cat("Running ROBU...\n")
k.blocks <- max(1, floor(p / 20))
t_robu <- system.time({
  fit_robu <- robu(x = X_cont, y = y_cont, k = k.blocks, 
                   robu.control = my.control, m.control = my.control)
})["elapsed"]
mse_robu <- calc_mse(fit_robu$coefficients)

# ______________________________
# 6. Summarize and Save Results
# ______________________________

results_table <- data.frame(
  Method = c("OLS", "Standard MM", "ROBU"),
  Time_Seconds = c(t_ols, t_mm, t_robu),
  MSE_vs_Baseline = c(mse_ols, mse_mm, mse_robu)
)

cat("\n========================================\n")
cat("      PROTEOGENOMICS APPLICATION RESULTS  \n")
cat("========================================\n")
print(results_table, row.names = FALSE)
cat("========================================\n")

# Save the table
saveRDS(results_table, "application/results/tcga_performance_table.rds")

# Save residuals for plotting (ROBU vs Baseline to show identical outlier detection)
plot_data <- data.frame(
  Patient_ID = rownames(X_cont),
  Is_Contaminated = is_outlier,
  Residuals_CleanBaseline = y_clean - (X_clean %*% beta_clean),
  Residuals_ROBU_Contam = fit_robu$residuals
)
saveRDS(plot_data, "application/results/tcga_residuals_plotdata.rds")

cat("\nAnalysis complete. Results and plot data saved to application/results/\n")