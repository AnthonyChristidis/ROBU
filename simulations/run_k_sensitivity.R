# ------------------------------------------------------
# ROBU: Sensitivity Analysis for Block Dimensions (k)
# ------------------------------------------------------

# Load necessary packages
library(robustbase)
library(RobStatTM)
library(mvnfast)

# Source algorithm and data generation scripts
source("R/irwls.R")
source("R/robu.R")
source("simulations/generate_data.R")

# __________________________________________________
# 1. Simulation Parameters (Matched to Section S5)
# __________________________________________________

set.seed(0)

n.obs <- 1000
p.vars <- 400
n.reps <- 50          # 50 reps is sufficient to show the trend
cont.level <- 0.15

# The grid of k values to test
# k=1 (Unblocked), k=5 (Blocks of 80), k=10 (Blocks of 40), 
# k=20 (Blocks of 20), k=40 (Blocks of 10), k=100 (Blocks of 4), k=400 (Univariate)
k.vec <- c(1, 5, 10, 20, 40, 100, 400)

# Strict fairness control 
my.control <- lmrob.control(
  method = "MM", 
  fast.s.large.n = Inf, 
  k.max = 200,           
  refine.tol = 1e-5      
)

# Create results directory if it doesn't exist
if (!dir.exists("simulations/results")) {
  dir.create("simulations/results", recursive = TRUE)
}

# Helper function to evaluate method and silently catch warnings
evaluate_robu_k <- function(x.dat, y.dat, k.val, beta.true) {
  s_conv <- TRUE
  m_conv <- TRUE
  
  start_time <- proc.time()["elapsed"]
  
  fit <- tryCatch({
    withCallingHandlers(
      robu(x = x.dat, y = y.dat, k = k.val, robu.control = my.control, m.control = my.control),
      warning = function(w) {
        msg <- w$message
        if (grepl("S refinements did not converge", msg) || grepl("S-step", msg)) {
          s_conv <<- FALSE
          invokeRestart("muffleWarning") 
        }
        if (grepl("not converged", msg) && !grepl("S refinements", msg)) {
          m_conv <<- FALSE
          invokeRestart("muffleWarning") 
        }
      }
    )
  }, error = function(e) {
    s_conv <<- FALSE
    m_conv <<- FALSE
    return(list(coefficients = rep(NA, p.vars)))
  })
  
  end_time <- proc.time()["elapsed"]
  time_elapsed <- end_time - start_time
  
  mse <- mean((fit$coefficients - beta.true)^2)
  
  return(list(
    MSE = mse,
    Time = unname(time_elapsed),
    S_Conv = s_conv,
    M_Conv = m_conv
  ))
}

# __________________________
# 2. Main Sensitivity Loop
# __________________________

cat("Starting k-Sensitivity Analysis...\n")
cat("Fixed p = 400 | Scenario = Concentrated Bad Leverage Points\n")
cat("Total Configurations:", length(k.vec) * n.reps, "\n\n")

for (k.val in k.vec) {
  
  block.size <- floor(p.vars / k.val)
  
  # Initialize results storage for this specific k
  # Using %03d adds leading zeros (e.g., k001, k020) so files sort perfectly in your folder
  results_file <- sprintf("simulations/results/k_sensitivity_k%03d.rds", k.val)
  
  results <- data.frame(
    K = integer(),
    BlockSize = integer(),
    Rep = integer(),
    MSE = numeric(),
    Time = numeric(),
    S_Conv = logical(),    
    M_Conv = logical()     
  )
  
  cat("\n======================================================\n")
  cat(sprintf("Starting Grid for k = %d (Block Size: %d)\n", k.val, block.size))
  cat(sprintf("Saved to: %s\n", results_file))
  cat("======================================================\n")
  
  for (rep in 1:n.reps) {
    
    cat(sprintf("Running: k = %-3d | Rep = %-3d \n", k.val, rep))
    
    # Always generate Scenario 3 (Concentrated Leverage Points)
    dat <- generate.data(n = n.obs, p = p.vars, cont.prop = cont.level, leverage = TRUE)
    
    # Run ROBU with the current k
    res.robu <- evaluate_robu_k(dat$x, dat$y, k.val, dat$beta.true)
    
    # Bind Results
    new_row <- data.frame(
      K = k.val,
      BlockSize = block.size,
      Rep = rep,
      MSE = res.robu$MSE,
      Time = res.robu$Time,
      S_Conv = res.robu$S_Conv,
      M_Conv = res.robu$M_Conv
    )
    
    results <- rbind(results, new_row)
    
    # Save intermediate results incrementally for this k
    saveRDS(results, results_file)
    
  } # End Rep Loop
  
  cat(sprintf("Finished all reps for k = %d.\n", k.val))
  
} # End K Loop

cat("\nSensitivity Analysis Complete!\n")