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

# _________________________
# 1. Simulation Parameters 
# _________________________

set.seed(0)

n.obs <- 1000
p.vars <- 400
n.reps <- 50          # 50 reps is sufficient to show the trend
cont.vec <- c(0.20)   # Updated to 0.20 to align with new main simulation grid

# The grid of k values to test
# k=1 (Unblocked), k=5 (Blocks of 80), k=10 (Blocks of 40), 
# k=20 (Blocks of 20), k=40 (Blocks of 10), k=100 (Blocks of 4), k=400 (Univariate)
k.vec <- c(1, 5, 10, 20, 40, 100, 400)

# Strict fairness control: Give it the massive budget so k=1 (Standard MM) has a fair fight
my.control <- lmrob.control(
  method = "MM", 
  fast.s.large.n = Inf, 
  k.max = 5000,          # Massive budget for p=400
  k.m_s = 2000,
  max.it = 1000,
  nResample = 5000,      # Massive budget for p=400
  refine.tol = 1e-5      
)

# Create results directory if it doesn't exist
if (!dir.exists("simulations/results")) {
  dir.create("simulations/results", recursive = TRUE)
}

# Helper function to evaluate method, silently catch warnings, and track outlier detection
evaluate_robu_k <- function(x.dat, y.dat, k.val, beta.true, out.ind) {
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
  
  # MSE: strictly the squared norm (sum of squared errors)
  mse <- sum((fit$coefficients - beta.true)^2)
  
  # Outlier Detection Tracking
  wt <- fit$rweights
  if (is.null(wt)) wt <- fit$weights # fallback
  
  if (!is.null(wt) && length(wt) == length(out.ind)) {
    is_detected <- wt < 0.5
    is_actual <- out.ind == 1
    tp <- sum(is_detected & is_actual)
    fp <- sum(is_detected & !is_actual)
  } else {
    tp <- NA
    fp <- NA
  }
  
  return(list(
    MSE = mse,
    Time = unname(time_elapsed),
    S_Conv = s_conv,
    M_Conv = m_conv,
    TP = tp,
    FP = fp
  ))
}

# __________________________
# 2. Main Sensitivity Loop
# __________________________

cat("Starting k-Sensitivity Analysis...\n")
cat("Fixed p = 400 | Scenario = Concentrated Bad Leverage Points\n")
cat("Total Configurations:", length(k.vec) * n.reps * length(cont.vec), "\n\n")

for (cont.level in cont.vec) {
  
  cat(sprintf("\n******************************************************\n"))
  cat(sprintf("   COMMENCING CONTAMINATION LEVEL: %.2f\n", cont.level))
  cat(sprintf("******************************************************\n\n"))
  
  for (k.val in k.vec) {
    
    block.size <- floor(p.vars / k.val)
    
    # Initialize results storage for this specific k and eps
    results_file <- sprintf("simulations/results/k_sensitivity_k%03d_eps%.2f.rds", k.val, cont.level)
    
    results <- data.frame(
      K = integer(),
      BlockSize = integer(),
      Contamination = numeric(),
      Rep = integer(),
      MSE = numeric(),
      Time = numeric(),
      S_Conv = logical(),    
      M_Conv = logical(),
      TP = integer(),
      FP = integer()
    )
    
    cat("\n======================================================\n")
    cat(sprintf("Starting Grid for k = %d (Block Size: %d) | eps = %.2f\n", k.val, block.size, cont.level))
    cat(sprintf("Saved to: %s\n", results_file))
    cat("======================================================\n")
    
    for (rep in 1:n.reps) {
      
      cat(sprintf("Running: eps = %.2f | k = %-3d | Rep = %-3d \n", cont.level, k.val, rep))
      
      # Always generate Scenario 3 (Concentrated Leverage Points)
      dat <- generate_data(n = n.obs, p = p.vars, cont.prop = cont.level, leverage = TRUE)
      
      # Run ROBU with the current k
      res.robu <- evaluate_robu_k(dat$x, dat$y, k.val, dat$beta.true, dat$out.ind)
      
      # Bind Results
      new_row <- data.frame(
        K = k.val,
        BlockSize = block.size,
        Contamination = cont.level,
        Rep = rep,
        MSE = res.robu$MSE,
        Time = res.robu$Time,
        S_Conv = res.robu$S_Conv,
        M_Conv = res.robu$M_Conv,
        TP = res.robu$TP,
        FP = res.robu$FP
      )
      
      results <- rbind(results, new_row)
      
      # Save intermediate results incrementally for this k and eps
      saveRDS(results, results_file)
      
    } # End Rep Loop
    
    cat(sprintf("Finished all reps for k = %d at eps = %.2f.\n", k.val, cont.level))
    
  } # End K Loop
} # End Contamination Loop

cat("\nSensitivity Analysis Complete!\n")