# ----------------------------------------------------------------
# ROBU: Supplementary S4 - Severe Contamination (35%) Experiment
# ----------------------------------------------------------------

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
p.vars <- 200        
n.reps <- 50         # 50 reps is enough to definitively show the breakdown/recovery
cont.level <- 0.35   # 35% Severe Contamination

# Create results directory if it doesn't exist
if (!dir.exists("simulations/results")) {
  dir.create("simulations/results", recursive = TRUE)
}
results_file <- "simulations/results/severe_contamination_eps0.35.rds"

# ______________________
# 2. Algorithmic Tuning 
# ______________________

# Standard MM Control: Budget of N_s = 2,500
# (Mathematically guaranteed to fail)
ctrl_mm <- lmrob.control(
  method = "MM", 
  fast.s.large.n = Inf, 
  k.max = 200,           
  refine.tol = 1e-5,
  nResample = 2500       # INFLATED BUDGET
)

# ROBU Control: Budget of N_s = 25,000 
# (Mathematically guaranteed to succeed at k=10, since block size is 20)
ctrl_robu <- lmrob.control(
  method = "MM", 
  fast.s.large.n = Inf, 
  k.max = 200,           
  refine.tol = 1e-5,
  nResample = 25000      # HIGHLY INFLATED BUDGET
)

# Initialize results storage
results <- data.frame(
  Method = character(),
  Rep = integer(),
  MSE = numeric(),
  Time = numeric()
)

# _________________________
# 3. Main Evaluation Loop
# _________________________

cat(sprintf("\n******************************************************\n"))
cat(sprintf(" COMMENCING SEVERE CONTAMINATION EXPERIMENT (eps=0.35)\n"))
cat(sprintf("******************************************************\n\n"))

for (rep in 1:n.reps) {
  
  cat(sprintf("Running Rep %d / %d ...\n", rep, n.reps))
  
  # Generate Scenario 3 (Concentrated Leverage Points) at 35%
  dat <- generate_data(n = n.obs, p = p.vars, cont.prop = cont.level, leverage = TRUE)
  
  # _________________________
  # Standard MM (Ns = 2500)
  # _________________________

  cat("  -> Running Standard MM (Ns = 2500)... ")
  t_start <- proc.time()["elapsed"]
  fit_mm <- tryCatch({
    lmrob(dat$y ~ dat$x - 1, control = ctrl_mm)
  }, error = function(e) list(coefficients = rep(NA, p.vars)))
  t_end <- proc.time()["elapsed"]
  
  mse_mm <- mean((coef(fit_mm) - dat$beta.true)^2, na.rm = TRUE)
  time_mm <- unname(t_end - t_start)
  cat(sprintf("Time: %.1fs | MSE: %.2f\n", time_mm, mse_mm))
  
  results <- rbind(results, data.frame(Method="Standard MM", Rep=rep, MSE=mse_mm, Time=time_mm))
  
  # ___________________________
  # ROBU (k = 10, Ns = 25000)
  # ___________________________

  cat("  -> Running ROBU (k = 10, Ns = 25000)... ")
  t_start <- proc.time()["elapsed"]
  fit_robu <- tryCatch({
    suppressWarnings(robu(x = dat$x, y = dat$y, k = 10, robu.control = ctrl_robu, m.control = ctrl_robu))
  }, error = function(e) list(coefficients = rep(NA, p.vars)))
  t_end <- proc.time()["elapsed"]
  
  mse_robu <- mean((fit_robu$coefficients - dat$beta.true)^2, na.rm = TRUE)
  time_robu <- unname(t_end - t_start)
  cat(sprintf("Time: %.1fs | MSE: %.2f\n", time_robu, mse_robu))
  
  results <- rbind(results, data.frame(Method="ROBU", Rep=rep, MSE=mse_robu, Time=time_robu))
  
  # Save incrementally
  saveRDS(results, results_file)
}

cat("\n Supplementary S4 Severe Contamination Experiment Complete!\n")
cat(sprintf("Results saved to: %s\n", results_file))
