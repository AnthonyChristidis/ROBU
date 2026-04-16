# --------------------------------
# ROBU: Master Simulation Script
# --------------------------------

# Load necessary packages
library(robustbase)
library(RobStatTM)
library(mvnfast)

# Source algorithm and data generation scripts
source("R/irwls.R")
source("R/robu.R")
source("simulations/generate_data.R")

# __________________________
# 1. Simulation Parameters 
# __________________________

set.seed(0)

n.obs <- 1000
p.vec <- c(100, 200, 400)
cont.vec <- c(0.10, 0.20, 0.30)   # Expanded epsilons
scenario.vec <- 1:3               # 1=Clean, 2=Vertical, 3=Leverage
n.reps <- 50

# Create results directory if it doesn't exist
if (!dir.exists("simulations/results")) {
  dir.create("simulations/results", recursive = TRUE)
}

# ____________________________________________________________________
# 2. Helper Function: Evaluate Method and Catch Convergence Warnings
# ____________________________________________________________________

evaluate_method <- function(method_name, expr, beta.true, p.vars, out.ind) {
  s_conv <- TRUE
  m_conv <- TRUE
  
  start_time <- proc.time()["elapsed"]
  
  fit <- tryCatch({
    withCallingHandlers(
      expr,
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
  
  if (inherits(fit, "lm")) {
    coefs <- coef(fit)
  } else {
    coefs <- fit$coefficients
  }
  
  # MSE: strictly the squared norm (sum of squared errors), NOT divided by p
  mse <- sum((coefs - beta.true)^2)
  
  # Outlier Detection Tracking (TP/FP)
  if (method_name == "OLS") {
    tp <- NA
    fp <- NA
  } else {
    # Extract weights (standardized to usually be < 0.5 for rejected outliers)
    wt <- fit$rweights
    if (is.null(wt)) wt <- fit$weights # fallback
    
    if (!is.null(wt)) {
      is_detected <- wt < 0.5
      is_actual <- out.ind == 1
      tp <- sum(is_detected & is_actual)
      fp <- sum(is_detected & !is_actual)
    } else {
      tp <- NA
      fp <- NA
    }
  }
  
  return(list(
    Method = method_name,
    MSE = mse,
    Time = unname(time_elapsed),
    S_Conv = s_conv,
    M_Conv = m_conv,
    TP = tp,
    FP = fp
  ))
}

# ________________________
# 3. Main Simulation Loop
# ________________________

cat("Starting Monte Carlo Simulations...\n")

for (p.vars in p.vec) {
  
  # Dynamic computational budget to give Standard MM a fair fight
  k.max_val <- ifelse(p.vars >= 400, 5000, 2000)
  nResample_val <- ifelse(p.vars >= 400, 5000, 2000)
  
  mm.control <- lmrob.control(
    method = "MM", 
    fast.s.large.n = Inf, 
    k.max = k.max_val,
    k.m_s = 2000,
    max.it = 1000,
    nResample = nResample_val,
    refine.tol = 1e-5      
  )
  
  # Calculate optimal block size for ROBU dynamically
  k.blocks <- max(1, floor(p.vars / 20)) 
  
  for (cont.level in cont.vec) {
    
    for (scenario in scenario.vec) {
      
      # Skip redundant "Clean Data" runs (only need to run 0% contamination once per p)
      if (scenario == 1 && cont.level != cont.vec[1]) {
        next
      }
      
      scenario.name <- switch(as.character(scenario),
                              "1" = "Clean",
                              "2" = "Vertical Outliers",
                              "3" = "Leverage Points")
                              
      scenario.file.suffix <- switch(as.character(scenario),
                                     "1" = "clean",
                                     "2" = "vertical",
                                     "3" = "leverage")
      
      actual.cont <- ifelse(scenario == 1, 0, cont.level)
      
      results_file <- sprintf("simulations/results/sim_n%d_p%d_eps%.2f_%s.rds", 
                              n.obs, p.vars, actual.cont, scenario.file.suffix)
      
      results <- data.frame()
      
      cat("\n======================================================\n")
      cat(sprintf("Starting Grid: p = %d | eps = %.2f | Scenario = %s\n", p.vars, actual.cont, scenario.name))
      cat(sprintf("Saved to: %s\n", results_file))
      cat("======================================================\n")
      
      for (rep in 1:n.reps) {
        
        cat(sprintf("Running: p = %-3d | Cont = %-4.2f | Scenario = %-17s | Rep = %-3d \n", 
                    p.vars, actual.cont, scenario.name, rep))
        
        # Generate Data
        if (scenario == 1) {
          dat <- generate_data(n = n.obs, p = p.vars, cont.prop = 0)
        } else if (scenario == 2) {
          dat <- generate_data(n = n.obs, p = p.vars, cont.prop = actual.cont, leverage = FALSE)
        } else {
          dat <- generate_data(n = n.obs, p = p.vars, cont.prop = actual.cont, leverage = TRUE)
        }
        
        # A. Standard OLS
        res.ols <- evaluate_method("OLS", lm(dat$y ~ dat$x - 1), dat$beta.true, p.vars, dat$out.ind)
        
        # B. Standard MM (with massive computational budget)
        res.mm <- evaluate_method("Standard MM", 
                                  robustbase::lmrob(dat$y ~ dat$x - 1, control = mm.control), 
                                  dat$beta.true, p.vars, dat$out.ind)
        
        # C. Deterministic MM (RobStatTM)
        res.det <- evaluate_method("Deterministic MM", 
                                   RobStatTM::lmrobdetMM(dat$y ~ dat$x - 1), 
                                   dat$beta.true, p.vars, dat$out.ind)
        
        # D. ROBU (Using identical budget/control for strict fairness)
        res.robu <- evaluate_method("ROBU", 
                                    robu(x = dat$x, y = dat$y, k = k.blocks, robu.control = mm.control, m.control = mm.control), 
                                    dat$beta.true, p.vars, dat$out.ind)
        
        # Bind Results
        new_rows <- data.frame(
          N = n.obs, P = p.vars, Contamination = actual.cont, Scenario = scenario.name, Rep = rep,
          Method = c("OLS", "Standard MM", "Deterministic MM", "ROBU"),
          MSE = c(res.ols$MSE, res.mm$MSE, res.det$MSE, res.robu$MSE),
          Time = c(res.ols$Time, res.mm$Time, res.det$Time, res.robu$Time),
          S_Conv = c(res.ols$S_Conv, res.mm$S_Conv, res.det$S_Conv, res.robu$S_Conv),
          M_Conv = c(res.ols$M_Conv, res.mm$M_Conv, res.det$M_Conv, res.robu$M_Conv),
          TP = c(res.ols$TP, res.mm$TP, res.det$TP, res.robu$TP),
          FP = c(res.ols$FP, res.mm$FP, res.det$FP, res.robu$FP)
        )
        
        results <- rbind(results, new_rows)
        saveRDS(results, results_file)
        
      } # End Rep Loop
    } # End Scenario Loop
  } # End Contamination Loop
} # End P Loop

cat("\nAll Simulations Complete!\n")