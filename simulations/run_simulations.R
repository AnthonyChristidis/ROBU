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
n.reps <- 100
cont.vec <- c(0.15, 0.35)   # Two contamination levels

# Strict fairness control for both Standard MM and ROBU
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

# ____________________________________________________________________
# 2. Helper Function: Evaluate Method and Catch Convergence Warnings
# ____________________________________________________________________

evaluate_method <- function(method_name, expr, beta.true, p.vars) {
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
  
  mse <- mean((coefs - beta.true)^2)
  
  return(list(
    Method = method_name,
    MSE = mse,
    Time = unname(time_elapsed),
    S_Conv = s_conv,
    M_Conv = m_conv
  ))
}

# ________________________
# 3. Main Simulation Loop
# ________________________

cat("Starting Monte Carlo Simulations...\n")

for (p.vars in p.vec) {
  
  # Calculate optimal block size for ROBU dynamically
  k.blocks <- max(1, floor(p.vars / 20)) 
  
  for (cont.level in cont.vec) {
    
    # Initialize results storage for this specific P and Epsilon combination
    results_file <- sprintf("simulations/results/sim_n%d_p%d_eps%.2f.rds", n.obs, p.vars, cont.level)
    
    results <- data.frame(
      N = integer(),
      P = integer(),
      Contamination = numeric(),
      Scenario = character(),
      Rep = integer(),
      Method = character(),
      MSE = numeric(),
      Time = numeric(),
      S_Conv = logical(),
      M_Conv = logical()
    )
    
    cat("\n======================================================\n")
    cat(sprintf("Starting Grid for p = %d | eps = %.2f\n", p.vars, cont.level))
    cat(sprintf("Saved to: %s\n", results_file))
    cat("======================================================\n")
    
    for (scenario in 1:3) {
      
      # Skip redundant "Clean Data" runs (only need to run 0% contamination once per p)
      if (scenario == 1 && cont.level != cont.vec[1]) {
        next
      }
      
      scenario.name <- switch(scenario,
                              "1" = "Clean",
                              "2" = "Vertical Outliers",
                              "3" = "Leverage Points")
      
      actual.cont <- ifelse(scenario == 1, 0, cont.level)
      
      for (rep in 1:n.reps) {
        
        cat(sprintf("Running: p = %-3d | Cont = %-4.2f | Scenario = %-17s | Rep = %-3d \n", 
                    p.vars, actual.cont, scenario.name, rep))
        
        # Generate Data based on scenario
        if (scenario == 1) {
          dat <- generate_data(n = n.obs, p = p.vars, cont.prop = 0)
        } else if (scenario == 2) {
          dat <- generate_data(n = n.obs, p = p.vars, cont.prop = actual.cont, leverage = FALSE)
        } else {
          dat <- generate_data(n = n.obs, p = p.vars, cont.prop = actual.cont, leverage = TRUE)
        }
        
        # A. Standard OLS
        res.ols <- evaluate_method("OLS", lm(dat$y ~ dat$x - 1), dat$beta.true, p.vars)
        res.ols$S_Conv <- NA
        res.ols$M_Conv <- NA
        
        # B. Standard MM
        res.mm <- evaluate_method("Standard MM", 
                                  robustbase::lmrob(dat$y ~ dat$x - 1, control = my.control), 
                                  dat$beta.true, p.vars)
        
        # C. ROBU
        res.robu <- evaluate_method("ROBU", 
                                    robu(x = dat$x, y = dat$y, k = k.blocks, robu.control = my.control, m.control = my.control), 
                                    dat$beta.true, p.vars)
        
        # Bind Results
        new_rows <- data.frame(
          N = n.obs,
          P = p.vars,
          Contamination = actual.cont,
          Scenario = scenario.name,
          Rep = rep,
          Method = c("OLS", "Standard MM", "ROBU"),
          MSE = c(res.ols$MSE, res.mm$MSE, res.robu$MSE),
          Time = c(res.ols$Time, res.mm$Time, res.robu$Time),
          S_Conv = c(res.ols$S_Conv, res.mm$S_Conv, res.robu$S_Conv),
          M_Conv = c(res.ols$M_Conv, res.mm$M_Conv, res.robu$M_Conv)
        )
        
        results <- rbind(results, new_rows)
        
        # Save intermediate results incrementally
        saveRDS(results, results_file)
        
      } # End Rep Loop
    } # End Scenario Loop
  } # End Contamination Loop
  cat(sprintf("\nFinished all scenarios for p = %d.\n", p.vars))
} # End P Loop

cat("\nAll Simulations Complete!\n")