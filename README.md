# ROBU: Fast Robust Regression via Orthogonal Block Updates

**Authors:** Anthony Christidis and Matías Salibián-Barrera

This repository contains the R implementation of the **ROBU** algorithm, as well as all code necessary to reproduce the simulation studies and real data applications from the manuscript *"Fast Robust Regression via Orthogonal Block Updates"* by Anthony Christidis and Matías Salibián-Barrera.

## Overview

Standard high-breakdown robust regression methods (such as MM-estimators via the Fast-S algorithm) scale poorly as the number of predictors ($p$) grows large, frequently experiencing empirical breakdown or severe computational bottlenecks. 

**ROBU** solves this by decoupling the predictors using orthogonalization (QR decomposition) and applying robust block-coordinate descent. This preserves the 50% breakdown point of the MM-estimator while reducing computation time by orders of magnitude for large $p$ (where $p < n$).

## Dependencies
The code relies on standard robust statistics packages in R. Please ensure the following are installed:
```R
install.packages(c("robustbase", "RobStatTM", "mvnfast"))
```

## Repository Structure

* **`R/`**: Contains the core algorithm functions (`robu.R` and `irwls.R`). These are all you need to apply ROBU to your own data.
* **`simulations/`**: Contains the data generation functions and the scripts to reproduce the Monte Carlo simulation grid (computation time and MSE) presented in the paper.
* **`application/`**: Contains the code for the real-world data application, including the artificial leverage point contamination and benchmarking.

## Quick Start

To use the ROBU algorithm on your own data, simply source the functions in the `R/` directory:

```R
# Load dependencies
library(robustbase)
library(RobStatTM)

# Source the algorithm
source("R/irwls.R")
source("R/robu.R")

# Define control settings for the block updates and final M-estimate
my.control <- lmrob.control(method = "MM", fast.s.large.n = Inf)

# Fit the ROBU model
# k is the number of blocks (defaults to block sizes of ~50 variables)
fit <- robu(x = design_matrix, 
            y = response_vector, 
            robu.control = my.control, 
            m.control = my.control)

# Access results
fit$coefficients
fit$scale
```

## Reproducing the Paper

To reproduce the numerical experiments from the manuscript:

1. Navigate to the `simulations/` folder.
2. Run `run_simulations.R` to execute the grid of scenarios comparing standard OLS, standard MM (`lmrob`), and ROBU.
3. Navigate to the `application/` folder and run the analysis script to reproduce the real data benchmark.