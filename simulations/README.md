# Simulations

This folder contains the scripts required to reproduce the Monte Carlo simulation study (Section 4) and the block-size sensitivity analysis from the manuscript.

## Files

* `generate_data.R`: Contains the data generation mechanism. It simulates high-dimensional predictors with a block-diagonal correlation structure and creates the three specific testing environments (Clean Data, Vertical Outliers, and Concentrated Bad Leverage Points).
* `run_simulations.R`: The master script that loops over $p \in \{100, 200, 400\}$ and the contamination levels $\epsilon \in \{0.15, 0.35\}$ to benchmark the computation time and Mean Squared Error (MSE) of Standard OLS, Standard MM, and ROBU.
* `run_k_sensitivity.R`: The script that evaluates the performance of the ROBU algorithm across a grid of block dimensions ($k$) to demonstrate the optimal trade-off between computational speed and empirical breakdown resistance.

## Dependencies

**1. Core Statistical Packages (Strict versions for reproducibility):**

```R
install.packages("remotes")
remotes::install_version("robustbase", version = "0.99.6")
remotes::install_version("RobStatTM", version = "1.0.11")
remotes::install_version("mvnfast", version = "0.2.8")
```

**2. Plotting and Formatting Packages:**

```R
remotes::install_version("ggplot2", version = "4.0.1")
remotes::install_version("dplyr", version = "1.1.4")
remotes::install_version("xtable", version = "1.8.4")
remotes::install_version("gridExtra", version = "2.3")
remotes::install_version("scales", version = "1.4.0")
```

## Usage

Ensure your working directory is set to the root of the repository. You can source the scripts directly from your R console:

```R
source("simulations/run_simulations.R")
source("simulations/run_k_sensitivity.R")
```

*Note: The scripts will automatically create a `results/` sub-directory and incrementally save the outputs as `.rds` files to prevent data loss during long runs.*
```