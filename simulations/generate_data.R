#' @title Simulate Contaminated Data for ROBU Monte Carlo Simulations
#'
#' @param n Total sample size (Default 1000 as per manuscript).
#' @param p Number of predictors.
#' @param snr Signal-to-noise ratio.
#' @param cont.prop Proportion of contamination (0 to 1).
#' @param leverage Logical; if TRUE, generates concentrated bad leverage points. 
#'                 If FALSE, generates vertical outliers only.
#' @param rho Correlation between variables within the same block.
#' @param cov.blocks Number of blocks for the correlation structure.
#'
#' @return A list containing x, y, beta.true, and out.ind (outlier indicators).
#' 
#' @export
#' 
library(mvnfast)

generate_data <- function(n = 1000, p = 200, snr = 3, cont.prop = 0.15, 
                          leverage = TRUE, rho = 0.5, cov.blocks = 10) {
  
  # ____________________________________
  # Step 1: Block-Correlation Structure
  # ____________________________________
  
  sigma.mat <- matrix(0, p, p)
  block.size <- floor(p / cov.blocks)
  
  for (i in 1:cov.blocks) {
    start.ind <- (i - 1) * block.size + 1
    end.ind <- if (i == cov.blocks) p else i * block.size
    
    size.current <- end.ind - start.ind + 1
    block.mat <- matrix(rho, size.current, size.current)
    diag(block.mat) <- 1
    
    sigma.mat[start.ind:end.ind, start.ind:end.ind] <- block.mat
  }
  
  # ______________________________
  # Step 2: Clean Data Generation
  # ______________________________
  
  n.cont <- floor(n * cont.prop)
  n.clean <- n - n.cont
  
  x.clean <- mvnfast::rmvn(n = n.clean, mu = rep(0, p), sigma = sigma.mat)
  beta.true <- rnorm(p)
  
  signal <- x.clean %*% beta.true
  var.signal <- var(drop(signal))
  var.error <- var.signal / snr
  sd.error <- sqrt(var.error)
  
  y.clean <- signal + rnorm(n.clean, mean = 0, sd = sd.error)
  
  # _________________________________
  # Step 3: Contamination Generation
  # _________________________________
  
  if (n.cont > 0) {
    
    if (leverage) {
      # CONCENTRATED BAD LEVERAGE POINTS (Scenario 3)
      x.cont <- mvnfast::rmvn(n = n.cont, mu = rep(0, p), sigma = sigma.mat)
      p.bad <- max(1, floor(p * 0.1)) # Corrupt the first 10% of variables
      
      shift.mat <- mvnfast::rmvn(n = n.cont, mu = rep(100, p.bad), sigma = diag(10^2, p.bad))
      x.cont[, 1:p.bad] <- x.cont[, 1:p.bad] + shift.mat
      
      beta.bad <- rnorm(p, mean = -5, sd = 2)
      y.cont <- x.cont %*% beta.bad + rnorm(n.cont, mean = 0, sd = sd.error)
      
    } else {
      # VERTICAL OUTLIERS (Scenario 2)
      x.cont <- mvnfast::rmvn(n = n.cont, mu = rep(0, p), sigma = sigma.mat)
      y.cont <- (x.cont %*% beta.true) + rnorm(n.cont, mean = 200 * sd.error, sd = sd.error)
    }
    
    x <- rbind(x.clean, x.cont)
    y <- c(drop(y.clean), drop(y.cont))
    out.ind <- c(rep(0, n.clean), rep(1, n.cont))
    
  } else {
    
    # CLEAN DATA (Scenario 1)
    x <- x.clean
    y <- drop(y.clean)
    out.ind <- rep(0, n)
    
  }
  
  return(list(
    x = x,
    y = y,
    beta.true = beta.true,
    out.ind = out.ind
  ))
}