#' @title Iteratively Reweighted Least Squares (IRWLS) for S-estimators
#'
#' @param x Design matrix (typically orthogonalized).
#' @param y Response vector.
#' @param initial.beta Initial coefficient vector.
#' @param initial.scale Initial scale estimate.
#' @param max.iter Maximum number of refining steps.
#' @param check.conv 1 to stop when convergence is detected, 0 to run all iterations.
#' @param b.const Tuning constant of the equation.
#' @param c.const Tuning constant of the equation.
#' @param tol.conv Tolerance for convergence detection.
#'
#' @return A list containing the refined beta and scale.
#' 
#' @export
#' 
irwls.s <- function(x, y, initial.beta, initial.scale, max.iter, 
                    check.conv, b.const, c.const, tol.conv) {
  
  # Pre-compute constants to avoid recalculating inside the loop
  cc2 <- c.const^2
  cc4 <- cc2^2
  tol.conv.sq <- tol.conv^2
  
  # _________________
  # Helper Functions
  # _________________
  
  # Optimized loss function (rho) for bisquare
  compute.rho <- function(u) {
    u2 <- u^2
    w <- u2 <= cc2
    v <- numeric(length(u))
    # Only compute polynomial for values within the cutoff
    v[w] <- (u2[w] / 2 * (1 - (u2[w] / cc2) + (u2[w]^2 / (3 * cc4))))
    v[!w] <- cc2 / 6
    return(v * 6 / cc2)
  }
  
  # Optimized weight function for IRWLS
  compute.weights <- function(u) {
    u2.cc2 <- (u^2) / cc2
    w <- u2.cc2 <= 1
    res <- numeric(length(u))
    # Only compute weights for values within the cutoff
    res[w] <- (1 - u2.cc2[w])^2
    return(res)
  }
  
  # ________________
  # Initialization
  # ________________
  
  res <- drop(y - x %*% initial.beta)
  
  if (missing(initial.scale)) {
    scale.val <- median(abs(res)) / 0.6745
  } else {
    scale.val <- initial.scale
  }
  
  beta.val <- initial.beta
  
  # __________________
  # IRWLS Iterations
  # __________________
  
  for (i in 1:max.iter) {
    
    # Step 1: Improve scale
    scale.val <- sqrt(scale.val^2 * mean(compute.rho(res / scale.val)) / b.const)
    
    # Step 2: IRWLS with improved scale
    weights <- compute.weights(res / scale.val)
    
    # Compute X^T W X and X^T W y using fast crossprod and vector recycling.
    # R recycles 'weights' column-by-column to match 'x', completely avoiding
    # the need to allocate an n x p weight matrix or calculate square roots.
    xtwx <- crossprod(x, weights * x)
    xtwy <- crossprod(x, weights * y)
    
    # Safely attempt matrix inversion
    beta.next <- try(solve(xtwx, xtwy), silent = TRUE)
    
    if (inherits(beta.next, "try-error")) {
      beta.next <- initial.beta
      scale.val <- initial.scale
      break
    }
    
    # Check for convergence using squared norms (avoids sqrt computations)
    if (check.conv == 1) {
      diff.sq <- sum((beta.val - beta.next)^2)
      beta.sq <- sum(beta.val^2)
      if ((diff.sq / beta.sq) < tol.conv.sq) {
        break
      }
    }
    
    res <- drop(y - x %*% beta.next)
    beta.val <- drop(beta.next)
  }
  
  return(list(beta.rw = beta.val, scale.rw = scale.val))
}