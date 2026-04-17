#' @title Robust Orthogonal Block Updates (ROBU) for Large-p Regression
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param k Number of groups for partitioning predictor variables. Defaults to block sizes of ~50.
#' @param robu.control Settings for the block updates (from lmrob.control).
#' @param m.control Settings for the final M-estimate (from lmrob.control).
#'
#' @return A list containing the final coefficients, scale, residuals, intermediate estimates, and robust weights.
#' 
#' @export
#' 
robu <- function(x, y, k = max(1, floor(ncol(x) / 50)), robu.control, m.control) {
  
  n <- nrow(x)
  p <- ncol(x)
  
  # Pre-compute constants
  b.const <- m.control$bb * (1 - p/n)
  chi.const <- m.control$tuning.chi
  
  # Partition indices for the predictor variables
  k.ind <- sort((1:p) %% k + 1)
  
  # _______________________________
  # Step 1: Orthogonal Decoupling
  # _______________________________
  
  # Standard QR decomposition 
  x.qr <- qr(x)
  q.mat <- qr.Q(x.qr)
  r.mat <- qr.R(x.qr)
  
  robu.coef.init <- numeric(p)
  y.robu <- y
  
  # ___________________________________________
  # Step 2: Vertical Updates for MM estimates
  # ___________________________________________

  for (k.it in 1:k) {
    x.dat <- q.mat[, k.ind == k.it, drop = FALSE]
    
    # Suppress internal Fast-S warnings caused by zero-intercept orthogonal fitting
    suppressWarnings({
      block.coef <- robustbase::lmrob(formula = y.robu ~ x.dat - 1, control = robu.control)$coef
    })
    
    y.robu <- drop(y.robu - x.dat %*% block.coef)  
    robu.coef.init[k.ind == k.it] <- block.coef
  }
  
  # ________________________________________
  # Step 3: Solving for the initial M-scale
  # ________________________________________

  robu.scale.init <- RobStatTM:::mscale(
    u = drop(y - q.mat %*% robu.coef.init),
    delta = b.const, 
    tuning.chi = chi.const, 
    family = 'bisquare'
  )
  
  # _________________________________________
  # Step 4: IRWLS iterations for S-estimator
  # _________________________________________

  tmp <- irwls.s(
    x = q.mat, 
    y = y, 
    initial.beta = robu.coef.init,
    initial.scale = robu.scale.init,
    max.iter = 200, 
    check.conv = 1, 
    b.const = b.const, 
    c.const = chi.const, 
    tol.conv = 1e-8
  ) 
  
  robu.scale.s <- RobStatTM:::mscale(
    u = drop(y - q.mat %*% tmp$beta.rw),
    delta = b.const, 
    tuning.chi = chi.const, 
    family = 'bisquare'
  )
  
  # _______________________________________
  # Step 5: Final M-estimate using M-scale
  # _______________________________________

  robu.final <- robustbase::lmrob..M..fit(
    x = q.mat, 
    y = y, 
    beta.initial = tmp$beta.rw, 
    scale = robu.scale.s, 
    control = m.control
  )
  
  # _______________________________________________
  # Step 6: Back-transformation to original scale
  # _______________________________________________

  robu.final$coefficients <- backsolve(r.mat, robu.final$coefficients)
  robu.s.coefficients <- backsolve(r.mat, tmp$beta.rw)
  
  robu.final$scale <- RobStatTM:::mscale(
    u = drop(y - x %*% robu.final$coefficients),
    delta = b.const, 
    tuning.chi = chi.const, 
    family = 'bisquare'
  )
  
  return(list(
    coefficients = robu.final$coefficients,
    scale = robu.final$scale,
    scale.initial = robu.scale.init,
    residuals = robu.final$residuals, 
    scale.s = robu.scale.s,
    s.coefficients = robu.s.coefficients,
    rweights = robu.final$rweights
  ))
}