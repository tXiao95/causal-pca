pCCA <- function(Y, X, C, center = TRUE, scale = FALSE) {
  # X: n x p, Y: n x 1 (or n x q), C: n x r
  X <- scale(as.matrix(X), center = center, scale = scale)
  Y <- scale(as.matrix(Y), center = center, scale = scale)
  C <- as.matrix(C)
  
  # Check sanity
  n <- nrow(X)
  stopifnot(n == nrow(Y), n == nrow(C))
  
  # Residual-maker for [1, C] using pivoted QR: Mz = I - Q Q'
  Z   <- cbind(1, C)
  qrZ <- qr(Z)
  Q   <- qr.Q(qrZ, complete = FALSE)         # n x rank(Z)
  Pz  <- if (ncol(Q) > 0) Q %*% t(Q) else matrix(0, n, n)
  Mz  <- diag(n) - Pz
  
  # Residualize X and Y on C
  Xr <- Mz %*% X
  Yr <- Mz %*% as.matrix(Y)
  
  # Canonical correlation analysis on residuals
  out <- cancor(Xr, Yr)
  
  # Return just what you asked for (loadings), plus rho
  return(out$xcoef)      # canonical weights for X|C
}