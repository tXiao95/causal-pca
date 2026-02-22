#' Fit a MultiGPS Estimator
#'
#' @param t_eval A numeric vector of length d_T representing the target treatment value.
#' @param X A numeric matrix (n x d_T) of observed continuous treatments.
#' @param C A numeric matrix (n x d_C) of observed confounders.
#' @param h1 A numeric scalar for the bandwidth.
#' @param mu_fitter A custom function(W, C) that trains a model and returns a fitted object.
#' @return An S3 object of class "multigps".

multigps <- function(t_eval, X, C, h1, mu_fitter) {
  n <- nrow(X)
  d_T <- ncol(X)
  
  # 1. Compute the pseudo-outcome W_i for all observations
  W <- rep(1, n)
  for (j in 1:d_T) {
    u <- (X[, j] - t_eval[j]) / h1
    W <- W * dnorm(u) # Standard Gaussian kernel
  }
  
  # 2. Fit the outcome regression using the user-supplied ML wrapper
  # This fit object should ideally support standard predict() methods later
  inner_fit <- mu_fitter(W, C)
  
  # 3. Store the components needed for future predictions
  res <- list(
    inner_fit = inner_fit,
    t_eval = t_eval,
    h1 = h1,
    d_T = d_T
  )
  
  # Assign the S3 class
  class(res) <- "multigps"
  
  return(res)
}