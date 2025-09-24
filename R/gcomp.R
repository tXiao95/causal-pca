gcomp <- function(Y, X, C, X.new = NULL, ...) {
  n <- length(Y); p <- ncol(X); q <- ncol(C)
  
  # Create full covariate matrix (X, C)
  C <- data.frame(C); colnames(C) <- paste0("C", 1:q)
  X <- data.frame(X)
  
  df <- cbind(X, C)
  
  # Fit outcome regression: E[Y | X, C]
  sl_fit <- SuperLearner::SuperLearner(Y = Y, X = df, ...)
  
  # If no new X values provided, just solve for the previous ones
  if (is.null(X.new)) {
    X.new <- X
  } else{
    X.new <- data.frame(X.new)
    if (ncol(X.new) != p) stop("X.new must have the same number of columns as X.")
  }
  m <- nrow(X.new)
  
  # Estimate each exposure value over entire confounder distribution 
  gcomp_est <- vapply(1:m, function(i){
    Xi.new <- X.new[i, , drop = FALSE]
    df.new <- cbind(Xi.new[rep(1, n), , drop = FALSE], C)
    mean( predict(sl_fit, newdata = df.new)$pred )
  }, numeric(1L) )
  
  names(gcomp_est) <- paste("x", 1:m)
  return(gcomp_est)
}
