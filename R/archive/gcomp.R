gcomp <- function(Y, X, C, X.new = NULL, ...) {
  n <- length(Y); p <- ncol(X); q <- ncol(C)
  
  # Need to relabel C since default names are X1,X2,...
  C <- data.frame(C); colnames(C) <- paste0("C", 1:q)
  X <- data.frame(X); colnames(X) <- paste0("X", 1:p)
  
  # Create full covariate matrix (X, C)
  df <- cbind(X, C)
  
  # Fit outcome regression: E[Y | X, C]
  sl_fit <- SuperLearner::SuperLearner(Y = Y, X = df, ...)
  
  # If no new X values provided, just solve for the previous ones
  if (is.null(X.new)) {
    X.new <- X
  } else{
    X.new <- data.frame(X.new); colnames(X.new) <- paste0("X", 1:p)
    if (ncol(X.new) != p) stop("X.new must have the same number of columns as X.")
  }
  m <- nrow(X.new)
  
  # Estimate each exposure value over entire confounder distribution 
  gcomp_est <- vapply(1:m, function(i){
    print(i)
    Xi.new <- X.new[i, , drop = FALSE]
    df.new <- cbind(Xi.new[rep(1, n), , drop = FALSE], C)
    mean( predict(sl_fit, newdata = df.new)$pred )
  }, numeric(1L) )
  
  names(gcomp_est) <- paste("x", 1:m)
  return(gcomp_est)
}

#' G-Computation Estimator   sum E(Y | X = x_new, C = C_new)
#' 
#' @param out_model An S3 object of class "outcome_model".
#' @param C Numeric matrix or data frame of observed confounders.
#' @param X_new Optional numeric matrix or data frame of new treatment values to evaluate. Defaults to NULL.
#' @return A numeric vector of marginalized expected outcomes.

gcomp2 <- function(outcome_model, C = NULL, X_new = NULL) {
  p <- outcome_model$p
  
  # Default to observed C if C_new is not provided
  if (is.null(C)) {
    C <- outcome_model$C_observed
  } else {
    C <- as.data.frame(C)
    colnames(C) <- outcome_model$C_names
  }
  n <- nrow(C)
  
  # Default to observed X if X_new is not provided
  if (is.null(X_new)) {
    X_new <- outcome_model$X_observed
  } else {
    X_new <- as.data.frame(X_new)
    colnames(X_new) <- outcome_model$X_names
  }
  
  if (ncol(X_new) != p) {
    stop("X_new must have the same number of columns as the training X.")
  }
  m <- nrow(X_new)
  
  # Estimate each exposure value over the entire observed confounder distribution 
  gcomp_est <- vapply(1:m, function(i) {
    print(i)
    Xi_new <- X_new[i, , drop = FALSE]
    
    # Replicate this treatment assignment for every individual
    Xi_rep <- Xi_new[rep(1, n), , drop = FALSE]
    df_new <- cbind(Xi_rep, C)
    
    # Predict and take the empirical mean
    mean(predict(outcome_model, newdata = df_new))
  }, numeric(1L))
  
  names(gcomp_est) <- paste("x", 1:m, sep = "_")
  return(gcomp_est)
}