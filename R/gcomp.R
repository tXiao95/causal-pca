gcomp <- function(Y, X, C, X.new = NULL, 
                               SL.library = c("SL.glm", 
                                              "SL.gam", 
                                              "SL.svm",
                                              "SL.xgboost",
                                              "SL.glmnet")) {
  n <- length(Y); p <- ncol(X); q <- ncol(C)
  
  # Create full covariate matrix (X, C)
  C <- data.frame(C); colnames(C) <- paste0("C", 1:q)
  X <- data.frame(X)
  
  df <- cbind(X, C)
  
  # Fit outcome regression: E[Y | X, C]
  sl_fit <- SuperLearner::SuperLearner(Y = Y, 
                                       X = df, 
                                       SL.library = SL.library, 
                                       family = gaussian())
  
  # If no new X values provided, just solve for the previous ones
  if (is.null(X.new)) {
    X.new <- X
  } else{
    X.new <- data.frame(X.new)
  }
  m <- nrow(X.new)
  
  # Create each combination of confounders and covariates. 
  C.block <- C[rep(1:n, times = m), , drop = FALSE]
  X.block <- X.new[rep(1:m, each = n), , drop = FALSE]
  
  df.new  <- cbind(X.block, C.block)
  
  # Predict and average the 'n' E[Y | X=x, C_i] over the confounders at each x value. 
  mu_hat_block <- predict(sl_fit, newdata = df.new)$pred
  mu_mat       <- matrix(mu_hat_block, nrow = n, ncol = m)
  gcomp_est    <- colMeans(mu_mat)
  
  names(gcomp_est) <- paste("x", 1:m)
  return(gcomp_est)
}
