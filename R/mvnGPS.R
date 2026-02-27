library(mvtnorm)

#' Fit a Global Parametric/Semiparametric MVN GPS Estimator
#'
#' @param X A numeric matrix or data frame of observed continuous treatments (n x p).
#' @param C A numeric matrix or data frame of observed confounders (n x q).
#' @param method String indicating how to model the conditional mean E[X|C]. 
#'               "linear" uses OLS; "SuperLearner" loops over X dimensions flexibly.
#' @param ... Additional arguments passed to SuperLearner (e.g., SL.library, cvControl).
#' @return An S3 object of class "mvnGPS".

mvnGPS <- function(X, C, method = c("linear", "SuperLearner"), ...) {
  method <- match.arg(method)
  
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  
  p <- ncol(X_df)
  q <- ncol(C_df)
  
  # Standardize column names
  colnames(X_df) <- paste0("X", 1:p)
  colnames(C_df) <- paste0("C", 1:q)
  
  X_mat <- as.matrix(X_df)
  
  if (method == "linear") {
    # Fit multivariate regression: X ~ C
    inner_fit <- lm(X_mat ~ ., data = C_df)
    resids <- residuals(inner_fit)
    
  } else if (method == "SuperLearner") {
    # Fit a separate SuperLearner model for each dimension of X
    inner_fit <- list()
    resids <- matrix(NA, nrow = nrow(X_mat), ncol = p)
    
    for (j in 1:p) {
      # SuperLearner requires the outcome to be a vector
      sl_fit <- SuperLearner::SuperLearner(Y = X_mat[, j], 
                                           X = C_df, 
                                           family = gaussian(), 
                                           ...)
      
      inner_fit[[paste0("X", j)]] <- sl_fit
      # SL.predict contains the predictions on the training data
      resids[, j] <- X_mat[, j] - sl_fit$SL.predict
    }
  }
  
  # Estimate the covariance matrix (Sigma) from the residuals
  sigma_hat <- as.matrix(cov(as.matrix(resids)))
  
  res <- list(
    inner_fit = inner_fit,
    sigma_hat = sigma_hat,
    method = method,
    p = p,
    X_names = colnames(X_df),
    C_names = colnames(C_df)
  )
  
  class(res) <- "mvnGPS"
  return(res)
}

#' Predict Method for MVN GPS Objects
#'
#' @param object An object of class "mvnGPS".
#' @param newdata A data frame containing BOTH target treatments (X) and confounders (C).
#' @param delta_n A small positive numeric threshold to floor the density estimates. Defaults to 1e-4.
#' @param ... Additional arguments (ignored, but kept for S3 consistency).
#' @return A numeric vector of estimated conditional densities f(X | C).

predict.mvnGPS <- function(object, newdata, delta_n = 1e-4, ...) {
  newdata_df <- as.data.frame(newdata)
  
  X_new <- newdata_df[, 1:object$p, drop = FALSE] |> as.matrix()
  C_new <- newdata_df[, (object$p+1):ncol(newdata_df), drop = FALSE]
  
  names(X_new) <- object$X_names
  names(C_new) <- object$C_names
  
  # Extract X and C columns
  # X_new <- as.matrix(newdata_df[, object$X_names, drop = FALSE])
  # C_new <- newdata_df[, object$C_names, drop = FALSE]
  
  # 1. Predict the conditional mean vector of X
  if (object$method == "linear") {
    mu_hat <- predict(object$inner_fit, newdata = C_new)
    mu_hat <- as.matrix(mu_hat)
    
  } else if (object$method == "SuperLearner") {
    # Initialize an empty matrix for the predictions
    mu_hat <- matrix(NA, nrow = nrow(C_new), ncol = object$p)
    
    # Predict from each univariate SuperLearner model
    for (j in 1:object$p) {
      mu_hat[, j] <- predict(object$inner_fit[[j]], newdata = C_new)$pred
    }
  }
  
  # 2. Evaluate the Multivariate Normal Density
  # Center X values by subtracting the predicted conditional means
  centered_x <- X_new - mu_hat
  
  raw_f_hat <- mvtnorm::dmvnorm(x = centered_x, 
                                mean = rep(0, object$p), 
                                sigma = object$sigma_hat)
  
  # 3. Apply flooring check
  f_hat_safe <- pmax(raw_f_hat, delta_n)
  
  return(f_hat_safe)
}