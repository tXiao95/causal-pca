library(SuperLearner)
library(parallel)

#' Fit a Global Outcome Regression Object E[Y | X, C]
#' 
#' @param Y Numeric vector of outcomes.
#' @param X Numeric matrix or data frame of observed treatments.
#' @param C Numeric matrix or data frame of observed confounders.
#' @param mu_fitter Function(Y, XC_df) that trains and returns a model.
#' @return An S3 object of class "outcome_model".

outcome_model <- function(Y, X, C, mu_fitter, ...) {
  n <- length(Y)
  p <- ncol(X)
  q <- ncol(C)
  
  # Standardize column names for safe prediction later
  C_df <- as.data.frame(C); colnames(C_df) <- paste0("C", 1:q)
  X_df <- as.data.frame(X); colnames(X_df) <- paste0("X", 1:p)
  
  df <- cbind(X_df, C_df)
  
  # Fit the model
  inner_fit <- mu_fitter(Y, df, ...)
  
  res <- list(
    inner_fit = inner_fit,
    p = p,
    q = q,
    X_names = colnames(X_df),
    C_names = colnames(C_df),
    X_observed = X_df,
    C_observed = C_df
  )
  
  class(res) <- "outcome_model"
  return(res)
}

#' Predict Method for Outcome Model
predict.outcome_model <- function(object, newdata, ...) {
  preds <- predict(object$inner_fit, newdata = newdata, ...)
  if (is.list(preds) && "pred" %in% names(preds)) {
    return(as.numeric(preds$pred))
  }
  return(as.numeric(preds))
}

#' G-Computation Estimator   sum E(Y | X = x_new, C = C_new)
#' 
#' @param out_model An S3 object of class "outcome_model".
#' @param C Numeric matrix or data frame of observed confounders.
#' @param X_new Optional numeric matrix or data frame of new treatment values to evaluate. Defaults to NULL.
#' @return A numeric vector of marginalized expected outcomes.

gcomp <- function(outcome_model, C = NULL, X_new = NULL) {
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

# Define the wrapper
SL_outcome_fitter <- function(Y, XC_df, ...) {
  SuperLearner::SuperLearner(Y = Y, X = XC_df, family = gaussian(), ...)
}