#' Causal Sufficient Dimension Reduction via csPCA (causally sufficient PCA)
#'
#' Implements a two-step procedure for causal sufficient dimension reduction:
#' (1) estimate the causal mean function \eqn{\mu(X)} using a user-specified method,
#' and (2) apply MAVE (minimum average variance estimation) to find low-dimensional
#' projections of the covariates that best capture variation in \eqn{\mu(X)}.
#'
#' @param Y A numeric vector of outcomes.
#' @param X A numeric matrix or data frame of covariates (rows = observations, columns = predictors).
#' @param C A numeric matrix or data frame of confounders (rows aligned with \code{X}).
#' @param mu_fun A function to estimate the causal mean \eqn{\mu(X)} given \code{Y}, \code{X}, and \code{C}.
#'   Defaults to \code{gcomp}, but may be replaced by any user-specified function with the same signature.
#' @param mu_args A named list of additional arguments passed to \code{mu_fun}.
#' @param mave_args A named list of additional arguments passed to \code{MAVE::mave}.
#'
#' @details
#' The function estimates \eqn{\mu(X)} by calling \code{mu_fun(Y, X, C, ...)}.
#' The fitted values are then treated as a univariate response in a mean-MAVE
#' regression of \eqn{\mu(X)} on \code{X}. This identifies low-dimensional
#' sufficient directions in the covariates for predicting the causal mean.
#'
#' By default, the procedure uses regression adjustment via \code{gcomp} for \eqn{\mu(X)},
#' and mean-MAVE from the \pkg{MAVE} package for sufficient dimension reduction.
#'
#' @return An object of class \code{"csPCA_fit"}, which is a list with components:
#' \describe{
#'   \item{mave}{The fitted MAVE object returned by \code{MAVE::mave}.}
#'   \item{mu_X}{The estimated causal mean vector \eqn{\mu(X)}.}
#' }
csPCA <- function(Y, X, C,
                  L = 1, folds = NULL,
                  mu_fun = gcomp,          # function to estimate mu(X)
                  mu_args = list(),        # args passed to mu_fun
                  mave_args = list()) {    # args passed to MAVE::mave
  
  stopifnot(length(Y) == nrow(X), nrow(X) == nrow(C))
  n <- nrow(X)
  # Cross-fitting setup -----------------------------------------------------
  if (is.null(folds)) 
    folds <- sample(rep(1:L, length.out = n))
  stopifnot(length(folds) == n, all(folds %in% 1:L))
  
  
  if(L > 1){
    message("Cross-fitting with ", L, " folds")
    mu_X <- numeric(n)
    for (l in 1:L) {
      message("Fold ", l)
      idx_tr <- which(folds != l)
      idx_te <- which(folds == l)
      
      # fit on training fold -> get predictor closure
      pred_mu <- do.call( mu_fun, c(list(Y = Y[idx_tr],
                                        X = X[idx_tr, , drop = FALSE], 
                                        C = C[idx_tr, , drop = FALSE], 
                                        X.new = X[idx_te, , drop = FALSE]), mu_args) )
      
      # predict mu(X_i) for test fold, integrating over C from training fold
      mu_X[idx_te] <- pred_mu
    }
  } else{
    # 1) Estimate causal mean mu(X)
    message("Super Learner")
    mu_X <- do.call(mu_fun, c(list(Y = Y, X = X, C = C), mu_args))
  }
  
  # 2) Prepare data for MAVE
  df <- data.frame(mu_X = as.numeric(mu_X), as.data.frame(X), check.names = FALSE)
  
  # 3) Build MAVE args with sensible defaults, allow user overrides
  message("MAVE")
  base_args <- list(formula = mu_X ~ ., data = df, method = "meanMAVE")
  fit <- do.call(MAVE::mave, utils::modifyList(base_args, mave_args))
  
  # 4) Postprocessing: ensure MAVE outputs for each 'd' are orthonormal.
  
  structure(list(mave = fit, mu_X = mu_X), class = "csPCA_fit")
}
