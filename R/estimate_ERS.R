#' Estimate Causal Mean using RA, IPW, or DR (Self-Normalized)
#'
#' @param Y Numeric vector of observed outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param x_eval Matrix or data frame of target treatment values to evaluate (m x p).
#' @param estimator String indicating the estimator to use: "DR", "RA", or "IPW". Defaults to "DR".
#' @param out_model An S3 object of class "outcome_model" for the outcome regression. Required for RA and DR.
#' @param gps_model An S3 object of class "gps_model" for the GPS. Required for IPW and DR.
#' @param h Optional numeric vector of bandwidths (length p). Defaults to rule-of-thumb.
#' @param c_multiplier Numeric scalar for rule-of-thumb bandwidth calculation. Defaults to 1.25.
#' @param delta_n Small positive threshold to floor propensity scores. Defaults to 1e-4.
#' @return A numeric vector of estimated causal means for each target row in x_eval.

estimate_ERS <- function(Y, X, C, x_eval = NULL, 
                                  estimator = c("DR", "RA", "IPW"), 
                                  out_model = NULL, 
                                  gps_model = NULL, 
                                  h = NULL, 
                                  c_multiplier = 1.25,
                                  delta_n = 1e-4) {
    
    estimator <- match.arg(estimator)
    
    X_df <- as.data.frame(X)
    C_df <- as.data.frame(C)
    
    n <- length(Y)
    p <- ncol(X_df)
    
    # ---------------------------------------------------------
    # Input Validation & Setup
    # ---------------------------------------------------------
    if (estimator %in% c("RA", "DR") && is.null(out_model)) {
      stop("An 'out_model' object is required for the RA and DR estimators.")
    }
    if (estimator %in% c("IPW", "DR") && is.null(gps_model)) {
      stop("A 'gps_model' object is required for the IPW and DR estimators.")
    }
    
    # Default evaluation points are the observed data
    if(is.null(x_eval)){
      x_eval <- X
    }
    x_eval_df <- as.data.frame(x_eval)
    if (ncol(x_eval_df) != p) {
      stop("x_eval must have the same number of columns as the training X.")
    }
    m <- nrow(x_eval_df)
    
    if (!is.null(out_model)) colnames(C_df) <- out_model$C_names
    
    # ---------------------------------------------------------
    # Pre-computations (Outside the evaluation loop!)
    # ---------------------------------------------------------
    
    if (estimator %in% c("IPW", "DR")) {
      # 1. Bandwidth calculations
      if (is.null(h)) {
        h <- c_multiplier * apply(X_df, 2, sd) * (n^(-0.2))
      }
      if (length(h) == 1 && p > 1) h <- rep(h, p)
      
      # Notice we no longer need h_prod because it cancels out in the self-normalization!
      
      # 2. Predict pi_hat: the density of the OBSERVED treatment given confounders
      X_obs_formatted <- X_df
      if (!is.null(gps_model$X_names)) colnames(X_obs_formatted) <- gps_model$X_names
      
      df_observed <- cbind(X_obs_formatted, C_df)
      pi_hat <- predict(gps_model, newdata = df_observed)
      pi_hat <- pmax(pi_hat, delta_n) # Flooring check
      
      # 3. SPEED OPTIMIZATION: Pre-compute the inverse propensity terms
      inv_pi <- 1 / pi_hat
      ipw_Y_weighted <- Y * inv_pi
    }
    
    # ---------------------------------------------------------
    # Main Loop over Target Treatment Values
    # ---------------------------------------------------------
    results <- vapply(1:m, function(i) {
      x_target <- as.numeric(x_eval_df[i, ])
      
      # --- Regression Adjustment Component ---
      if (estimator %in% c("RA", "DR")) {
        Xi_rep <- x_eval_df[rep(i, n), , drop = FALSE]
        colnames(Xi_rep) <- out_model$X_names
        df_new <- cbind(Xi_rep, C_df)
        
        m_hat <- predict(out_model, newdata = df_new)
      }
      
      if (estimator == "RA") {
        return(mean(m_hat))
      }
      
      # --- Kernel Distance Component ---
      K_weights <- rep(1, n)
      for (j in 1:p) {
        K_weights <- K_weights * dnorm((X_df[, j] - x_target[j]) / h[j])
      }
      
      # --- IPW Estimator (Self-Normalized) ---
      if (estimator == "IPW") {
        num <- sum(K_weights * ipw_Y_weighted)
        den <- sum(K_weights * inv_pi)
        
        # Safety check: if target is completely outside observed support
        if (den < 1e-12) return(NA_real_) 
        
        return(num / den)
      }
      
      # --- DR Estimator (Self-Normalized) ---
      if (estimator == "DR") {
        ra_est <- mean(m_hat)
        
        num <- sum(K_weights * inv_pi * (Y - m_hat))
        den <- sum(K_weights * inv_pi)
        
        # Safety check: fallback to RA if kernel weights sum to effectively zero
        if (den < 1e-12) return(ra_est)
        
        return(ra_est + (num / den))
      }
      
    }, numeric(1L))
    
    names(results) <- paste("x", 1:m, sep = "_")
    return(results)
  }