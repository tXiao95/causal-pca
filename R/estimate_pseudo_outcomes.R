#' Compute Doubly Robust Pseudo-Outcomes from Kennedy et al. (2017) JRSS-B
#'
#' @param Y Numeric vector of observed outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param out_model An S3 object of class "outcome_model".
#' @param gps_model An S3 object representing a global conditional density model.
#' @param delta_n Small positive threshold to floor propensity scores. Defaults to 1e-4.
#' @return A numeric vector of pseudo-outcomes (length n).

estimate_pseudo_outcomes <- function(Y, X, C, out_model, gps_model, delta_n = 1e-4) {
  
  # ---------------------------------------------------------
  # Input Validation & Safe Column Naming
  # ---------------------------------------------------------
  # Capture original names BEFORE coercion
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL
  
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  n <- length(Y)
  
  # Safely apply names only if the ORIGINAL input lacked them
  if (is.null(orig_X_names) && !is.null(out_model$X_names)) colnames(X_df) <- out_model$X_names
  if (is.null(orig_C_names) && !is.null(out_model$C_names)) colnames(C_df) <- out_model$C_names
  
  # ---------------------------------------------------------
  # Pre-computations (Batch predicting observed data)
  # ---------------------------------------------------------
  df_observed <- cbind(X_df, C_df)
  
  # Pre-calculate the "diagonal" elements: m(X_j, C_j) and pi(X_j | C_j)
  # Doing this outside the loop saves overhead and is much safer
  m_obs  <- predict(out_model, newdata = df_observed)
  pi_obs <- predict(gps_model, newdata = df_observed)
  pi_obs <- pmax(pi_obs, delta_n) # Apply safety flooring to the denominator
  
  # ---------------------------------------------------------
  # Main Loop over Individuals
  # ---------------------------------------------------------
  pseudo_outcomes <- vapply(1:n, function(j) {
    
    # Create a grid where individual j's treatment is repeated n times,
    # paired with EVERY individual's confounders (C_1 to C_n)
    X_j_rep <- X_df[rep(j, n), , drop = FALSE]
    df_grid <- cbind(X_j_rep, C_df)
    
    # Predict m(X_j, C_i) and pi(X_j | C_i) across all i = 1...n
    m_grid  <- predict(out_model, newdata = df_grid)
    pi_grid <- predict(gps_model, newdata = df_grid)
    
    # Calculate the empirical expectations (marginalized over C)
    mean_pi <- mean(pi_grid)
    mean_m  <- mean(m_grid)
    
    # Assemble the pseudo-outcome using the pre-computed observed values
    xi_j <- ((Y[j] - m_obs[j]) / pi_obs[j]) * mean_pi + mean_m
    
    return(xi_j)
    
  }, numeric(1L))
  
  return(pseudo_outcomes)
}