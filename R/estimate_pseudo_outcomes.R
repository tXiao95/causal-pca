#' Compute Doubly Robust Pseudo-Outcomes
#'
#' @param Y Numeric vector of observed outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param out_model An S3 object of class "outcome_model".
#' @param gps_model An S3 object representing a global conditional density model.
#' @param delta_n Small positive threshold to floor propensity scores. Defaults to 1e-4.
#' @return A numeric vector of pseudo-outcomes (length n).

estimate_pseudo_outcomes <- function(Y, X, C, out_model, gps_model, delta_n = 1e-4) {
  n <- length(Y)
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  
  # Ensure C and X match the names used during model training
  colnames(C_df) <- out_model$C_names
  colnames(X_df) <- out_model$X_names
  
  # Main loop: calculate the pseudo-outcome for each individual j
  pseudo_outcomes <- vapply(1:n, function(j) {
    
    # Create a grid where individual j's treatment is repeated n times,
    # paired with EVERY individual's confounders (C_1 to C_n)
    X_j_rep <- X_df[rep(j, n), , drop = FALSE]
    df_grid <- cbind(X_j_rep, C_df)
    
    # Predict m(X_j, C_i) and pi(X_j | C_i) across all i = 1...n
    m_grid <- predict(out_model, newdata = df_grid)
    pi_grid <- predict(gps_model, newdata = df_grid)
    
    # We also need m(X_j, C_j) and pi(X_j | C_j). 
    # Because the j-th row of df_grid pairs X_j exactly with C_j, 
    # these values are simply the j-th elements of the grids we just predicted.
    m_jj <- m_grid[j]
    pi_jj <- pmax(pi_grid[j], delta_n) # Apply safety flooring
    
    # Calculate the empirical expectations (means) over i=1...n
    mean_pi <- mean(pi_grid)
    mean_m  <- mean(m_grid)
    
    # Assemble the pseudo-outcome xi(O_j)
    xi_j <- ((Y[j] - m_jj) / pi_jj) * mean_pi + mean_m
    
    return(xi_j)
    
  }, numeric(1L))
  
  return(pseudo_outcomes)
}