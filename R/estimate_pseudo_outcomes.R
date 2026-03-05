library(data.table)

#' Compute Doubly Robust Pseudo-Outcomes from Kennedy et al. (2017) JRSS-B
#'
#' @param Y Numeric vector of observed outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param out_model An S3 object of class "outcome_model".
#' @param gps_model An S3 object representing a global conditional density model.
#' @param delta_n Small positive threshold to floor propensity scores. Defaults to 1e-4.
#' @return A numeric vector of pseudo-outcomes (length n).

# estimate_pseudo_outcomes_old <- function(Y, X, C, out_model, gps_model, delta_n = 1e-16) {
#   
#   # ---------------------------------------------------------
#   # Input Validation & Safe Column Naming
#   # ---------------------------------------------------------
#   # Capture original names BEFORE coercion
#   orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
#   orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL
#   
#   X_df <- as.data.frame(X)
#   C_df <- as.data.frame(C)
#   n <- length(Y)
#   
#   # Safely apply names only if the ORIGINAL input lacked them
#   if (is.null(orig_X_names) && !is.null(out_model$X_names)) colnames(X_df) <- out_model$X_names
#   if (is.null(orig_C_names) && !is.null(out_model$C_names)) colnames(C_df) <- out_model$C_names
#   
#   # ---------------------------------------------------------
#   # Pre-computations (Batch predicting observed data)
#   # ---------------------------------------------------------
#   df_observed <- cbind(X_df, C_df)
#   
#   # Pre-calculate the "diagonal" elements: m(X_j, C_j) and pi(X_j | C_j)
#   # Doing this outside the loop saves overhead and is much safer
#   m_obs  <- predict(out_model, newdata = df_observed)
#   pi_obs <- predict(gps_model, newdata = df_observed)
#   pi_obs <- pmax(pi_obs, delta_n) # Apply safety flooring to the denominator
#   
#   # ---------------------------------------------------------
#   # Main Loop over Individuals
#   # ---------------------------------------------------------
#   pseudo_outcomes <- vapply(1:n, function(j) {
#     message(j)
#     # Create a grid where individual j's treatment is repeated n times,
#     # paired with EVERY individual's confounders (C_1 to C_n)
#     X_j_rep <- X_df[rep(j, n), , drop = FALSE]
#     df_grid <- cbind(X_j_rep, C_df)
#     
#     # Predict m(X_j, C_i) and pi(X_j | C_i) across all i = 1...n
#     m_grid  <- predict(out_model, newdata = df_grid)
#     pi_grid <- predict(gps_model, newdata = df_grid)
#     
#     # Calculate the empirical expectations (marginalized over C)
#     mean_pi <- mean(pi_grid)
#     mean_m  <- mean(m_grid)
#     
#     # Assemble the pseudo-outcome using the pre-computed observed values
#     xi_j <- ((Y[j] - m_obs[j]) / pi_obs[j]) * mean_pi + mean_m
#     
#     return(xi_j)
#     
#   }, numeric(1L))
#   
#   return(pseudo_outcomes)
# }

#' #' Compute Doubly Robust Pseudo-Outcomes from Kennedy et al. (2017) JRSS-B
#' #'
#' #' @param Y Numeric vector of observed outcomes (length n).
#' #' @param X Numeric matrix or data frame of observed treatments (n x p).
#' #' @param C Numeric matrix or data frame of observed confounders (n x q).
#' #' @param out_model An S3 object of class "outcome_model".
#' #' @param gps_model An S3 object representing a global conditional density model.
#' #' @param delta_n Small positive threshold to floor propensity scores. Defaults to 1e-16.
#' #' @param max_predict_size Maximum number of rows to predict at once to balance RAM and speed.
#' #' @return A numeric vector of pseudo-outcomes (length n).
#' 
#' estimate_pseudo_outcomes <- function(Y, X, C, out_model, gps_model, delta_n = 1e-16, max_predict_size = 500000) {
#'   
#'   # ---------------------------------------------------------
#'   # Input Validation & Safe Column Naming
#'   # ---------------------------------------------------------
#'   orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
#'   orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL
#'   
#'   X_df <- as.data.frame(X)
#'   C_df <- as.data.frame(C)
#'   n <- length(Y)
#'   
#'   if (is.null(orig_X_names) && !is.null(out_model$X_names)) colnames(X_df) <- out_model$X_names
#'   if (is.null(orig_C_names) && !is.null(out_model$C_names)) colnames(C_df) <- out_model$C_names
#'   
#'   # ---------------------------------------------------------
#'   # Pre-computations (Batch predicting observed data)
#'   # ---------------------------------------------------------
#'   df_observed <- cbind(X_df, C_df)
#'   
#'   m_obs  <- predict(out_model, newdata = df_observed)
#'   pi_obs <- predict(gps_model, newdata = df_observed)
#'   pi_obs <- pmax(pi_obs, delta_n) 
#'   
#'   # ---------------------------------------------------------
#'   # Chunked Grid Prediction (Balances Speed and Memory)
#'   # ---------------------------------------------------------
#'   # Determine how many individuals to process in one batch
#'   chunk_size <- max(1, floor(max_predict_size / n))
#'   
#'   # Split individuals 1:n into chunks
#'   indices <- 1:n
#'   chunks <- split(indices, ceiling(seq_along(indices) / chunk_size))
#'   
#'   pseudo_outcomes <- numeric(n)
#'   
#'   for (chk in chunks) {
#'     message("Chunk ", chk, " out of ", length(chunks), " chunks")
#'     chunk_len <- length(chk)
#'     
#'     # 1. Expand X and C for the current chunk
#'     # X repeats each individual's row 'n' times
#'     idx_X_rep <- rep(chk, each = n)
#'     X_chunk <- X_df[idx_X_rep, , drop = FALSE]
#'     
#'     # C tiles the entire dataset 'chunk_len' times
#'     idx_C_rep <- rep(1:n, times = chunk_len)
#'     C_chunk <- C_df[idx_C_rep, , drop = FALSE]
#'     
#'     df_grid <- cbind(X_chunk, C_chunk)
#'     
#'     # 2. Batch Predict
#'     m_grid  <- predict(out_model, newdata = df_grid)
#'     pi_grid <- predict(gps_model, newdata = df_grid)
#'     
#'     # 3. Reshape and calculate means instantly
#'     # Converting the 1D prediction vector into an (n x chunk_len) matrix 
#'     # allows colMeans to calculate the marginalized mean for each individual instantly.
#'     mean_m  <- colMeans(matrix(m_grid, nrow = n, ncol = chunk_len))
#'     mean_pi <- colMeans(matrix(pi_grid, nrow = n, ncol = chunk_len))
#'     
#'     # 4. Assemble the pseudo-outcomes for this chunk
#'     xi_chunk <- ((Y[chk] - m_obs[chk]) / pi_obs[chk]) * mean_pi + mean_m
#'     
#'     pseudo_outcomes[chk] <- xi_chunk
#'   }
#'   
#'   return(pseudo_outcomes)
#' }

#' Compute Doubly Robust Pseudo-Outcomes from Kennedy et al. (2017) JRSS-B
#'
#' @param Y Numeric vector of observed outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param out_model An S3 object of class "outcome_model".
#' @param gps_model An S3 object representing a global conditional density model.
#' @param delta_n Small positive threshold to floor propensity scores. Defaults to 1e-16.
#' @return A numeric vector of pseudo-outcomes (length n).

estimate_pseudo_outcomes <- function(Y, X, C, out_model, gps_model, delta_n = 1e-16) {
  
  # ---------------------------------------------------------
  # Input Validation & Safe Column Naming
  # ---------------------------------------------------------
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL
  
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  n <- length(Y)
  
  if (is.null(orig_X_names) && !is.null(out_model$X_names)) colnames(X_df) <- out_model$X_names
  if (is.null(orig_C_names) && !is.null(out_model$C_names)) colnames(C_df) <- out_model$C_names
  
  # ---------------------------------------------------------
  # Pre-computations (Batch predicting observed data)
  # ---------------------------------------------------------
  df_observed <- cbind(X_df, C_df)
  
  m_obs  <- predict(out_model, newdata = df_observed)
  pi_obs <- predict(gps_model, newdata = df_observed)
  pi_obs <- pmax(pi_obs, delta_n) 
  
  # ---------------------------------------------------------
  # Zero-Copy Evaluation Grid Setup
  # ---------------------------------------------------------
  
  # 1. Base the grid entirely on C (which never changes for the marginalization)
  dt_grid <- as.data.table(C_df)
  
  # 2. Pre-allocate the X columns with dummy values (0.0)
  x_names <- colnames(X_df)
  for (x_col in x_names) {
    dt_grid[, (x_col) := 0.0]
  }
  
  # 3. Force the exact column order the models expect
  setcolorder(dt_grid, c(x_names, colnames(C_df)))
  
  # ---------------------------------------------------------
  # Main Loop over Individuals (Hyper-Optimized)
  # ---------------------------------------------------------
  pseudo_outcomes <- vapply(1:n, function(j) {
    message(j)
    
    # ZERO-COPY IN-PLACE UPDATE:
    # Instead of `rep()` and `cbind()`, we instantly overwrite the memory 
    # of the X columns with the values for individual 'j'.
    for (x_col in x_names) {
      data.table::set(dt_grid, j = x_col, value = X_df[[x_col]][j])
    }
    
    # Predict 
    m_grid  <- predict(out_model, newdata = dt_grid)
    pi_grid <- predict(gps_model, newdata = dt_grid)
    
    # Calculate the empirical expectations
    mean_pi <- mean(pi_grid)
    mean_m  <- mean(m_grid)
    
    # Assemble the pseudo-outcome
    xi_j <- ((Y[j] - m_obs[j]) / pi_obs[j]) * mean_pi + mean_m
    
    return(xi_j)
    
  }, numeric(1L))
  
  return(pseudo_outcomes)
}
