library(MAVE)

#' Causal Sufficient Dimension Reduction via MAVE (csMAVE)
#'
#' @param Y Numeric vector of outcomes (length n).
#' @param X Numeric matrix or data frame of observed treatments (n x p).
#' @param C Numeric matrix or data frame of observed confounders (n x q).
#' @param method String indicating which causal transformation to use ("RA", "DR", "PO", "RP").
#' @param args_compute_new_response List of top-level arguments for `compute_new_response_and_exposure` 
#'        (e.g., list(L = 10, outcome_fitter = SL_outcome_fitter, seed = 123)).
#' @param args_outcome List of tuning parameters for the outcome fitter.
#' @param args_C List of tuning parameters for the C fitter.
#' @param args_gps List of tuning parameters for the GPS fitter.
#' @param args_ers List of tuning parameters for the ERS estimation.
#' @param args_MAVE List of arguments to pass to `MAVE::mave` (e.g., list(method = "meanOPG")).
#' @return A list containing the MAVE fit, dimension estimation object, selected dimension, 
#'         the generated pseudo-data, and run metadata.

csMAVE <- function(Y, X, C, 
                   method = c("RA", "DR", "PO", "RP"),
                   args_compute_new_response = list(), 
                   args_outcome = list(),
                   args_C = list(),
                   args_gps = list(),
                   args_ers = list(),
                   args_MAVE = list()) {
  
  method <- match.arg(method)
  
  # ---------------------------------------------------------
  # 1. Compute the New Response and Exposure (Causal Transformation)
  # ---------------------------------------------------------
  
  # Build the base argument list mapping to the exact arguments required
  cre_base_args <- list(
    Y = Y, 
    X = X, 
    C = C, 
    method = method,
    args_outcome = args_outcome,
    args_C = args_C,
    args_gps = args_gps,
    args_ers = args_ers
  )
  
  # Merge with any user-supplied top-level args (like L, seed, fitters).
  # utils::modifyList safely overwrites defaults without causing duplicates.
  cre_final_args <- utils::modifyList(cre_base_args, args_compute_new_response)
  
  # Execute the cross-fitting pipeline
  message("Computing new response and exposure...")
  new_data_obj <- do.call(compute_new_response_and_exposure, cre_final_args)
  
  new_Y <- new_data_obj$new_Y
  new_X <- new_data_obj$new_X
  p <- ncol(new_X)
  
  # ---------------------------------------------------------
  # 2. Prepare Data for MAVE
  # ---------------------------------------------------------
  
  # Combine into a single dataframe. We name the response 'newY' for the formula
  df <- data.frame(newY = new_Y, new_X)
  
  # ---------------------------------------------------------
  # 3. Fit MAVE
  # ---------------------------------------------------------
  
  # Base arguments for MAVE
  mave_base_args <- list(formula = newY ~ ., data = df, method = "meanMAVE")
  
  # Merge with user-supplied MAVE arguments (allows user to easily override "meanOPG")
  mave_final_args <- utils::modifyList(mave_base_args, args_MAVE)
  
  message("Running MAVE...")
  fit_mave <- do.call(MAVE::mave, mave_final_args)
  
  # ---------------------------------------------------------
  # 4. Estimate Structural Dimension
  # ---------------------------------------------------------
  
  # max.dim is bounded by the number of columns in the exposure matrix
  message("Estimating the structural dimension...")
  dhat_obj <- MAVE::mave.dim(fit_mave, max.dim = p)
  
  # Extract the specific chosen dimension (MAVE typically stores this in $dim)
  d_hat <- if (!is.null(dhat_obj$dim)) dhat_obj$dim.min else NA
  
  # ---------------------------------------------------------
  # 5. Format and Return Results
  # ---------------------------------------------------------
  
  message("DONE!")
  return(list(
    mave_fit     = fit_mave,
    mave_dim_obj = dhat_obj,
    d_hat        = d_hat,
    new_data     = list(new_Y = new_Y, new_X = new_X),
    metadata     = list(
      causal_method   = method,
      mave_method     = mave_final_args$method,
      n_observations  = nrow(new_X),
      p_exposures     = p,
      cre_pipeline    = new_data_obj$metadata
    )
  ))
}