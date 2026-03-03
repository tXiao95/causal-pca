library(testthat)

source("R/compute_new_response_and_exposure.R")
source("R/estimate_ERS.R")
source("R/estimate_pseudo_outcomes.R")
source("R/estimate_residualized_pair.R")

source("R/nuisance_gps.R")
source("R/nuisance_outcome_regression.R")
source("R/nuisance_outcome_C.R")

# ---------------------------------------------------------
# Define Lightweight Fitters for Fast Testing
# ---------------------------------------------------------
fast_lm_fitter <- function(target_var, predictors_df, ...) {
  # Dynamically bind target and predictors for standard lm()
  df <- cbind(target_var = target_var, predictors_df)
  lm(target_var ~ ., data = df)
}

# ---------------------------------------------------------
# Test Suite
# ---------------------------------------------------------
test_that("compute_new_response_and_exposure correctly cross-fits all methods", {
  
  # 1. Setup Data
  set.seed(42)
  n <- 500 # Multiple of L=5 for even folds
  p <- 2
  q <- 2
  
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- c("X1", "X2")
  
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- c("C1", "C2")
  
  Y <- as.numeric(C %*% c(0.5, -0.5) + X %*% c(1, -1) + rnorm(n))
  
  # =======================================================
  # Test 1: Regression Adjustment (RA)
  # =======================================================
  suppressWarnings({
    res_ra <- compute_new_response_and_exposure(Y = Y, X = X, C = C, 
                                                method = "RA", 
                                                L = 5,
                                                outcome_fitter = SL_outcome_fitter, 
                                                args_outcome = list(SL.lib = c("SL.glm", "SL.gam", "SL.nnet"),
                                                                    cvControl = list(V = 2)))
  })
  
  expect_type(res_ra, "list")
  expect_length(res_ra$new_Y, n)
  expect_equal(dim(res_ra$new_X), c(n, p))
  expect_equal(res_ra$metadata$method, "RA")
  
  # X should remain unchanged for RA
  expect_equal(as.numeric(res_ra$new_X), as.numeric(X))
  
  # =======================================================
  # Test 2: Doubly Robust (DR)
  # =======================================================
  suppressWarnings({
    res_dr <- compute_new_response_and_exposure(Y = Y, X = X, C = C, 
                                                method = "DR", 
                                                L = 5,
                                                outcome_fitter = fast_lm_fitter,
                                                gps_fitter = mvn_fitter,
                                                args_ers = list(optimize_bw = FALSE)) # Passed via ... to mvn_fitter
  })
  
  expect_type(res_dr, "list")
  expect_length(res_dr$new_Y, n)
  expect_true(all(!is.na(res_dr$new_Y)))
  
  # =======================================================
  # Test 3: Pseudo-Outcome (PO)
  # =======================================================
  suppressWarnings({
    res_po <- compute_new_response_and_exposure(Y = Y, X = X, C = C, 
                                                method = "PO", 
                                                L = 5,
                                                outcome_fitter = fast_lm_fitter,
                                                gps_fitter = mvn_fitter)
  })
  
  expect_type(res_po, "list")
  expect_length(res_po$new_Y, n)
  expect_true(all(!is.na(res_po$new_Y)))
  
  # =======================================================
  # Test 4: Residualized Pair (RP)
  # =======================================================
  res_rp <- compute_new_response_and_exposure(Y = Y, X = X, C = C, 
                                              method = "RP", 
                                              L = 5,
                                              C_fitter = fast_lm_fitter)
  
  expect_type(res_rp, "list")
  expect_length(res_rp$new_Y, n)
  
  # X MUST be modified for the Residualized Pair method
  expect_false(isTRUE(all.equal(as.numeric(res_rp$new_X), as.numeric(X))))
  
  # =======================================================
  # Test 5: Seed Reproducibility
  # =======================================================
  res_rp_2 <- compute_new_response_and_exposure(Y = Y, X = X, C = C, 
                                                method = "RP", 
                                                L = 5,
                                                C_fitter = fast_lm_fitter,
                                                seed = 42) # Same default seed
  
  # The outputs should be exactly identical because the folds were identical
  expect_equal(res_rp$new_Y, res_rp_2$new_Y)
  expect_equal(res_rp$new_X, res_rp_2$new_X)
  
  # =======================================================
  # Test 6: Error Handling for Missing Fitters
  # =======================================================
  expect_error(
    compute_new_response_and_exposure(Y = Y, X = X, C = C, method = "RA", outcome_fitter = NULL),
    "A valid 'outcome_fitter' function must be provided"
  )
  expect_error(
    compute_new_response_and_exposure(Y = Y, X = X, C = C, method = "RP", C_fitter = NULL),
    "A valid 'C_fitter' function must be provided"
  )
})
