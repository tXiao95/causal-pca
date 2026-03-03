library(testthat)
library(SuperLearner)
library(mvtnorm)

# (Assume your outcome_model, mvnGPS, and estimate_ERS functions are loaded)
source("R/estimate_ERS.R")
source("R/nuisance_outcome_regression.R")
source("R/nuisance_mvnGPS.R")

test_that("estimate_ERS correctly recovers the true causal mean and handles bandwidth optimization", {
  
  # ---------------------------------------------------------
  # 1. Setup Data with a Known Causal Truth
  # ---------------------------------------------------------
  set.seed(42)
  n <- 1000 # Increased slightly to ensure kernel smoothing has enough local data
  p <- 2
  q <- 2
  
  # Confounders
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- c("C1", "C2")
  
  # Treatments (Confounded by C)
  X <- C + rmvnorm(n, sigma = diag(p))
  colnames(X) <- c("X1", "X2")
  
  # Outcome (Y = C1 + C2 + X1 + X2 + noise)
  Y <- as.numeric(C[,1] + C[,2] + X[,1] + X[,2] + rnorm(n, sd = 0.5))
  
  # THE MATHEMATICAL TRUTH:
  # Because E[C] = 0, the true causal dose-response at any x is exactly x1 + x2.
  true_ERS <- X[, 1] + X[, 2]
  
  # ---------------------------------------------------------
  # 2. Fit Nuisance Models
  # ---------------------------------------------------------
  # We use GLM/Linear to perfectly match the linear DGP, isolating the 
  # test to verify the causal integration logic, not the machine learning algorithms.
  
  suppressWarnings({
    out_mod <- outcome_model(Y = Y, X = X, C = C, 
                             mu_fitter = SL_outcome_fitter, 
                             SL.library = "SL.glm")
  })
  gps_mod <- mvnGPS(X = X, C = C, method = "linear")
  
  # ---------------------------------------------------------
  # 3. Test Estimators against the True Causal Mean
  # ---------------------------------------------------------
  
  # A. Regression Adjustment (RA)
  suppressWarnings({
    res_ra <- estimate_ERS(Y, X, C, estimator = "RA", out_model = out_mod)
  })
  # RA should be highly accurate here because the outcome model is correctly specified
  mae_ra <- mean(abs(res_ra - true_ERS))
  expect_lt(mae_ra, 0.1) 
  
  # B. Inverse Probability Weighting (IPW)
  res_ipw <- estimate_ERS(Y, X, C, estimator = "IPW", gps_model = gps_mod)
  # IPW relies on kernel smoothing, so it naturally has higher variance/error than RA in finite samples
  mae_ipw <- mean(abs(res_ipw - true_ERS))
  expect_lt(mae_ipw, 0.6)
  
  # C. Doubly Robust (DR) - Standard Bandwidth
  res_dr_vec <- estimate_ERS(Y, X, C, estimator = "DR", out_model = out_mod, 
                             gps_model = gps_mod, return_vector = TRUE)
  mae_dr <- mean(abs(res_dr_vec - true_ERS))
  # DR should stabilize the IPW kernel variance, performing better
  expect_lt(mae_dr, 0.4)
  
  # ---------------------------------------------------------
  # 4. Test Bandwidth Optimization (optimize_bw = TRUE)
  # ---------------------------------------------------------
  res_dr_opt <- estimate_ERS(Y, X, C, 
                             estimator = "DR", 
                             out_model = out_mod, 
                             gps_model = gps_mod,
                             optimize_bw = TRUE,
                             return_vector = FALSE) # Return the full list to check metadata
  
  # Check Structural Integrity
  expect_type(res_dr_opt, "list")
  expect_named(res_dr_opt, c("results", "metadata"))
  expect_true(res_dr_opt$metadata$optimized_bw)
  
  # Verify the optimized estimates are also statistically valid
  mae_dr_opt <- mean(abs(res_dr_opt$results$estimate - true_ERS))
  expect_lt(mae_dr_opt, 0.4)
  
  # Check that bandwidths were dynamically updated and saved in the output
  # They should not be NA and should be numeric
  expect_true(all(!is.na(res_dr_opt$results$h_used_1)))
  expect_type(res_dr_opt$results$h_used_1, "double")
})