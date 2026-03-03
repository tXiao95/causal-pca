library(testthat)
library(MASS)
library(mvtnorm)

# (Insert your estimate_pseudo_outcomes function here)
source("R/estimate_pseudo_outcomes.R")
source("R/nuisance_outcome_regression.R")
source("R/nuisance_gps.R")
# ---------------------------------------------------------
# Define Lightweight Mock S3 Models
# ---------------------------------------------------------

# Mock Outcome Model: Predicts X1 + C1
mock_out_model <- function(X_names, C_names) {
  obj <- list(X_names = X_names, C_names = C_names)
  class(obj) <- "mock_out_model"
  return(obj)
}
predict.mock_out_model <- function(object, newdata, ...) {
  # Simply adds the first treatment and first confounder
  return(newdata[[object$X_names[1]]] + newdata[[object$C_names[1]]])
}

# Mock GPS Model: Predicts a constant density 
mock_gps_model <- function(X_names, C_names, dens_val = 0.5) {
  obj <- list(X_names = X_names, C_names = C_names, dens_val = dens_val)
  class(obj) <- "mock_gps_model"
  return(obj)
}
predict.mock_gps_model <- function(object, newdata, ...) {
  return(rep(object$dens_val, nrow(newdata)))
}

# ---------------------------------------------------------
# Test Suite
# ---------------------------------------------------------

test_that("estimate_pseudo_outcomes computes the exact Kennedy formula", {
  
  # 1. Setup Data
  set.seed(42)
  n <- 50
  p <- 2
  q <- 2
  
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- c("X1", "X2")
  
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- c("C1", "C2")
  
  Y <- rnorm(n)
  
  # Initialize Mock Models
  out_mod <- mock_out_model(colnames(X), colnames(C))
  gps_mod <- mock_gps_model(colnames(X), colnames(C), dens_val = 0.5)
  
  # 2. Run the Function
  computed_xi <- estimate_pseudo_outcomes(Y = Y, 
                                          X = X, 
                                          C = C, 
                                          out_model = out_mod, 
                                          gps_model = gps_mod)
  
  # 3. Calculate the Exact Mathematical Expectation
  # Based on the mock models: m(X,C) = X1 + C1, and pi(X|C) = 0.5
  # xi_j = [(Y_j - (X_j1 + C_j1)) / 0.5] * 0.5 + (X_j1 + mean(C1))
  # xi_j = Y_j - C_j1 + mean(C1)
  
  expected_xi <- Y - C[, "C1"] + mean(C[, "C1"])
  
  # ---------------------------------------------------------
  # Assertions
  # ---------------------------------------------------------
  
  # A. Structure and Types
  expect_type(computed_xi, "double")
  expect_length(computed_xi, n)
  expect_null(dim(computed_xi)) # Ensures it is a flat vector, not a matrix
  
  # B. Mathematical Accuracy
  # Tolerance is set very tight because it should be an exact algebraic match
  expect_equal(computed_xi, expected_xi, tolerance = 1e-10)
  
  # C. Test Propensity Score Flooring (delta_n)
  # If the GPS predicts 0, the denominator should hit the delta_n floor to prevent Inf
  gps_zero <- mock_gps_model(colnames(X), colnames(C), dens_val = 0)
  
  computed_xi_floored <- estimate_pseudo_outcomes(Y = Y, 
                                                  X = X, 
                                                  C = C, 
                                                  out_model = out_mod, 
                                                  gps_model = gps_zero, 
                                                  delta_n = 1e-5)
  
  # Manual calculation with floored denominator (0.5 becomes 1e-5, but mean_pi is 0)
  # Since pi_grid is identically 0, mean_pi = 0. 
  # xi_j = [(Y - m) / 1e-5] * 0 + mean_m = mean_m
  expected_floored <- X[, "X1"] + mean(C[, "C1"])
  
  expect_equal(computed_xi_floored, expected_floored, tolerance = 1e-10)
  
  # D. Missing Column Names Fallback
  # Strip names from X to test the inheritance block
  X_noname <- X
  colnames(X_noname) <- NULL
  
  computed_xi_noname <- estimate_pseudo_outcomes(Y = Y, 
                                                 X = X_noname, 
                                                 C = C, 
                                                 out_model = out_mod, 
                                                 gps_model = gps_mod)
  
  expect_equal(computed_xi_noname, expected_xi, tolerance = 1e-10)
})

test_that("estimate_pseudo_outcomes integrates flawlessly with actual nuisance models", {
  
  # ---------------------------------------------------------
  # 1. Setup Data
  # ---------------------------------------------------------
  set.seed(123)
  n <- 100
  p <- 2
  q <- 2
  
  # Confounders
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- c("C1", "C2")
  
  # Treatments (Multivariate Normal given C)
  beta_X <- matrix(c(0.5, -0.2, 
                     -0.8,  0.4), nrow = q, ncol = p)
  X <- C %*% beta_X + rmvnorm(n, sigma = diag(p))
  colnames(X) <- c("X1", "X2")
  
  # Outcome (Y = C1 - C2 + X1 - X2 + noise)
  Y <- as.numeric(C %*% c(0.5, -0.5) + X %*% c(1, -1) + rnorm(n, sd = 0.5))
  
  # ---------------------------------------------------------
  # 2. Fit Actual Nuisance Models
  # ---------------------------------------------------------
  
  # Fit Outcome Regression E[Y | X, C]
  # Using SL.glm to keep the test instantaneous
  suppressWarnings({
    out_mod <- outcome_model(Y = Y, 
                             X = X, 
                             C = C, 
                             mu_fitter = SL_outcome_fitter, 
                             SL.library = "SL.glm")
  })
  
  # Fit Propensity Score Density f(X | C) using the new gps_model gateway
  gps_mod <- gps_model(X = X, 
                       C = C, 
                       pi_fitter = mvn_fitter, 
                       method = "linear")
  
  # ---------------------------------------------------------
  # 3. Compute Pseudo-Outcomes
  # ---------------------------------------------------------
  
  pseudo_outcomes <- estimate_pseudo_outcomes(Y = Y, 
                                              X = X, 
                                              C = C, 
                                              out_model = out_mod, 
                                              gps_model = gps_mod)
  
  # ---------------------------------------------------------
  # 4. Assertions
  # ---------------------------------------------------------
  
  # A. Structural Integrity
  expect_type(pseudo_outcomes, "double")
  expect_length(pseudo_outcomes, n)
  expect_null(dim(pseudo_outcomes)) # Must be a flat numeric vector
  
  # B. Safety Checks
  # Ensure no NAs, NaNs, or Infs were produced by exploding inverse probability weights
  expect_true(all(!is.na(pseudo_outcomes)))
  expect_true(all(is.finite(pseudo_outcomes)))
  
  # C. Statistical Sanity Check
  # By the Law of Total Expectation, the mean of the pseudo-outcomes should 
  # be roughly in the same neighborhood as the mean of the observed outcomes 
  # (though exact equality doesn't hold in finite samples due to the grid loop).
  mean_Y <- mean(Y)
  mean_po <- mean(pseudo_outcomes)
  
  # They should be on the same scale (within a reasonable finite-sample bound)
  expect_lt(abs(mean_po - mean_Y), 1.0)
})
