library(testthat)
library(here)
library(mvtnorm)

# (Assume gps_model, predict.gps_model, mvn_fitter, and predict.mvn_inner are loaded)
source(here("R/nuisance_gps.R"))

test_that("gps_model and mvn_fitter integrate correctly and estimate accurate densities", {
  
  # ---------------------------------------------------------
  # 1. Simulate Data with a Known MVN GPS
  # ---------------------------------------------------------
  set.seed(42)
  n <- 500
  p <- 2
  q <- 2
  
  # Confounders
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- c("C1", "C2")
  
  # Treatments X | C ~ MVN(C %*% beta, Sigma)
  beta_true <- matrix(c(0.5, -0.2, 
                        -0.8,  0.4), nrow = q, ncol = p)
  
  # Create a positive definite true covariance matrix
  Sigma_true <- matrix(c(1.0, 0.3, 
                         0.3, 1.0), nrow = p, ncol = p)
  
  mu_true <- C %*% beta_true
  X <- mu_true + rmvnorm(n, sigma = Sigma_true)
  colnames(X) <- c("X1", "X2")
  
  # ---------------------------------------------------------
  # 2. Fit the Abstracted GPS Model
  # ---------------------------------------------------------
  # Pass the inner fitter function and its arguments (method = "linear")
  fit <- gps_model(X = X, C = C, pi_fitter = mvn_fitter, method = "linear")
  
  # --- Assertions for Object Structure ---
  expect_s3_class(fit, "gps_model")
  expect_s3_class(fit$inner_fit, "mvn_inner")
  expect_equal(fit$X_names, c("X1", "X2"))
  expect_equal(fit$C_names, c("C1", "C2"))
  
  # ---------------------------------------------------------
  # 3. Test Predictions and Statistical Accuracy
  # ---------------------------------------------------------
  newdata <- cbind(as.data.frame(X), as.data.frame(C))
  predicted_densities <- predict(fit, newdata = newdata)
  
  # Calculate theoretical true densities
  true_densities <- numeric(n)
  for(i in 1:n) {
    true_densities[i] <- dmvnorm(x = X[i, ], mean = mu_true[i, ], sigma = Sigma_true)
  }
  
  # --- Assertions for Accuracy ---
  expect_type(predicted_densities, "double")
  expect_length(predicted_densities, n)
  
  # The correlation between estimated and true densities should be highly positive
  expect_gt(cor(predicted_densities, true_densities), 0.98)
  
  # Mean Absolute Error should be small
  expect_lt(mean(abs(predicted_densities - true_densities)), 0.05)
  
  # ---------------------------------------------------------
  # 4. Test API Safety Contracts
  # ---------------------------------------------------------
  
  # A. Missing Columns in Predict
  newdata_missing <- newdata[, c("X1", "C1", "C2")] # Dropped X2
  expect_error(
    predict(fit, newdata = newdata_missing),
    "The following required columns are missing from 'newdata': X2"
  )
  
  # B. Scrambled Columns in Predict (Should be forcefully reordered)
  newdata_scrambled <- newdata[, c("C2", "X1", "C1", "X2")]
  preds_scrambled <- predict(fit, newdata = newdata_scrambled)
  expect_equal(predicted_densities, preds_scrambled)
  
  # C. Overlapping Column Names in Setup
  C_bad <- C
  colnames(C_bad) <- c("X1", "C2") # Force a name overlap
  expect_error(
    gps_model(X = X, C = C_bad, pi_fitter = mvn_fitter, method = "linear"),
    "Overlapping column names detected: X1"
  )
  
  # D. Missing Column Names Fallback
  X_noname <- X
  colnames(X_noname) <- NULL
  fit_noname <- gps_model(X = X_noname, C = C, pi_fitter = mvn_fitter, method = "linear")
  expect_equal(fit_noname$X_names, c("X1", "X2"))
})
