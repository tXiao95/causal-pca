library(testthat)
library(SuperLearner)

# (Insert your nuisance_C_model, predict.nuisance_C_model, and SL_nuisance_fitter here)
source("R/residualized_pair.R")

test_that("nuisance_C_model fits and predicts correctly", {
  
  # 1. Setup: Simulate simple confounding data
  set.seed(123)
  n_train <- 10000
  n_test <- 500
  q <- 3
  
  C_train <- matrix(rnorm(n_train * q), n_train, q)
  colnames(C_train) <- c("C1", "C2", "C3")
  
  # True linear relationship
  beta <- c(1.5, -0.8, 2.2)
  Y_train <- as.numeric(C_train %*% beta + rnorm(n_train, sd = 0.5))
  
  C_test <- matrix(rnorm(n_test * q), n_test, q)
  colnames(C_test) <- c("C1", "C2", "C3")
  
  # 2. Define a fast wrapper for SuperLearner to keep tests speedy
  fast_SL_fitter <- function(target, C_df, ...) {
    SuperLearner::SuperLearner(Y = target, X = C_df, 
                               family = gaussian(), 
                               SL.library = "SL.glm", ...) # Use only GLM for speed
  }
  
  # 3. Fit the Model
  # We suppress warnings just in case SL complains about identical CV folds in small data
  suppressWarnings({
    fit <- nuisance_C_model(target = Y_train, 
                            C = C_train, 
                            fitter = fast_SL_fitter)
  })
  
  # --- Assertions for the Constructor ---
  expect_s3_class(fit, "nuisance_C_model")
  expect_equal(fit$q, 3)
  expect_equal(fit$C_names, c("C1", "C2", "C3"))
  
  # 4. Generate Predictions
  preds <- predict(fit, newdata = C_test)
  
  # --- Assertions for Standard Prediction ---
  # Check that SuperLearner's list output was successfully unwrapped into a numeric vector
  expect_type(preds, "double")
  expect_length(preds, n_test)
  expect_null(dim(preds)) # Should be a flat vector, not a matrix
  
  # 5. Test the Column Safety Contract
  
  # Test A: Missing Columns
  C_test_missing <- C_test[, c("C1", "C2")] # Drop C3
  expect_error(
    predict(fit, newdata = C_test_missing),
    "The following required columns are missing from 'newdata': C3"
  )
  
  # Test B: Scrambled Columns
  # The predict method should forcefully reorder these to match training
  C_test_scrambled <- C_test[, c("C3", "C1", "C2")]
  preds_scrambled <- predict(fit, newdata = C_test_scrambled)
  
  # The predictions should be identical to the correctly ordered data
  expect_equal(preds, preds_scrambled)
})
