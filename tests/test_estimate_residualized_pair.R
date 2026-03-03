library(testthat)
library(here)

# (Insert your nuisance_C_model, predict.nuisance_C_model, and estimate_residualized_pair here)
source("R/nuisance_outcome_C.R")
source("R/estimate_residualized_pair.R")

test_that("Pipeline integration: train_nuisance_models feeds estimate_residualized_pair", {
  
  # 1. Simulate data
  set.seed(42)
  n <- 200
  p <- 3
  q <- 2
  
  C <- matrix(rnorm(n * q), n, q)
  colnames(C) <- c("C1", "C2")
  
  # Generate X with a linear dependence on C
  beta_X <- matrix(c(0.5, -0.2, 
                     -0.5,  0.4, 
                     0.1,  0.8), nrow = q, ncol = p, byrow = TRUE)
  X <- C %*% beta_X + matrix(rnorm(n * p), n, p)
  colnames(X) <- c("Exposure_1", "Exposure_2", "Exposure_3")
  
  # Generate Y with a linear dependence on C and X
  Y <- as.numeric(C %*% c(1, -1) + X %*% c(0.5, 1.5, -0.5) + rnorm(n))
  
  # 2. Define a fast fitter compatible with nuisance_C_model
  fast_lm_fitter <- function(target, C_df, ...) {
    # Combine into one dataframe for lm()
    df <- cbind(target = target, C_df)
    lm(target ~ ., data = df)
  }
  
  # ---------------------------------------------------------
  # Step A: Train the models using your wrapper
  # ---------------------------------------------------------
  trained_models <- train_nuisance_models(Y = Y, 
                                          X = X, 
                                          C = C, 
                                          fitter = fast_lm_fitter)
  
  # Assert the training function built the correct S3 structure
  expect_type(trained_models, "list")
  expect_named(trained_models, c("Y_model", "X_models"))
  expect_s3_class(trained_models$Y_model, "nuisance_C_model")
  expect_length(trained_models$X_models, p)
  expect_s3_class(trained_models$X_models[[1]], "nuisance_C_model")
  
  # ---------------------------------------------------------
  # Step B: Pass the trained models to the residualization function
  # ---------------------------------------------------------
  res <- estimate_residualized_pair(Y = Y, 
                                    X = X, 
                                    C = C, 
                                    C_models = trained_models)
  
  # ---------------------------------------------------------
  # Assertions on the final output
  # ---------------------------------------------------------
  
  # Check Structure
  expect_type(res, "list")
  expect_named(res, c("Ytilde", "Xtilde"))
  expect_length(res$Ytilde, n)
  expect_equal(dim(res$Xtilde), c(n, p))
  
  # Check Column Naming Contract
  expect_true(!is.null(colnames(res$Xtilde)))
  expect_equal(colnames(res$Xtilde), c("Exposure_1", "Exposure_2", "Exposure_3"))
  
  # Check Mean-Zero Residuals (In-sample OLS property)
  expect_lt(abs(mean(res$Ytilde)), 1e-10)
  
  col_means_Xtilde <- colMeans(res$Xtilde)
  for (i in 1:p) {
    expect_lt(abs(col_means_Xtilde[i]), 1e-10)
  }
})