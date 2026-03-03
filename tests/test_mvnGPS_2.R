library(testthat)
library(here)
library(mvtnorm)

source(here("R","mvnGPS.R"))

test_that("mvnGPS accurately predicts true MVN density", {
  
  # 1. Set up the true parameters and simulate data
  set.seed(42)
  n_train <- 100000
  n_test <- 1000
  p <- 12  # Number of treatments
  q <- 20  # Number of confounders
  
  # True coefficients (q x p)
  # Dynamically generated to match dimensions p and q
  beta_true <- matrix(runif(q * p, min = -1, max = 1), nrow = q, ncol = p)
  
  # True Covariance matrix (Sigma)
  # Create a positive definite matrix
  A <- matrix(runif(p^2, -1, 1), p, p)
  Sigma_true <- t(A) %*% A + diag(p) 
  
  # Generate Training Data
  C_train <- matrix(rnorm(n_train * q), n_train, q)
  mu_train <- C_train %*% beta_true
  X_train <- mu_train + rmvnorm(n_train, sigma = Sigma_true)
  
  # Generate Testing Data
  C_test <- matrix(rnorm(n_test * q), n_test, q)
  mu_test <- C_test %*% beta_true
  X_test <- mu_test + rmvnorm(n_test, sigma = Sigma_true)
  
  # 2. Fit the mvnGPS model
  fit <- mvnGPS(X = X_train, C = C_train, method = "linear")
  
  # 3. Predict densities using the custom predict method
  newdata <- as.data.frame(cbind(X_test, C_test))
  predicted_densities <- predict(fit, newdata = newdata, delta_n = 1e-16)
  
  # 4. Calculate the true theoretical densities
  true_densities <- numeric(n_test)
  for (i in 1:n_test) {
    true_densities[i] <- dmvnorm(x = X_test[i, ], 
                                 mean = mu_test[i, ], 
                                 sigma = Sigma_true)
  }
  
  # 5. Evaluate Accuracy and Consistency
  
  # Assert that the correlation between predicted and true densities is highly positive
  pearson_cor <- cor(predicted_densities, true_densities)
  expect_gt(pearson_cor, 0.99) 
  
  # Assert that the Mean Absolute Error (MAE) is sufficiently small
  # (Note: With high dimensions like p=12, multivariate densities become extremely small, 
  # so this MAE check is virtually guaranteed to pass. You might want to evaluate log-densities instead).
  mae <- mean(abs(predicted_densities - true_densities))
  expect_lt(mae, 0.005)
  
  # Assert that the estimated Sigma is close to the true Sigma
  expect_equal(as.numeric(fit$sigma_hat), as.numeric(Sigma_true), tolerance = 0.05)
})
