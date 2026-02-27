library(testthat)

test_that("mvnGPS returns correct object structure and predictions", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n), ncol = 1)
  C <- matrix(rnorm(n), ncol = 1)
  colnames(C) <- "C1"
  
  mod <- mvnGPS(x_eval = c(0), X = X, C = C)
  
  # Check class and structure
  expect_s3_class(mod, "mvnGPS")
  expect_equal(mod$p, 1)
  expect_true(!is.null(mod$sigma_hat))
  
  # Check prediction shape
  C_new <- data.frame(C1 = c(0.5, -0.5))
  preds <- predict(mod, newdata = C_new)
  
  expect_type(preds, "double")
  expect_length(preds, 2)
})

test_that("mvnGPS estimates obvious independent marginal density (Univariate)", {
  # Scenario: X is completely independent of C. X ~ N(0, 1).
  set.seed(123)
  n <- 5000 # Large N to ensure sample variance is very close to 1
  C <- matrix(runif(n, -1, 1), ncol = 1)
  X <- matrix(rnorm(n, mean = 0, sd = 1), ncol = 1)
  colnames(C) <- "C1"
  
  # Target 1: x_eval = 0. True density = dnorm(0) ≈ 0.3989
  mod_0 <- mvnGPS(x_eval = c(0), X = X, C = C)
  pred_0 <- unname( predict(mod_0, newdata = data.frame(C1 = 0)) )
  
  # Target 2: x_eval = 1. True density = dnorm(1) ≈ 0.2419
  mod_1 <- mvnGPS(x_eval = c(1), X = X, C = C)
  pred_1 <- unname( predict(mod_1, newdata = data.frame(C1 = 0)) )
  
  # We can use a much tighter tolerance here (0.01) because there is no kernel smoothing bias
  expect_equal(pred_0[1], dnorm(0), tolerance = 0.01)
  expect_equal(pred_1[1], dnorm(1), tolerance = 0.01)
})

test_that("mvnGPS estimates obvious conditional density (Linear Shift)", {
  # Scenario: X depends on C exactly. X ~ N(C, 1).
  set.seed(456)
  n <- 5000
  C <- matrix(rnorm(n, mean = 0, sd = 2), ncol = 1)
  colnames(C) <- "C1"
  X <- matrix(C + rnorm(n, mean = 0, sd = 1), ncol = 1)
  
  # Target: Evaluate the density of X=2 given C=2.
  # Since X|C=2 ~ N(2, 1), the density at X=2 should be the peak: dnorm(2, 2, 1) = dnorm(0)
  mod <- mvnGPS(x_eval = c(2), X = X, C = C)
  pred <- unname(predict(mod, newdata = data.frame(C1 = 2)))
  
  expect_equal(pred[1], dnorm(0), tolerance = 0.01)
  
  # Conversely, the density of X=0 given C=2 is dnorm(0, 2, 1) = dnorm(2)
  mod_tail <- mvnGPS(x_eval = c(0), X = X, C = C)
  pred_tail <- unname(predict(mod_tail, newdata = data.frame(C1 = 2)))
  
  expect_equal(pred_tail[1], dnorm(2), tolerance = 0.01)
})

test_that("mvnGPS works for bivariate independent treatments (p = 2)", {
  # Scenario: X1 and X2 are independent standard normals.
  set.seed(789)
  n <- 5000
  C <- matrix(runif(n), ncol = 1)
  colnames(C) <- "C1"
  X <- cbind(rnorm(n), rnorm(n))
  
  mod <- mvnGPS(x_eval = c(0, 0), X = X, C = C)
  pred <- unname(predict(mod, newdata = data.frame(C1 = 0.5)))
  
  true_bivariate_density <- dnorm(0) * dnorm(0) # ~ 0.1591
  expect_equal(pred[1], true_bivariate_density, tolerance = 0.01)
})

test_that("mvnGPS works for correlated bivariate treatments", {
  # Scenario: X1 and X2 are correlated. X2 = 0.5*X1 + noise
  set.seed(101)
  n <- 5000
  C <- matrix(runif(n), ncol = 1)
  colnames(C) <- "C1"
  
  X1 <- rnorm(n)
  X2 <- 0.5 * X1 + rnorm(n, sd = sqrt(0.75)) # Var(X2) = 1, Cov(X1,X2) = 0.5
  X <- cbind(X1, X2)
  
  # True covariance matrix
  true_sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  
  mod <- mvnGPS(x_eval = c(0, 0), X = X, C = C)
  pred <- predict(mod, newdata = data.frame(C1 = 0))
  
  true_correlated_density <- mvtnorm::dmvnorm(c(0, 0), mean = c(0, 0), sigma = true_sigma)
  expect_equal(pred[1], true_correlated_density, tolerance = 0.01)
})

test_that("predict.mvnGPS correctly applies flooring (delta_n)", {
  set.seed(999)
  n <- 500
  C <- matrix(runif(n), ncol = 1)
  colnames(C) <- "C1"
  X <- matrix(rnorm(n), ncol = 1)
  
  # Target an extreme outlier: x_eval = 20. 
  # The true normal PDF for an event 20 standard deviations away is virtually 0 (approx 5.5e-89).
  mod <- mvnGPS(x_eval = c(20), X = X, C = C)
  
  custom_delta <- 1e-4
  pred <- predict(mod, newdata = data.frame(C1 = 0.5), delta_n = custom_delta)
  
  # Because the raw density is astronomically small, pmax() should replace it with custom_delta
  expect_equal(pred[1], custom_delta)
})
