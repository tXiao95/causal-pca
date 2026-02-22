library(testthat)

# We will use SL.glm for fast testing. In production, you would use a broader SL.library.
fast_SL_fitter <- function(W, C) {
  SL_fitter(W, C, SL.library = c("SL.glm", "SL.mean", "SL.gam"))
}

test_that("MultiGPS returns correct object structure and predictions", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n), ncol = 1)
  C <- matrix(rnorm(n), ncol = 1)
  colnames(C) <- "C1"
  
  mod <- multigps(x_eval = c(0), X = X, C = C, mu_fitter = fast_SL_fitter)
  
  # Check class
  expect_s3_class(mod, "multigps")
  expect_equal(mod$p, 1)
  
  # Check prediction shape
  C_new <- data.frame(C1 = c(0.5, -0.5))
  preds <- predict(mod, newdata = C_new)
  
  expect_type(preds, "double")
  expect_length(preds, 2)
})

test_that("MultiGPS estimates obvious independent marginal density", {
  # Scenario: X is completely independent of C. X ~ N(0, 1).
  # Therefore, f(x | c) = f(x) = dnorm(x).
  set.seed(123)
  n <- 2000
  C <- matrix(runif(n, -1, 1), ncol = 1)
  colnames(C) <- "C1"
  X <- matrix(rnorm(n, mean = 0, sd = 1), ncol = 1)
  
  # Target 1: x_eval = 0. True density = dnorm(0) ≈ 0.3989
  mod_0 <- multigps(x_eval = c(0), X = X, C = C, mu_fitter = fast_SL_fitter)
  pred_0 <- predict(mod_0, newdata = data.frame(C1 = 0))
  
  # Target 2: x_eval = 1. True density = dnorm(1) ≈ 0.2419
  mod_1 <- multigps(x_eval = c(1), X = X, C = C, mu_fitter = fast_SL_fitter)
  pred_1 <- predict(mod_1, newdata = data.frame(C1 = 0))
  
  # Allow a tolerance of 0.05 for kernel smoothing bias
  expect_equal(pred_0[1], dnorm(0), tolerance = 0.05)
  expect_equal(pred_1[1], dnorm(1), tolerance = 0.05)
})

test_that("MultiGPS estimates obvious conditional density (Linear Shift)", {
  # Scenario: X depends on C exactly. X ~ N(C, 1).
  # Therefore, f(x | C=c) = dnorm(x, mean = c, sd = 1).
  set.seed(456)
  n <- 2000
  C <- matrix(rnorm(n, mean = 0, sd = 2), ncol = 1)
  colnames(C) <- "C1"
  X <- matrix(C + rnorm(n, mean = 0, sd = 1), ncol = 1)
  
  # Target: Evaluate the density of X=2 given C=2.
  # Since X|C=2 ~ N(2, 1), the density at X=2 should be the peak: dnorm(2, 2, 1) = dnorm(0) ≈ 0.3989
  mod <- multigps(x_eval = c(2), X = X, C = C, mu_fitter = fast_SL_fitter)
  pred <- predict(mod, newdata = data.frame(C1 = 2))
  
  expect_equal(pred[1], dnorm(0), tolerance = 0.05)
  
  # Conversely, the density of X=0 given C=2 is dnorm(0, 2, 1) = dnorm(2) ≈ 0.054
  mod_tail <- multigps(x_eval = c(0), X = X, C = C, mu_fitter = fast_SL_fitter)
  pred_tail <- predict(mod_tail, newdata = data.frame(C1 = 2))
  
  expect_equal(pred_tail[1], dnorm(2), tolerance = 0.05)
})

test_that("MultiGPS works for bivariate treatments (p = 2)", {
  # Scenario: X1 and X2 are independent standard normals.
  # True joint density at (0,0) is dnorm(0) * dnorm(0) ≈ 0.159
  set.seed(789)
  n <- 3000
  C <- matrix(runif(n), ncol = 1)
  colnames(C) <- "C1"
  X <- cbind(rnorm(n), rnorm(n))
  
  mod <- multigps(x_eval = c(0, 0), X = X, C = C, mu_fitter = fast_SL_fitter)
  pred <- predict(mod, newdata = data.frame(C1 = 0.5))
  
  true_bivariate_density <- dnorm(0) * dnorm(0)
  expect_equal(pred[1], true_bivariate_density, tolerance = 0.05)
})

test_that("MultiGPS predict correctly applies flooring (delta_n)", {
  set.seed(999)
  n <- 500
  C <- matrix(runif(n), ncol = 1)
  colnames(C) <- "C1"
  X <- matrix(rnorm(n), ncol = 1)
  
  # Target an extreme outlier: x_eval = 20. 
  # There are no observations near 20, so the pseudo-outcome W will be effectively 0 everywhere.
  # The raw predicted density will be extremely tiny (or slightly negative due to linear model extrapolation).
  mod <- multigps(x_eval = c(20), X = X, C = C, mu_fitter = fast_SL_fitter)
  
  custom_delta <- 1e-3
  pred <- predict(mod, newdata = data.frame(C1 = 0.5), delta_n = custom_delta)
  
  # It should be floored exactly to the delta_n value
  expect_equal(pred[1], custom_delta)
})
