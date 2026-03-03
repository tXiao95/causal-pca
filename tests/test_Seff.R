library(testthat)
library(MASS)

source("R/Seff.R")

# -------------------------------------------------------------------
# Assuming estimate_m_gradient_cpp is already compiled and in your environment
# -------------------------------------------------------------------

test_that("estimate_m_gradient perfectly recovers a linear function", {
  n <- 150
  p <- 5
  d <- 2
  set.seed(111)
  
  X <- matrix(rnorm(n * p), n, p)
  # Create a random orthogonal projection matrix
  beta <- qr.Q(qr(matrix(rnorm(p * d), p, d))) 
  U <- X %*% beta
  
  # True function: Y = 3 * U1 - 1.5 * U2
  gamma <- c(3, -1.5)
  Y <- U %*% gamma
  
  # Bandwidth
  b <- 1.0 
  
  res <- estimate_m_gradient_cpp(X, Y, beta, b)
  
  # 1. The estimated mean should exactly equal the true Y
  expect_equal(as.vector(res$m_est), as.vector(Y), tolerance = 1e-5)
  
  # 2. The estimated gradient should exactly equal gamma for every point
  expected_grad <- matrix(rep(gamma, each = n), n, d)
  expect_equal(as.matrix(res$m_prime_est), expected_grad, tolerance = 1e-5)
})


test_that("estimate_m_gradient perfectly recovers a quadratic function", {
  n <- 200
  p <- 3
  d <- 1
  set.seed(222)
  
  X <- matrix(rnorm(n * p), n, p)
  beta <- matrix(c(1, 0, 0), p, d)
  U <- X %*% beta
  
  # True function: Y = 2 * U^2 - 0.5 * U + 1
  Y <- 2 * U^2 - 0.5 * U + 1
  
  # Use a large bandwidth so neighborhoods are full and well-conditioned
  b <- 2.5 
  
  res <- estimate_m_gradient_cpp(X, Y, beta, b)
  
  # 1. The estimated mean should exactly equal the true quadratic Y
  expect_equal(as.vector(res$m_est), as.vector(Y), tolerance = 1e-4)
  
  # 2. The true analytical derivative is 4 * U - 0.5
  expected_grad <- 4 * U - 0.5
  expect_equal(as.matrix(res$m_prime_est), expected_grad, tolerance = 1e-4)
})


test_that("estimate_m_gradient matches base R matrix algebra exactly", {
  n <- 100
  p <- 4
  d <- 2
  set.seed(333)
  
  X <- matrix(rnorm(n * p), n, p)
  beta <- qr.Q(qr(matrix(rnorm(p * d), p, d)))
  U <- X %*% beta
  
  # Non-linear function with noise
  Y <- sin(U[, 1]) + cos(U[, 2]) + rnorm(n, 0, 0.5)
  b <- 0.75
  
  # Run the C++ function
  res <- estimate_m_gradient_cpp(X, Y, beta, b)
  
  # ---------------------------------------------------------
  # Manually compute the WLS solution for the FIRST observation (i = 1)
  # ---------------------------------------------------------
  u1 <- U[1, ]
  diff_U <- sweep(U, 2, u1, "-")
  
  # Kernel Weights
  dists <- sqrt(rowSums(diff_U^2))
  h_pow <- b^d
  inv_sqrt_2pi <- 1.0 / sqrt(2.0 * pi)
  weights <- (1.0 / h_pow) * inv_sqrt_2pi * exp(-0.5 * (dists / b)^2)
  
  # Design Matrix Z (Constant, Linear, and Quadratic terms)
  Z <- cbind(
    1, 
    diff_U, 
    diff_U[, 1] * diff_U[, 1], 
    diff_U[, 1] * diff_U[, 2], 
    diff_U[, 2] * diff_U[, 2]
  )
  
  # Weighted Least Squares with Ridge Penalty
  ZW <- Z * weights
  ZTWZ <- t(Z) %*% ZW
  ZTWY <- t(ZW) %*% Y
  ridge <- diag(1e-8, ncol(Z))
  
  coeffs_manual <- solve(ZTWZ + ridge, ZTWY)
  
  # Compare the C++ output to the manual Base R math
  expect_equal(res$m_est[1], coeffs_manual[1], tolerance = 1e-8)
  expect_equal(res$m_prime_est[1, 1], coeffs_manual[2], tolerance = 1e-8)
  expect_equal(res$m_prime_est[1, 2], coeffs_manual[3], tolerance = 1e-8)
})


# Sigma2 function ---------------------------------------------------------

test_that("estimate_sigma2 recovers constant variance perfectly", {
  n <- 100
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  
  # If every residual is 3, the squared residual is 9.
  # A weighted average of 9s must always be exactly 9.
  constant_res <- rep(3, n)
  h <- 1.5
  
  sigma2_est <- estimate_sigma2(X, constant_res, h)
  
  expect_equal(sigma2_est, rep(9, n), tolerance = 1e-10)
})


test_that("estimate_sigma2 approaches global mean as h -> infinity", {
  n <- 50
  p <- 2
  X <- matrix(rnorm(n * p), n, p)
  residuals <- rnorm(n, mean = 0, sd = 2)
  
  # A massive bandwidth forces all weights to be nearly identical
  h_inf <- 1e10 
  
  sigma2_est <- estimate_sigma2(X, residuals, h_inf)
  global_mean_sq_res <- mean(residuals^2)
  
  # Every point's estimated variance should be the global average
  expect_equal(sigma2_est, rep(global_mean_sq_res, n), tolerance = 1e-5)
})


test_that("estimate_sigma2 approaches exact interpolation as h -> 0", {
  n <- 50
  p <- 2
  # Use a uniform grid/spread to ensure no two points are identical
  X <- matrix(runif(n * p, -10, 10), n, p)
  residuals <- rnorm(n, mean = 0, sd = 2)
  
  # A microscopic bandwidth isolates the point itself (dist = 0)
  # All off-diagonal weights will underflow to 0
  h_zero <- 1e-6 
  
  sigma2_est <- estimate_sigma2(X, residuals, h_zero)
  
  # The estimated variance should just be the point's own squared residual
  expect_equal(sigma2_est, residuals^2, tolerance = 1e-5)
})


test_that("estimate_sigma2 exactly matches manual loop-based Nadaraya-Watson", {
  n <- 80
  p <- 4
  #set.seed(42)
  X <- matrix(rnorm(n * p), n, p)
  residuals <- rnorm(n, mean = 0, sd = 1.5)
  h <- 0.8
  
  # 1. Run the fast vectorized function
  fast_est <- estimate_sigma2(X, residuals, h)
  
  # 2. Run a mathematically explicit, slow for-loop
  slow_est <- numeric(n)
  for (i in 1:n) {
    x_i <- X[i, ]
    # Calculate Euclidean distance from point i to all points
    dists <- sqrt(rowSums(sweep(X, 2, x_i, "-")^2))
    
    # Kernel weights
    w <- (1 / h^p) * dnorm(dists / h)
    
    # NW formula
    slow_est[i] <- sum(w * residuals^2) / (sum(w) + 1e-10)
  }
  
  # 3. They must be mathematically identical
  expect_equal(fast_est, slow_est, tolerance = 1e-10)
})
