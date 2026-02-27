# tests/testthat/test-compute_score_eq7.R

test_that("compute_score_eq7 returns correct shapes and matches a reference implementation (n > p)", {
  
  set.seed(123)
  
  n <- 50
  p <- 8
  d <- 2
  
  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)
  
  m_hat <- rnorm(n)
  mprime_hat <- matrix(rnorm(n * d), n, d)
  
  sigma2_hat <- abs(rnorm(n)) + 0.5            # positive
  EX_over_sigma2 <- matrix(rnorm(n * p), n, p)
  E_inv_sigma2 <- abs(rnorm(n)) + 0.1          # positive
  
  out <- compute_score_eq7(
    X, Y, m_hat, mprime_hat,
    sigma2_hat, EX_over_sigma2, E_inv_sigma2
  )
  
  # ---- dimension checks ----
  expect_equal(dim(out$S_i), c(n, p, d))
  expect_equal(dim(out$S_bar), c(p, d))
  expect_equal(length(out$vec_S_bar), p * d)
  
  # ---- independent reference implementation ----
  ref_S_i <- array(0, dim = c(n, p, d))
  for (i in 1:n) {
    ratio_i <- EX_over_sigma2[i, ] / E_inv_sigma2[i]
    bracket_i <- X[i, ] - ratio_i
    g_i <- bracket_i / sigma2_hat[i]
    r_i <- Y[i] - m_hat[i]
    ref_S_i[i, , ] <- (g_i %o% mprime_hat[i, ]) * r_i
  }
  ref_S_bar <- apply(ref_S_i, c(2, 3), mean)
  ref_vec <- as.vector(ref_S_bar)
  
  expect_equal(out$S_i, ref_S_i, tolerance = 1e-12)
  expect_equal(out$S_bar, ref_S_bar, tolerance = 1e-12)
  expect_equal(out$vec_S_bar, ref_vec, tolerance = 1e-12)
})

test_that("compute_score_eq7 returns (near) zero when residuals are identically zero", {
  
  set.seed(456)
  
  n <- 40
  p <- 6
  d <- 3
  
  X <- matrix(rnorm(n * p), n, p)
  m_hat <- rnorm(n)
  Y <- m_hat                                 # forces residuals = 0
  
  mprime_hat <- matrix(rnorm(n * d), n, d)
  sigma2_hat <- abs(rnorm(n)) + 0.5
  EX_over_sigma2 <- matrix(rnorm(n * p), n, p)
  E_inv_sigma2 <- abs(rnorm(n)) + 0.1
  
  out <- compute_score_eq7(
    X, Y, m_hat, mprime_hat,
    sigma2_hat, EX_over_sigma2, E_inv_sigma2
  )
  
  expect_true(max(abs(out$vec_S_bar)) < 1e-12)
  expect_true(max(abs(out$S_bar)) < 1e-12)
})

test_that("compute_score_eq7 is invariant to common scaling of EX_over_sigma2 and E_inv_sigma2", {
  
  set.seed(789)
  
  n <- 60
  p <- 7
  d <- 2
  
  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)
  
  m_hat <- rnorm(n)
  mprime_hat <- matrix(rnorm(n * d), n, d)
  
  sigma2_hat <- abs(rnorm(n)) + 0.5
  EX_over_sigma2 <- matrix(rnorm(n * p), n, p)
  E_inv_sigma2 <- abs(rnorm(n)) + 0.1
  
  out1 <- compute_score_eq7(
    X, Y, m_hat, mprime_hat,
    sigma2_hat, EX_over_sigma2, E_inv_sigma2
  )
  
  cscale <- 10
  out2 <- compute_score_eq7(
    X, Y, m_hat, mprime_hat,
    sigma2_hat, EX_over_sigma2 * cscale, E_inv_sigma2 * cscale
  )
  
  expect_equal(out1$vec_S_bar, out2$vec_S_bar, tolerance = 1e-12)
  expect_equal(out1$S_bar, out2$S_bar, tolerance = 1e-12)
})
