library(MAVE)
library(MASS)

test_mave_vs_efficient <- function() {
  cat("\n--- Comparing MAVE vs. Efficient Estimator ---\n")
  #set.seed(999)
  n <- 1000
  p <- 6
  d <- 2
  
  # 1. Generate Data (Heteroscedastic Multi-Index Model)
  X <- matrix(runif(n * p, -2, 2), n, p)
  
  beta_true <- matrix(0, nrow = p, ncol = d)
  beta_true[1, 1] <- 1
  beta_true[2, 2] <- 1
  
  m_X <- X[, 1] / (0.5 + (1.5 + X[, 2])^2)
  sigma_X <- sqrt(0.5 * m_X^2 + 0.1)
  Y <- m_X + sigma_X * rnorm(n)
  
  # ---------------------------------------------------------
  # 2. Run the CRAN MAVE Package
  # ---------------------------------------------------------
  cat("Running MAVE::mave (Refined MAVE)...\n")
  # Use meanMAVE to target the Central Mean Subspace
  mave_fit <- MAVE::mave(Y ~ X, method = "meanMAVE")
  
  # Extract the d-dimensional estimate
  beta_mave <- mave_fit$dir[[d]] 
  
  dist_mave <- norm(Pi(beta_mave) - Pi(beta_true) )
  cat(sprintf("Distance of MAVE to Truth:               %f\n", dist_mave))
  
  # ---------------------------------------------------------
  # 3. Run Your Vectorized Efficient Estimator
  # ---------------------------------------------------------
  cat("\nRunning Vectorized Efficient Estimator...\n")
  
  # Calculate dynamic bandwidths based on the MAVE projection
  b_val <- 1.06 * sd(X) * n^(-1 / (d + 4))
  h_val <- 1.06 * sd(X) * n^(-1 / (4 * p))
  
  # Polish the MAVE estimate using the Efficient Score
  beta_eff <- run_efficient_estimator(X, Y, beta_init = beta_mave, 
                                      b = b_val, h = h_val, max_iters = 100, SL = FALSE)
  
  dist_eff <- norm(Pi(beta_eff) -  Pi(beta_true))
  cat(sprintf("Distance of EFFICIENT ESTIMATOR to Truth:  %f\n", dist_eff))
  
  print(beta_mave)
  print(beta_eff)
  
  # ---------------------------------------------------------
  # 4. Conclusion
  # ---------------------------------------------------------
  cat("\n--- Results ---\n")
  if (dist_eff < dist_mave) {
    improvement <- (dist_mave - dist_eff) / dist_mave * 100
    cat(sprintf("SUCCESS: Efficient Estimator improved MAVE by %.1f%%!\n", improvement))
  } else {
    cat("FAIL: Efficient Estimator did not improve upon MAVE.\n")
  }
}

# Run the comparison
test_mave_vs_efficient()
