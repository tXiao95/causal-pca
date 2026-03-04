# 1. Robust Gaussian Weights (Unnormalized to prevent underflow)
get_gaussian_weights <- function(diff_mat, bandwidth) {
  # diff_mat should already be X_j - x_i
  u <- sqrt(rowSums(diff_mat^2)) / bandwidth
  # We drop the 1/h^p term because W is proportional and it cancels out in WLS
  weights <- exp(-0.5 * u^2) 
  
  # Prevent absolute zero weights for numerical stability
  pmax(weights, 1e-10)
}

# 2. Stabilized Local Quadratic Regression (Gaussian)
estimate_m_and_gradient <- function(X, Y, beta, b) {
  n <- nrow(X)
  d <- ncol(beta)
  betaX <- X %*% beta 
  
  m_est <- numeric(n)
  m_prime_est <- matrix(0, nrow = n, ncol = d)
  
  for (i in 1:n) {
    beta_x_i <- betaX[i, ]
    
    diff_betaX <- sweep(betaX, 2, beta_x_i, "-")
    weights <- get_gaussian_weights(diff_betaX, b)
    
    # Scale distances by bandwidth for numerical stability
    u_mat <- diff_betaX / b 
    
    Z <- matrix(1, nrow = n, ncol = 1)
    Z <- cbind(Z, u_mat)
    
    for (k in 1:d) {
      for (l in k:d) {
        Z <- cbind(Z, u_mat[, k] * u_mat[, l])
      }
    }
    
    W <- diag(weights)
    
    # Heavier ridge penalty to ensure invertibility
    coeffs <- MASS::ginv(t(Z) %*% W %*% Z + diag(1e-5, ncol(Z))) %*% t(Z) %*% W %*% Y
    
    m_est[i] <- coeffs[1]
    # Divide linear coefficients by b to recover true gradient scale
    m_prime_est[i, ] <- coeffs[2:(d + 1)] / b 
  }
  
  list(m_est = m_est, m_prime_est = m_prime_est)
}

# 3. Stabilized Variance Estimation (Gaussian)
estimate_sigma2 <- function(X, residuals, h) {
  n <- nrow(X)
  sigma2_est <- numeric(n)
  
  for (i in 1:n) {
    diff_X <- sweep(X, 2, X[i, ], "-")
    weights <- get_gaussian_weights(diff_X, h)
    
    sigma2_est[i] <- sum(weights * (residuals^2)) / sum(weights)
  }
  
  # Keep variance strictly positive
  pmax(sigma2_est, 1e-6) 
}

# 4. Routine Update (Ensure compute_efficient_score_and_update calls the new weights)
compute_efficient_score_and_update <- function(X, Y, beta, b, h) {
  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(beta)
  betaX <- X %*% beta
  
  m_results <- estimate_m_and_gradient(X, Y, beta, b)
  residuals <- Y - m_results$m_est
  sigma2 <- estimate_sigma2(X, residuals, h)
  
  ratio_term <- matrix(0, nrow = n, ncol = p)
  for (i in 1:n) {
    beta_x_i <- betaX[i, ]
    diff_betaX <- sweep(betaX, 2, beta_x_i, "-")
    weights <- get_gaussian_weights(diff_betaX, b)
    
    combined_weights <- weights / sigma2
    num <- colSums(sweep(X, 1, combined_weights, "*"))
    den <- sum(weights * (1 / sigma2))
    ratio_term[i, ] <- num / den
  }
  
  S_eff_list <- list()
  for (i in 1:n) {
    term1 <- (X[i, ] - ratio_term[i, ]) / sigma2[i] 
    term2 <- m_results$m_prime_est[i, ]             
    term3 <- residuals[i]                           
    S_eff_list[[i]] <- term1 %*% t(term2) * term3 
  }
  
  vec_S <- sapply(S_eff_list, as.vector) 
  mean_vec_S <- rowMeans(vec_S) 
  avg_outer_S <- (vec_S %*% t(vec_S)) / n
  
  vec_beta_k <- as.vector(beta)
  # Ridge added to Fisher Information matrix to prevent singular steps
  vec_beta_next <- vec_beta_k + MASS::ginv(avg_outer_S + diag(1e-6, p*d)) %*% mean_vec_S
  
  matrix(vec_beta_next, nrow = p, ncol = d)
}