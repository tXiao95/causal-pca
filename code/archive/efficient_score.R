library(MASS) # For ginv (Moore-Penrose generalized inverse)

# Kernel function K_{h, p}(u)
# Uses a Gaussian kernel for simplicity
kernel_weight <- function(u, h, p) {
  norm_u <- sqrt(sum(u^2))
  (1 / h^p) * dnorm(norm_u / h)
}

# Computes the projection matrix Pi(beta) = beta(beta^T beta)^{-1} beta^T
projection_matrix <- function(beta) {
  beta %*% solve(t(beta) %*% beta) %*% t(beta)
}

frob_norm <- function(A, B){
  sqrt(sum((A - B)^2))
}

# Estimates m(beta^T x) and its gradient m'(beta^T x) for all X_i
estimate_m_and_gradient <- function(X, Y, beta, b) {
  n <- nrow(X)
  d <- ncol(beta)
  betaX <- X %*% beta 
  
  m_est <- numeric(n)
  m_prime_est <- matrix(0, nrow = n, ncol = d)
  
  for (i in 1:n) {
    x_i <- X[i, ]
    beta_x_i <- betaX[i, ]
    
    # Calculate weights for all j
    weights <- apply(betaX, 1, function(bx_j) kernel_weight(bx_j - beta_x_i, b, d))
    
    # Construct design matrix for local quadratic regression
    # intercept + linear terms + quadratic terms (flattened)
    Z <- matrix(1, nrow = n, ncol = 1)
    diff_betaX <- sweep(betaX, 2, beta_x_i, "-")
    Z <- cbind(Z, diff_betaX)
    
    for (k in 1:d) {
      for (l in k:d) {
        Z <- cbind(Z, diff_betaX[, k] * diff_betaX[, l])
      }
    }
    
    # Weighted least squares
    W <- diag(weights)
    # Adding small ridge penalty for numerical stability in inversion
    coeffs <- MASS::ginv(t(Z) %*% W %*% Z + diag(1e-8, ncol(Z))) %*% t(Z) %*% W %*% Y
    
    m_est[i] <- coeffs[1]
    m_prime_est[i, ] <- coeffs[2:(d + 1)]
  }
  
  list(m_est = m_est, m_prime_est = m_prime_est)
}

# Estimates the conditional variance sigma^2(X) 
estimate_sigma2 <- function(X, residuals, h) {
  n <- nrow(X)
  p <- ncol(X)
  sigma2_est <- numeric(n)
  
  for (i in 1:n) {
    x_i <- X[i, ]
    weights <- apply(X, 1, function(x_j) kernel_weight(x_j - x_i, h, p))
    # OPTIMIZATION: Added 1e-10 to denominator to prevent division by zero
    sigma2_est[i] <- sum(weights * (residuals^2)) / (sum(weights) + 1e-10)
  }
  
  # Ensure variance is strictly positive to prevent division by zero downstream
  pmax(sigma2_est, 1e-6) 
}

# Evaluates Step 2 and returns the updated vec(beta)
compute_efficient_score_and_update <- function(X, Y, beta, b, h) {
  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(beta)
  betaX <- X %*% beta
  
  # 1. Get m and m'
  m_results <- estimate_m_and_gradient(X, Y, beta, b)
  residuals <- Y - m_results$m_est
  
  # 2. Get sigma^2
  sigma2 <- estimate_sigma2(X, residuals, h)
  
  # 3. Compute Conditional Expectations Eq (5)
  ratio_term <- matrix(0, nrow = n, ncol = p)
  for (i in 1:n) {
    beta_x_i <- betaX[i, ]
    weights <- apply(betaX, 1, function(bx_j) kernel_weight(bx_j - beta_x_i, b, d))
    
    # OPTIMIZATION: Safer row-wise multiplication using sweep
    combined_weights <- weights / sigma2
    num <- colSums(sweep(X, 1, combined_weights, "*"))
    
    # OPTIMIZATION: Prevent division by zero
    den <- sum(weights * (1 / sigma2)) + 1e-10
    ratio_term[i, ] <- num / den
  }
  
  # 4. Compute Efficient Score Eq (7)
  S_eff_list <- list()
  for (i in 1:n) {
    term1 <- (X[i, ] - ratio_term[i, ]) / sigma2[i] # vector of length p
    term2 <- m_results$m_prime_est[i, ]             # vector of length d
    term3 <- residuals[i]                           # scalar
    
    # S_eff is p x d matrix
    S_eff_list[[i]] <- term1 %*% t(term2) * term3 
  }
  
  # 5. Update beta using Newton-Raphson update
  # vec(beta^(k+1)) = vec(beta^(k)) + E_n^+[vec(S_eff)^otimes 2] E_n[vec(S_eff)]
  vec_S <- sapply(S_eff_list, as.vector) # pd x n matrix
  
  mean_vec_S <- rowMeans(vec_S) # E_n[vec(S)]
  
  # OPTIMIZATION: Vectorized outer product (much faster than a for loop)
  avg_outer_S <- (vec_S %*% t(vec_S)) / n
  
  # Update step (Addition is correct due to gradient/Fisher information equivalence)
  vec_beta_k <- as.vector(beta)
  vec_beta_next <- vec_beta_k + MASS::ginv(avg_outer_S) %*% mean_vec_S
  
  matrix(vec_beta_next, nrow = p, ncol = d)
}

# Step 3: Main driver function iterating until convergence
run_efficient_estimator <- function(X, Y, beta_init, b, h, max_iters = 100) {
  n <- nrow(X)
  p <- ncol(X)
  
  beta_current <- beta_init
  threshold <- p / n 
  
  for (k in 1:max_iters) {
    cat(sprintf("Starting iteration %d...\n", k))
    
    # Run Step 2
    beta_next <- compute_efficient_score_and_update(X, Y, beta_current, b, h)
    
    # Standard practice to maintain orthonormal basis to prevent matrix drift
    beta_next <- qr.Q(qr(beta_next))
    
    # Calculate projection distance using Frobenius norm
    pi_current <- projection_matrix(beta_current)
    pi_next <- projection_matrix(beta_next)
    
    # Frobenius norm is the sqrt of the sum of squared elements
    dist <- sqrt(sum((pi_current - pi_next)^2))
    #dist <- frob_norm
    
    cat(sprintf("Distance after iteration %d: %f (Threshold: %f)\n", k, dist, threshold))
    
    if (dist < threshold) {
      cat("Convergence threshold reached.\n")
      return(beta_next)
    }
    
    beta_current <- beta_next
  }
  
  warning("Maximum iterations reached without falling below threshold.")
  return(beta_current)
}