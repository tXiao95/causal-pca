library(MASS) # For ginv

# Vectorized Gaussian kernel matching the paper's K_{h, p}(u)
# Retains the 1/h^p scaling to keep the ridge penalty balanced
kernel_weight_vec <- function(diff_mat, h, p) {
  norm_u <- sqrt(rowSums(diff_mat^2))
  (1 / h^p) * dnorm(norm_u / h)
}

# Computes the projection matrix Pi(beta)
projection_matrix <- function(beta) {
  beta %*% solve(t(beta) %*% beta) %*% t(beta)
}

frob_norm <- function(A,B){
  sqrt(sum((A-B)^2))
}

# 1. Fast Vectorized Variance Estimator (Eliminates the 500ms estimate_sigma2 bottleneck)
estimate_sigma2 <- function(X, residuals, h) {
  p <- ncol(X)
  
  # Calculate all pairwise distances at once in C
  dist_mat <- as.matrix(dist(X, method = "euclidean"))
  
  # Apply kernel to the entire matrix
  weight_mat <- (1 / h^p) * dnorm(dist_mat / h)
  
  num <- weight_mat %*% (residuals^2)
  den <- rowSums(weight_mat) + 1e-10
  
  as.vector(num / den)
}

# 2. Fully Vectorized Main Update (Eliminates the 2400ms ratio_term loop bottlenecks)
compute_efficient_score_and_update <- function(X, Y, beta, b, h) {
  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(beta)
  betaX <- X %*% beta
  
  # Fast C++ m and gradient calculation
  m_results <- estimate_m_gradient_cpp(X, Y, beta, b)
  residuals <- Y - m_results$m_est
  
  # Fast Variance Calculation
  sigma2 <- estimate_sigma2(X, residuals, h)
  
  # --- THE MAJOR SPEEDUP: Vectorized Ratio Term ---
  # Calculate pairwise distances of projected X for all points at once
  dist_mat_betaX <- as.matrix(dist(betaX, method = "euclidean"))
  
  # Calculate the full n x n weight matrix simultaneously
  W <- (1 / b^d) * dnorm(dist_mat_betaX / b)
  
  # Scale columns by sigma2 using an optimized C-level sweep
  W_combined <- sweep(W, 2, sigma2, "/")
  
  # The entire ratio_term loop condenses into one BLAS-optimized matrix multiplication!
  ratio_num <- W_combined %*% X
  ratio_den <- rowSums(W_combined) + 1e-10
  ratio_term <- ratio_num / ratio_den
  # ------------------------------------------------
  
  # Assembly of the efficient score vectors
  term1_mat <- (X - ratio_term) / sigma2 
  term3_vec <- residuals
  
  # lapply is highly optimized for list creation in R
  S_eff_list <- lapply(1:n, function(i) {
    term1 <- term1_mat[i, ]
    term2 <- m_results$m_prime_est[i, ]
    as.vector(term1 %*% t(term2) * term3_vec[i])
  })
  
  # Fast matrix assembly and Newton-Raphson math
  vec_S <- do.call(cbind, S_eff_list) 
  mean_vec_S <- rowMeans(vec_S) 
  avg_outer_S <- (vec_S %*% t(vec_S)) / n
  
  vec_beta_k <- as.vector(beta)
  
  # Slight ridge to prevent singular inversions on sparse steps
  vec_beta_next <- vec_beta_k + MASS::ginv(avg_outer_S + diag(1e-8, p*d)) %*% mean_vec_S
  
  matrix(vec_beta_next, nrow = p, ncol = d)
}

# Main iterative loop (No QR decomposition inside the loop)
run_efficient_estimator <- function(X, Y, beta_init, b, h, max_iters = 100) {
  n <- nrow(X)
  p <- ncol(X)
  
  beta_current <- beta_init
  threshold <- p / n 
  
  for (k in 1:max_iters) {
    cat(sprintf("Starting iteration %d...\n", k))
    
    beta_next <- compute_efficient_score_and_update(X, Y, beta_current, b, h)
    
    pi_current <- projection_matrix(beta_current)
    pi_next <- projection_matrix(beta_next)
    
    dist <- sqrt(sum((pi_current - pi_next)^2))
    cat(sprintf("Distance after iteration %d: %f (Threshold: %f)\n", k, dist, threshold))
    
    if (dist < threshold) {
      cat("Convergence threshold reached.\n")
      # Return the final estimate
      beta_next <- qr.Q(qr(beta_next))
      return(beta_next)
    }
    beta_current <- beta_next
  }
  
  warning("Maximum iterations reached without falling below threshold.")
  beta_current <- qr.Q(qr(beta_current))
  return(beta_current)
}