# Estimate mean and gradient
# estimate_m_and_gradient <- function(X, Y, beta, b) {
#   n <- nrow(X)
#   d <- ncol(beta)
#   betaX <- X %*% beta
# 
#   m_est <- numeric(n)
#   m_prime_est <- matrix(0, nrow = n, ncol = d)
# 
#   # 1. Pre-allocate the Z matrix memory once
#   p_z <- 1 + d + (d * (d + 1)) / 2
#   Z <- matrix(1, nrow = n, ncol = p_z)
# 
#   # Pre-allocate the ridge penalty matrix
#   ridge_pen <- diag(1e-8, p_z)
# 
#   for (i in 1:n) {
#     beta_x_i <- betaX[i, ]
# 
#     diff_betaX <- sweep(betaX, 2, beta_x_i, "-")
#     weights <- kernel_weight_vec(diff_betaX, b, d)
# 
#     # 2. Fill Z in-place (No cbind copying overhead)
#     Z[, 2:(d + 1)] <- diff_betaX
#     col_idx <- d + 2
#     for (k in 1:d) {
#       for (l in k:d) {
#         Z[, col_idx] <- diff_betaX[, k] * diff_betaX[, l]
#         col_idx <- col_idx + 1
#       }
#     }
# 
#     # 3. Vectorized Row-wise Multiplication
#     # In R, `Z * weights` naturally recycles the weights vector column-by-column.
#     # This entirely avoids creating the massive n x n diagonal matrix W.
#     ZW <- Z * weights
# 
#     # 4. Use crossprod (optimized BLAS) instead of t(A) %*% B
#     ZTWZ <- crossprod(Z, ZW)
#     ZTWY <- crossprod(ZW, Y)
# 
#     # 5. Use solve() instead of ginv(). The ridge penalty ensures invertibility.
#     coeffs <- solve(ZTWZ + ridge_pen, ZTWY)
# 
#     m_est[i] <- coeffs[1]
#     m_prime_est[i, ] <- coeffs[2:(d + 1)]
#   }
# 
#   list(m_est = m_est, m_prime_est = m_prime_est)
# }

# # Fast conditional variance estimation
# estimate_sigma2 <- function(X, residuals, h) {
#   p <- ncol(X)
#   
#   # 1. Compute the n x n matrix of all pairwise Euclidean distances at once.
#   # dist() is written in C and is extremely fast. as.matrix() makes it n x n.
#   dist_mat <- as.matrix(dist(X, method = "euclidean"))
#   
#   # 2. Apply the kernel function to the entire distance matrix simultaneously.
#   # This replaces the need for the kernel_weight_vec function inside a loop.
#   weight_mat <- (1 / h^p) * dnorm(dist_mat / h)
#   
#   # 3. Compute the numerators: sum(weights * residuals^2) for all i
#   # Matrix multiplication %*% inherently does this sum-product for every row
#   num <- weight_mat %*% (residuals^2)
#   
#   # 4. Compute the denominators: sum(weights) for all i
#   den <- rowSums(weight_mat) + 1e-10
#   
#   # 5. Divide and floor the variances
#   # as.vector flattens the resulting n x 1 matrix back to a standard numeric vector
#   sigma2_est <- as.vector(num / den)
#   
#   pmax(sigma2_est, 1e-6)
# }
# 
# # Fast Efficient Score Calculation
# compute_efficient_score_and_update <- function(X, Y, beta, b, h) {
#   n <- nrow(X)
#   p <- ncol(X)
#   d <- ncol(beta)
#   betaX <- X %*% beta
#   
#   m_results <- estimate_m_and_gradient(X, Y, beta, b)
#   residuals <- Y - m_results$m_est
#   sigma2 <- estimate_sigma2(X, residuals, h)
#   
#   ratio_term <- matrix(0, nrow = n, ncol = p)
#   for (i in 1:n) {
#     beta_x_i <- betaX[i, ]
#     diff_betaX <- sweep(betaX, 2, beta_x_i, "-")
#     weights <- kernel_weight_vec(diff_betaX, b, d)
#     
#     combined_weights <- weights / sigma2
#     num <- colSums(sweep(X, 1, combined_weights, "*"))
#     den <- sum(weights * (1 / sigma2)) + 1e-10
#     ratio_term[i, ] <- num / den
#   }
#   
#   S_eff_list <- vector("list", n)
#   for (i in 1:n) {
#     term1 <- (X[i, ] - ratio_term[i, ]) / sigma2[i] 
#     term2 <- m_results$m_prime_est[i, ]             
#     term3 <- residuals[i]                           
#     S_eff_list[[i]] <- term1 %*% t(term2) * term3 
#   }
#   
#   vec_S <- sapply(S_eff_list, as.vector) 
#   mean_vec_S <- rowMeans(vec_S) 
#   
#   # Fast vectorized outer product
#   avg_outer_S <- (vec_S %*% t(vec_S)) / n
#   
#   vec_beta_k <- as.vector(beta)
#   # Standard Newton step
#   vec_beta_next <- vec_beta_k + MASS::ginv(avg_outer_S) %*% mean_vec_S
#   
#   matrix(vec_beta_next, nrow = p, ncol = d)
# }