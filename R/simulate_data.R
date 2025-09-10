simulate_data <- function(n,
                          p = 20, 
                          q = 10,
                          rho = 0.7) {
  # Confounders [C] ---------------------------------------------------------
  Sigma_C <- toeplitz(c(1, .5, rep(0, q-2)))
  C       <- MASS::mvrnorm(n = n, mu = rep(0, q), Sigma = Sigma_C)

  # Exposures [X | C] -------------------------------------------------------
  Sigma_X <- toeplitz(rho^(0:((p)-1)))
  Theta   <- matrix(1 / (1:q)^2, nrow = q, ncol = p)
  X_mean  <- pnorm(C %*% Theta) + 0.75*rnorm(n) - 0.5
  
  X <- t(apply(X_mean, 1, function(mu_i)
    MASS::mvrnorm(1, mu = mu_i, Sigma = Sigma_X)
  ))  # n x (p-4)
  
  # 2 dimensions, first direction is 4, second is 2. 
  beta <- matrix(c(rep(1 / sqrt(2), 2),rep(0,18), 
                 c(rep(0,16), rep(1 / sqrt(4), 4)) ),nrow = p, ncol = d)
  # the unique Projection matrix from R^p to R^d
  P_beta <- beta %*% solve(t(beta) %*% beta) %*% t(beta)
  
  # --- Outcome ---
  Z   <- X %*% beta                   # n x 2
  g_C <- C %*% Theta[,1]  
  h_Z <- Z[, 1] + 0.1*Z[,2]^2
  
  eps_Y <- rnorm(n)
  
  Y     <- h_Z + g_C + 0.1* Z[,1]*C[,1] + eps_Y * (sqrt( 0.5 + pnorm(rowSums(C[,1:2])) )) 
  
  list(C = C, X = X, Y = Y, beta = beta, P_beta= P_beta)
}
