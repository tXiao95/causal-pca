simulate_data <- function(n,
                          p = 20, 
                          q = 10,
                          rho = 0.8,
                          Z1_coef = 10,
                          Z2_coef = 2,
                          Z12_coef = 1,
                          h_Z_coef = 0.1,
                          g_C_coef = 1,
                          var_scale = 3,
                          interaction_coef = 1) {
  if(p < 6){
    stop("Number of exposure variables must be greater than or equal to 6")
  }
  # Confounders [C] ---------------------------------------------------------
  Sigma_C <- toeplitz(c(1, .5, rep(0, q-2)))
  C       <- MASS::mvrnorm(n = n, mu = rep(0, q), Sigma = Sigma_C)

  # Inflate variance in X's that don't feed into beta to distract PCA 
  idx_beta  <- c(1:2, 17:20)
  scale_vec <- rep(1, p)
  scale_vec[-idx_beta] <- var_scale
  D         <- diag(scale_vec, p, p)
  
  # Exposures [X | C] -------------------------------------------------------
  Sigma_X <- D %*% toeplitz(rho^(0:((p)-1))) %*% D
  Theta   <- matrix(  ( 1:p / (1:q)^2 ), nrow = q, ncol = p)
  X_mean  <- pnorm(C %*% Theta) 
  #X_mean  <- (C %*% Theta) 
  
  X <- t(apply(X_mean, 1, function(mu_i)
    MASS::mvrnorm(1, mu = mu_i, Sigma = Sigma_X)
  ))  # n x (p-4)
  
  # 2 dimensions, first direction is 4, second is 2. 
  beta <- matrix(c(rep(1 / sqrt(2), 2),rep(0, p - 2), 
                 c(rep(0, p - 4), rep(1 / sqrt(4), 4)) ),nrow = p, ncol = 2)
  # the unique Projection matrix from R^p to R^d
  P_beta <- beta %*% solve(t(beta) %*% beta) %*% t(beta)
  
  # Dimension reduction
  Z   <- X %*% beta                   # n x 2
  
  # --- Outcome ---
  g_C   <- sin(C %*% Theta[,1]) 
  h_Z   <- Z1_coef * Z[,1] + Z2_coef * Z[,2]^2 + Z12_coef * Z[,1]*Z[,2]
  eps_Y <- rnorm(n, mean = 0, sd = 0.5)
  
  # Y response model
  Y     <- h_Z_coef * h_Z + g_C_coef * g_C + interaction_coef * (Z[,1]*C[,1]) + eps_Y * (sqrt( 0.5 + pnorm(rowSums(C[,1:2])) )) 
  
  # Causal mean (g_C is mean zero)
  mu_X <- h_Z_coef * h_Z
  
  list(Y = Y, C = C, X = X, mu_X = mu_X, beta = beta, P_beta= P_beta)
}
