simulate_data <- function(n,
                          p = 20, 
                          q = 10,
                          rho = 0.8,
                          Z1_coef = 1,
                          Z2_coef = 0.5,
                          Z12_coef = 0.5,
                          h_Z_coef = 0.1,
                          g_C_coef = 0.1,
                          var_scale = 3,
                          interaction_coef = 0.1) {
  if(p < 6){
    stop("Number of exposure variables must be greater than or equal to 6")
  }
  # Confounders [C] ---------------------------------------------------------
  Sigma_C <- toeplitz(c(1, .5, rep(0, q-2)))
  C       <- MASS::mvrnorm(n = n, mu = rep(0, q), Sigma = Sigma_C)

  # Inflate variance in X's that don't feed into beta to distract PCA 
  idx_beta  <- c(1:2, 7:10)
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
  g_C   <- sin( C %*% Theta[,1] ) 
  h_Z   <- Z1_coef * Z[,1] + Z2_coef * Z[,2]^2 + Z12_coef * Z[,1]*Z[,2]
  eps_Y <- rnorm(n, mean = 0, sd = 0.5)
  
  # Y response model
  Y     <- h_Z_coef * h_Z + g_C_coef * g_C + interaction_coef * (Z[,1]*C[,1]) + eps_Y * (sqrt( 0.5 + pnorm(rowSums(C[,1:2])) )) 
  
  # Causal mean at observed points (g_C is mean zero)
  mu_X <- h_Z_coef * h_Z
  
  # Causal mean function
  h_Z_fun <- function(Z1, Z2){
    return( h_Z_coef * (Z1_coef * Z1 + Z2_coef * Z2^2 + Z12_coef * Z1 * Z2 ) )
  }
  
  list(Y = Y, C = C, X = X, mu_X = mu_X, beta_true = beta, P_beta= P_beta, h_Z_fun = h_Z_fun, d = ncol(beta))
}

simulate_weak_dim_signal <- function(n = 1000,
                                     p = 10,
                                     q = 5,
                                     rho_X = 0.3,
                                     confounding_strength = 3,
                                     spurious_strength = 6,
                                     var_scale = 5,             # <-- New: Inflates variance for PCA trap
                                     interaction_coef = 1,
                                     sigma2 = 0.25,
                                     heteroskedastic = FALSE) {

  if (p < 7) stop("p must be >= 7 to accommodate causal, spurious, and noise directions.")
  if (q < 2) stop("q must be >= 2 for the confounding structure.")

  # 1. Confounders [C] ------------------------------------------------------
  C <- matrix(rnorm(n * q), nrow = n, ncol = q)

  # 2. Structural Directions in X -------------------------------------------
  # True causal directions (beta): depends on X1-X4
  beta <- matrix(0, nrow = p, ncol = 2)
  beta[1:2, 1] <- 1 / sqrt(2)
  beta[3:4, 2] <- 1 / sqrt(2)

  # Spurious direction (gamma): depends on X5-X6
  gamma <- rep(0, p)
  gamma[5:6] <- 1 / sqrt(2)

  # 3. Exposures [X | C] ----------------------------------------------------
  # PCA Trap: Inflate variance in X's that don't feed into beta or gamma
  idx_signal <- 1:6
  scale_vec <- rep(1, p)
  scale_vec[-idx_signal] <- var_scale
  D <- diag(scale_vec, p, p)

  # Apply variance inflation to the covariance matrix
  Sigma_X <- D %*% toeplitz(rho_X^(0:(p-1))) %*% D

  # Mean of X given C
  mu_X <- matrix(0, nrow = n, ncol = p)
  # C2 confounds the causal variables (X1-X4)
  mu_X[, 1:4] <- C[, 2] * confounding_strength
  # C1 heavily drives the spurious variables (X5, X6)
  mu_X[, 5:6] <- C[, 1] * spurious_strength

  # Generate X
  X <- mu_X + MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma_X)

  # 4. Projections ----------------------------------------------------------
  Z <- X %*% beta          # n x 2 (Causal dimensions)

  # 5. Causal Mean Function h(Z) --------------------------------------------
  h_Z_fun <- function(Z1, Z2) {
    return( 3 * tanh(Z1) + 8 * (pnorm(2 * Z2) - 0.5) + 0.2 * Z2*Z1 )
  }

  mu_causal <- h_Z_fun(Z[, 1], Z[, 2])

  # 6. Outcome Y | X, C -----------------------------------------------------
  # Spurious association: Y strongly depends on C1.
  g_C <- 5 * C[, 1] + 2 * C[, 2]

  # Interaction
  C_interact <- if(q >= 3) C[, 3] else C[, 1]
  interaction_term <- interaction_coef * Z[, 1] * C_interact

  # Variance of Y
  if (heteroskedastic) {
    sigma_Y <- sigma2 * exp(0.2 * Z[, 1] + 0.2 * C[, 2])
  } else {
    sigma_Y <- sigma2 * rep(1, n)
  }

  eps_Y <- rnorm(n, mean = 0, sd = sigma_Y)

  # Final Response
  Y <- mu_causal + g_C + interaction_term + eps_Y

  # True Projection Matrix
  P_beta <- beta %*% solve(t(beta) %*% beta) %*% t(beta)

  list(Y = Y,
       C = C,
       X = X,
       Z = Z,
       mu_X = mu_causal,
       beta_true = beta,
       P_beta = P_beta,
       h_Z_fun = h_Z_fun,
       d = ncol(beta))
}

simulate_causal_sdr <- function(n = 1000,
                                p = 10,
                                q = 5,
                                rho_X = 0.7,
                                causal_conf_strength = 1.0,   # Keeps Z bounded in [-3, 3] grid
                                spurious_strength = 2.0,      # Keeps the MAVE trap strong
                                var_scale = 5,                # Keeps the PCA trap strong
                                signal_multiplier = 2.0,      # Knob 1: ERS Signal strength
                                noise_sd = 0.5,               # Knob 2: Noise strength (SNR control)
                                interaction_coef = 10,
                                heteroskedastic = FALSE) {
  
  if (p < 7) stop("p must be >= 7 to accommodate causal, spurious, and noise directions.")
  if (q < 2) stop("q must be >= 2 for the confounding structure.")
  
  # 1. Confounders [C] ------------------------------------------------------
  C <- matrix(rnorm(n * q), nrow = n, ncol = q)
  
  # 2. Structural Directions in X -------------------------------------------
  # True causal directions (beta): depends on X1-X4
  beta <- matrix(0, nrow = p, ncol = 2)
  beta[1:2, 1] <- 1 / sqrt(2)
  beta[3:4, 2] <- 1 / sqrt(2)
  
  # Spurious direction (gamma): depends on X5-X6
  gamma <- rep(0, p)
  gamma[5:6] <- 1 / sqrt(2)
  
  # 3. Exposures [X | C] (Linear to preserve DGP requirements) --------------
  # PCA Trap: Inflate variance in noise X's
  idx_signal <- 1:6
  scale_vec <- rep(1, p)
  scale_vec[-idx_signal] <- var_scale
  D <- diag(scale_vec, p, p)
  
  Sigma_X <- D %*% toeplitz(rho_X^(0:(p-1))) %*% D
  
  # Mean of X given C
  mu_X <- matrix(0, nrow = n, ncol = p)
  # C2 confounds the causal variables lightly
  mu_X[, 1:4] <- C[, 2] * causal_conf_strength
  
  # NEW: Create TWO independent spurious directions
  mu_X[, 5] <- C[, 1] * spurious_strength 
  mu_X[, 6] <- C[, 3] * spurious_strength
  
  X <- mu_X + MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma_X)
  
  # 4. Projections ----------------------------------------------------------
  Z <- X %*% beta          # n x 2 
  
  # 5. Causal Mean Function h(Z) --------------------------------------------
  # Highly distinct: Z1 is an odd function (Sine), Z2 is an even function (Quadratic)
  h_Z_fun <- function(Z1, Z2) {
    return( signal_multiplier * ( 2 * sin(Z1) + 0.5 * Z2^2 ) )
  }
  
  mu_causal <- h_Z_fun(Z[, 1], Z[, 2])
  
  # 6. Outcome Y | X, C -----------------------------------------------------
  # Spurious association relies heavily on BOTH C1 and C3
  # By increasing the multipliers to 10, we ensure the spurious variance 
  # completely dominates the causal signal, baiting MAVE perfectly.
  g_C <- 10 * sin(1.5 * C[, 1]) + 10 * (C[, 3]^2 - 1)
  
  # Interaction
  C_interact <- if(q >= 3) C[, 3] else C[, 1]
  #interaction_term <- interaction_coef * Z[, 1] * C_interact
  interaction_term <- interaction_coef * ( X[,5] + X[,6] + Z[,1] ) * C_interact
  
  # Variance of Y 
  if (heteroskedastic) {
    # Variance grows with Z1 and C2. 
    sigma_Y <- noise_sd * exp(0.2 * Z[, 1] + 0.2 * C[, 2])
  } else {
    sigma_Y <- noise_sd * rep(1, n)
  }
  
  eps_Y <- rnorm(n, mean = 0, sd = sigma_Y)
  
  # Final Response
  Y <- mu_causal + g_C + interaction_term + eps_Y
  
  # True Projection Matrix
  P_beta <- beta %*% solve(t(beta) %*% beta) %*% t(beta)
  
  list(Y = Y, 
       C = C, 
       X = X, 
       Z = Z, 
       mu_X = mu_causal, 
       beta_true = beta, 
       P_beta = P_beta, 
       h_Z_fun = h_Z_fun, 
       d = ncol(beta))
}