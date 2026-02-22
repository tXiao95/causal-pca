
# Luo and Cai (2016) ------------------------------------------------------

sim_luo_cai_mod1 <- function(n = 200, p = 10, seed = 123) {
  set.seed(seed)
  d <- 1
  
  # Predictors: Independent Uniform(-2, 2)
  X <- matrix(runif(n * p, -2, 2), n, p)
  
  # True Basis (m1 depends on the sum of the first 4 components)
  beta_true <- matrix(0, p, d)
  beta_true[1:4, 1] <- 1
  beta_true <- qr.Q(qr(beta_true))
  
  # Conditional Mean m1(X)
  sum_x4 <- rowSums(X[, 1:4, drop = FALSE])
  m_1 <- sin(0.3 * sum_x4)
  
  # Conditional Variance Function
  # The multiplier is strictly 0.5. Var(eps_U) is 1/3.
  sigma2_true <- rep((0.5^2) * (1/3), n)
  
  # Error Term: eps_U is Uniform(-1, 1)
  eps <- runif(n, -1, 1)
  Y <- m_1 + 0.5 * eps
  
  # Return standardized list
  list(
    X = X, Y = Y, beta_true = beta_true, sigma2_true = sigma2_true, 
    d = d, n = n, p = p, paper = "Luo and Cai (2016) - Model I"
  )
}

sim_luo_cai_mod2 <- function(n = 200, p = 10, seed = 123) {
  set.seed(seed)
  d <- 1
  
  # Predictors: Normal with AR(1) covariance 0.5^|i-j|
  Sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  
  # True Basis
  beta_true <- matrix(0, p, d)
  beta_true[1:4, 1] <- 1
  beta_true <- qr.Q(qr(beta_true))
  
  # Conditional Mean m1(X)
  sum_x4 <- rowSums(X[, 1:4, drop = FALSE])
  m_1 <- sin(0.3 * sum_x4)
  
  # Conditional Variance Function
  # The multiplier relies on m1. Var(eps_N) is 0.5^2 = 0.25.
  sigma_mult <- 0.1 + (m_1^2) / 5
  sigma2_true <- (sigma_mult^2) * (0.5^2)
  
  # Error Term: eps_N is Normal(0, 0.5^2)
  eps <- rnorm(n, mean = 0, sd = 0.5)
  Y <- m_1 + sigma_mult * eps
  
  # Return standardized list
  list(
    X = X, Y = Y, beta_true = beta_true, sigma2_true = sigma2_true, 
    d = d, n = n, p = p, paper = "Luo and Cai (2016) - Model II"
  )
}

sim_luo_cai_mod3 <- function(n = 200, p = 10, seed = 123) {
  set.seed(seed)
  # Default matches Luo and Cai (2016) Model III settings
  d <- 2       
  
  # Predictors: Independent Uniform(-2, 2)
  X <- matrix(runif(n * p, -2, 2), n, p)
  
  # True Basis
  beta_true <- matrix(0, p, d)
  beta_true[1, 1] <- 1
  beta_true[2, 2] <- 1
  beta_true <- qr.Q(qr(beta_true))
  
  # Conditional Mean
  m_2 <- X[, 1] / (0.5 + (1.5 + X[, 2])^2)
  
  # Conditional Variance Function
  sigma_X_true <- sqrt((X[, 5]^2) * (m_2^2) + 0.01)
  
  # Error Term: Uniform(-1, 1) has variance 1/3
  eps <- runif(n, -1, 1)
  Y <- m_2 + sigma_X_true * eps
  
  # Return standardized list
  list(
    X = X,
    Y = Y,
    beta_true = beta_true,
    sigma2_true = (sigma_X_true^2) * (1/3), 
    d = d,
    n = n,
    p = p,
    paper = "Luo and Cai (2016) - Model III"
  )
}

sim_luo_cai_mod4 <- function(n = 200, p = 10, seed = 123) {
  set.seed(seed)
  d <- 2
  
  # Predictors: Independent Uniform(-2, 2)
  X <- matrix(runif(n * p, -2, 2), n, p)
  
  # True Basis (m2 depends only on X1 and X2)
  beta_true <- matrix(0, p, d)
  beta_true[1, 1] <- 1
  beta_true[2, 2] <- 1
  beta_true <- qr.Q(qr(beta_true))
  
  # Conditional Mean m2(X)
  m_2 <- X[, 1] / (0.5 + (1.5 + X[, 2])^2)
  
  # Conditional Variance Function
  # The multiplier relies on X5. Var(eps_N) is 0.5^2 = 0.25.
  sigma_mult <- exp(X[, 5])
  sigma2_true <- (sigma_mult^2) * (0.5^2)
  
  # Error Term: eps_N is Normal(0, 0.5^2)
  eps <- rnorm(n, mean = 0, sd = 0.5)
  Y <- m_2 + sigma_mult * eps
  
  # Return standardized list
  list(
    X = X, Y = Y, beta_true = beta_true, sigma2_true = sigma2_true, 
    d = d, n = n, p = p, paper = "Luo and Cai (2016) - Model IV"
  )
}

sim_luo_cai_mod5 <- function(n = 200, p = 10, seed = 123) {
  set.seed(seed)
  d <- 2
  
  # Predictors: Normal with AR(1) covariance 0.5^|i-j|
  Sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  
  # True Basis (m2 depends only on X1 and X2)
  beta_true <- matrix(0, p, d)
  beta_true[1, 1] <- 1
  beta_true[2, 2] <- 1
  beta_true <- qr.Q(qr(beta_true))
  
  # Conditional Mean m2(X)
  m_2 <- X[, 1] / (0.5 + (1.5 + X[, 2])^2)
  
  # Conditional Variance Function
  # The multiplier is a step function based on the squared Euclidean norm of X.
  norm_sq_X <- rowSums(X^2)
  sigma_mult <- ifelse(norm_sq_X < p, 1/3, 3)
  
  # Var(eps_U) is 1/3.
  sigma2_true <- (sigma_mult^2) * (1/3)
  
  # Error Term: eps_U is Uniform(-1, 1)
  eps <- runif(n, -1, 1)
  Y <- m_2 + sigma_mult * eps
  
  # Return standardized list
  list(
    X = X, Y = Y, beta_true = beta_true, sigma2_true = sigma2_true, 
    d = d, n = n, p = p, paper = "Luo and Cai (2016) - Model V"
  )
}

# Xia et al. (2002) -------------------------------------------------------
sim_xia_2002_ex2 <- function(n = 200, p = 10, seed = 123) {
  set.seed(seed)
  d <- 4
  
  # Predictors: Independent Standard Normal N(0, 1)
  X <- matrix(rnorm(n * p), n, p)
  
  # Explicit true coefficients from the text
  b1 <- c(1, 2, 3, 4, 0, 0, 0, 0, 0, 0) / sqrt(30)
  b2 <- c(-2, 1, -4, 3, 1, 2, 0, 0, 0, 0) / sqrt(35)
  b3 <- c(0, 0, 0, 0, 2, -1, 2, 1, 2, 1) / sqrt(15)
  b4 <- c(0, 0, 0, 0, 0, 0, -1, -1, 1, 1) / 2
  
  # True Basis Matrix B0
  beta_true <- cbind(b1, b2, b3, b4)
  beta_true <- qr.Q(qr(beta_true))
  
  # Project X onto the true directions to calculate the mean
  X_b1 <- as.vector(X %*% b1)
  X_b2 <- as.vector(X %*% b2)
  X_b3 <- as.vector(X %*% b3)
  X_b4 <- as.vector(X %*% b4)
  
  # Conditional Mean: y = (X'b1)(X'b2)^2 + (X'b3)(X'b4)
  m_X <- X_b1 * (X_b2^2) + (X_b3 * X_b4)
  
  # Conditional Variance: (0.5^2) * Var(eps_N) = 0.25
  sigma2_true <- rep(0.25, n)
  
  # Error Term: eps_N is Standard Normal N(0, 1)
  eps <- rnorm(n, mean = 0, sd = 1)
  Y <- m_X + 0.5 * eps
  
  list(
    X = X, Y = Y, beta_true = beta_true, sigma2_true = sigma2_true, 
    d = d, n = n, p = p, paper = "Xia et al. (2002) - Ex 2"
  )
}

# Nabi et al. (2022) ------------------------------------------------------

sim_nabi_case1_p6 <- function(n = 200, seed = 123) {
  set.seed(seed)
  p <- 6
  d <- 2
  
  # Confounders: Standard Multivariate Normal
  C <- matrix(rnorm(n * 4), n, 4)
  
  X <- matrix(0, n, p)
  
  # X1 and X2 generated from MVN with specific means and AR(1) covariance
  mu_12 <- cbind(C[, 1], -C[, 1])
  Sigma_12 <- 0.5^abs(outer(1:2, 1:2, "-"))
  X[, 1:2] <- mu_12 + MASS::mvrnorm(n, mu = rep(0, 2), Sigma = Sigma_12)
  
  # X3 and X4 are heteroscedastic normal distributions
  X[, 3] <- rnorm(n, mean = abs(X[, 1] + X[, 2]), sd = sqrt(abs(X[, 1])))
  X[, 4] <- rnorm(n, mean = sqrt(abs(X[, 1] + X[, 2])), sd = sqrt(abs(X[, 2])))
  
  # X5 and X6 are Bernoulli (Discrete)
  X[, 5] <- rbinom(n, 1, prob = exp(X[, 2]) / (1 + exp(X[, 2])))
  X[, 6] <- rbinom(n, 1, prob = pnorm(X[, 2]))
  
  # True Basis
  beta_1 <- rep(1 / sqrt(6), 6)
  beta_2 <- c(1, -1, 1, -1, 1, -1) / sqrt(6)
  beta_true <- qr.Q(qr(cbind(beta_1, beta_2)))
  
  # Outcome Y
  sum_C <- rowSums(C)
  sum_X <- rowSums(X)
  m_X <- as.vector(X %*% beta_1 + (X %*% beta_2)^2 + sum_C + sum_X * sum_C)
  eps <- rnorm(n, mean = 0, sd = 1)
  Y <- m_X + eps
  
  list(
    X = X, C = C, Y = Y, beta_true = beta_true, sigma2_true = rep(1, n),
    d = d, n = n, p = p, paper = "Nabi et al. (2022) - Case 1, p=6"
  )
}

sim_nabi_case1_p12 <- function(n = 200, seed = 123) {
  set.seed(seed)
  p <- 12
  d <- 2
  
  C <- matrix(rnorm(n * 4), n, 4)
  X <- matrix(0, n, p)
  
  # X1, X2, and X7 to X12 are generated from a (p-4) dimensional MVN
  mu_8 <- cbind(
    C[, 1], -C[, 1],                     # X1, X2
    C[, 1], C[, 2], C[, 3],              # X7, X8, X9
    -C[, 1] + C[, 2],                    # X10
    -C[, 2] + C[, 3],                    # X11
    -C[, 3] + C[, 4]                     # X12
  )
  Sigma_8 <- 0.5^abs(outer(1:8, 1:8, "-"))
  err_8 <- MASS::mvrnorm(n, mu = rep(0, 8), Sigma = Sigma_8)
  
  X[, 1:2] <- mu_8[, 1:2] + err_8[, 1:2]
  X[, 7:12] <- mu_8[, 3:8] + err_8[, 3:8]
  
  # Non-linear and discrete variables remain the same structure based on X1 and X2
  X[, 3] <- rnorm(n, mean = abs(X[, 1] + X[, 2]), sd = sqrt(abs(X[, 1])))
  X[, 4] <- rnorm(n, mean = sqrt(abs(X[, 1] + X[, 2])), sd = sqrt(abs(X[, 2])))
  X[, 5] <- rbinom(n, 1, prob = exp(X[, 2]) / (1 + exp(X[, 2])))
  X[, 6] <- rbinom(n, 1, prob = pnorm(X[, 2]))
  
  # True Basis (Last 6 components are zero)
  beta_1 <- c(rep(1 / sqrt(6), 6), rep(0, 6))
  beta_2 <- c(c(1, -1, 1, -1, 1, -1) / sqrt(6), rep(0, 6))
  beta_true <- qr.Q(qr(cbind(beta_1, beta_2)))
  
  sum_C <- rowSums(C)
  sum_X <- rowSums(X)
  m_X <- as.vector(X %*% beta_1 + (X %*% beta_2)^2 + sum_C + sum_X * sum_C)
  eps <- rnorm(n, mean = 0, sd = 1)
  Y <- m_X + eps
  
  list(
    X = X, C = C, Y = Y, beta_true = beta_true, sigma2_true = rep(1, n),
    d = d, n = n, p = p, paper = "Nabi et al. (2022) - Case 1, p=12"
  )
}

sim_nabi_case2_p6 <- function(n = 200, seed = 123) {
  set.seed(seed)
  p <- 6
  d <- 2
  
  C <- matrix(rnorm(n * 4), n, 4)
  sum_C <- rowSums(C)
  
  # All predictors generated from MVN
  mu <- cbind(
    sum_C,
    -C[, 1] + C[, 2] - C[, 3] + C[, 4],   # Sum_i (-1)^i C_i
    C[, 1] - C[, 2] + C[, 3] - C[, 4],
    -C[, 1] + C[, 2] + C[, 3] - C[, 4],
    sum_C - 2 * C[, 3],
    sum_C - 2 * C[, 1]
  )
  
  Sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
  X <- mu + MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  
  # True Basis
  beta_1 <- rep(1 / sqrt(6), 6)
  beta_2 <- c(1, -1, 1, -1, 1, -1) / sqrt(6)
  beta_true <- qr.Q(qr(cbind(beta_1, beta_2)))
  
  sum_X <- rowSums(X)
  m_X <- as.vector(X %*% beta_1 + (X %*% beta_2)^2 + sum_C + sum_X * sum_C)
  eps <- rnorm(n, mean = 0, sd = 1)
  Y <- m_X + eps
  
  list(
    X = X, C = C, Y = Y, beta_true = beta_true, sigma2_true = rep(1, n),
    d = d, n = n, p = p, paper = "Nabi et al. (2022) - Case 2, p=6"
  )
}

sim_nabi_case2_p12 <- function(n = 200, seed = 123) {
  set.seed(seed)
  p <- 12
  d <- 2
  
  C <- matrix(rnorm(n * 4), n, 4)
  sum_C <- rowSums(C)
  
  #X1 to X12
  mu <- cbind(
    sum_C,
    -C[, 1] + C[, 2] - C[, 3] + C[, 4], 
    C[, 1] - C[, 2] + C[, 3] - C[, 4],
    -C[, 1] + C[, 2] + C[, 3] - C[, 4],
    sum_C - 2 * C[, 3],
    sum_C - 2 * C[, 1],
    C[, 1], C[, 2], C[, 3],               # X7 to X9
    -C[, 1], -C[, 2], -C[, 3]             # X10 to X12
  )
  
  Sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
  X <- mu + MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  
  # True Basis (Last 6 components are zero)
  beta_1 <- c(rep(1 / sqrt(6), 6), rep(0, 6))
  beta_2 <- c(c(1, -1, 1, -1, 1, -1) / sqrt(6), rep(0, 6))
  beta_true <- qr.Q(qr(cbind(beta_1, beta_2)))
  
  sum_X <- rowSums(X)
  m_X <- as.vector(X %*% beta_1 + (X %*% beta_2)^2 + sum_C + sum_X * sum_C)
  eps <- rnorm(n, mean = 0, sd = 1)
  Y <- m_X + eps
  
  list(
    X = X, C = C, Y = Y, beta_true = beta_true, sigma2_true = rep(1, n),
    d = d, n = n, p = p, paper = "Nabi et al. (2022) - Case 2, p=12"
  )
}