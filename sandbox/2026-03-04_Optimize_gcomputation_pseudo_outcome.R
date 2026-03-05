library(data.table)
library(mvtnorm)
source("R/nuisance_outcome_regression.R")
source("R/nuisance_gps.R")

# Can maybe apply this same fix to gcomputation for the ERS estimation

# Functions ---------------------------------------------------------------
estimate_pseudo_outcomes_old <- function(Y, X, C, out_model, gps_model, delta_n = 1e-16) {

  # ---------------------------------------------------------
  # Input Validation & Safe Column Naming
  # ---------------------------------------------------------
  # Capture original names BEFORE coercion
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL

  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  n <- length(Y)

  # Safely apply names only if the ORIGINAL input lacked them
  if (is.null(orig_X_names) && !is.null(out_model$X_names)) colnames(X_df) <- out_model$X_names
  if (is.null(orig_C_names) && !is.null(out_model$C_names)) colnames(C_df) <- out_model$C_names

  # ---------------------------------------------------------
  # Pre-computations (Batch predicting observed data)
  # ---------------------------------------------------------
  df_observed <- cbind(X_df, C_df)

  # Pre-calculate the "diagonal" elements: m(X_j, C_j) and pi(X_j | C_j)
  # Doing this outside the loop saves overhead and is much safer
  m_obs  <- predict(out_model, newdata = df_observed)
  pi_obs <- predict(gps_model, newdata = df_observed)
  pi_obs <- pmax(pi_obs, delta_n) # Apply safety flooring to the denominator

  # ---------------------------------------------------------
  # Main Loop over Individuals
  # ---------------------------------------------------------
  pseudo_outcomes <- vapply(1:n, function(j) {
    message(j)
    # Create a grid where individual j's treatment is repeated n times,
    # paired with EVERY individual's confounders (C_1 to C_n)
    X_j_rep <- X_df[rep(j, n), , drop = FALSE]
    df_grid <- cbind(X_j_rep, C_df)

    # Predict m(X_j, C_i) and pi(X_j | C_i) across all i = 1...n
    m_grid  <- predict(out_model, newdata = df_grid)
    pi_grid <- predict(gps_model, newdata = df_grid)

    # Calculate the empirical expectations (marginalized over C)
    mean_pi <- mean(pi_grid)
    mean_m  <- mean(m_grid)

    # Assemble the pseudo-outcome using the pre-computed observed values
    xi_j <- ((Y[j] - m_obs[j]) / pi_obs[j]) * mean_pi + mean_m

    return(xi_j)

  }, numeric(1L))

  return(pseudo_outcomes)
}

estimate_pseudo_outcomes_new <- function(Y, X, C, out_model, gps_model, delta_n = 1e-16) {
  
  # ---------------------------------------------------------
  # Input Validation & Safe Column Naming
  # ---------------------------------------------------------
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL
  
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  n <- length(Y)
  
  if (is.null(orig_X_names) && !is.null(out_model$X_names)) colnames(X_df) <- out_model$X_names
  if (is.null(orig_C_names) && !is.null(out_model$C_names)) colnames(C_df) <- out_model$C_names
  
  # ---------------------------------------------------------
  # Pre-computations (Batch predicting observed data)
  # ---------------------------------------------------------
  df_observed <- cbind(X_df, C_df)
  
  m_obs  <- predict(out_model, newdata = df_observed)
  pi_obs <- predict(gps_model, newdata = df_observed)
  pi_obs <- pmax(pi_obs, delta_n) 
  
  # ---------------------------------------------------------
  # Zero-Copy Evaluation Grid Setup
  # ---------------------------------------------------------
  
  # 1. Base the grid entirely on C (which never changes for the marginalization)
  dt_grid <- as.data.table(C_df)
  
  # 2. Pre-allocate the X columns with dummy values (0.0)
  x_names <- colnames(X_df)
  for (x_col in x_names) {
    dt_grid[, (x_col) := 0.0]
  }
  
  # 3. Force the exact column order the models expect
  setcolorder(dt_grid, c(x_names, colnames(C_df)))
  
  # ---------------------------------------------------------
  # Main Loop over Individuals (Hyper-Optimized)
  # ---------------------------------------------------------
  pseudo_outcomes <- vapply(1:n, function(j) {
    message(j)
    
    # ZERO-COPY IN-PLACE UPDATE:
    # Instead of `rep()` and `cbind()`, we instantly overwrite the memory 
    # of the X columns with the values for individual 'j'.
    for (x_col in x_names) {
      data.table::set(dt_grid, j = x_col, value = X_df[[x_col]][j])
    }
    
    # Predict 
    m_grid  <- predict(out_model, newdata = dt_grid)
    pi_grid <- predict(gps_model, newdata = dt_grid)
    
    # Calculate the empirical expectations
    mean_pi <- mean(pi_grid)
    mean_m  <- mean(m_grid)
    
    # Assemble the pseudo-outcome
    xi_j <- ((Y[j] - m_obs[j]) / pi_obs[j]) * mean_pi + mean_m
    
    return(xi_j)
    
  }, numeric(1L))
  
  return(pseudo_outcomes)
}

# ---------------------------------------------------------
# 1. Setup Data
# ---------------------------------------------------------
set.seed(123)
n <- 5000
p <- 2
q <- 2

# Confounders
C <- matrix(rnorm(n * q), n, q)
colnames(C) <- c("C1", "C2")

# Treatments (Multivariate Normal given C)
beta_X <- matrix(c(0.5, -0.2, 
                   -0.8,  0.4), nrow = q, ncol = p)
X <- C %*% beta_X + rmvnorm(n, sigma = diag(p))
colnames(X) <- c("X1", "X2")

# Outcome (Y = C1 - C2 + X1 - X2 + noise)
Y <- as.numeric(C %*% c(0.5, -0.5) + X %*% c(1, -1) + rnorm(n, sd = 0.5))

# ---------------------------------------------------------
# 2. Fit Actual Nuisance Models
# ---------------------------------------------------------

# Fit Outcome Regression E[Y | X, C]
# Using SL.glm to keep the test instantaneous
suppressWarnings({
  out_mod <- outcome_model(Y = Y, 
                           X = X, 
                           C = C, 
                           mu_fitter = SL_outcome_fitter, 
                           SL.lib = "SL.glm")
})

# Fit Propensity Score Density f(X | C) using the new gps_model gateway
gps_mod <- gps_model(X = X, 
                     C = C, 
                     pi_fitter = mvn_fitter, 
                     method = "linear")

# ---------------------------------------------------------
# 3. Compute Pseudo-Outcomes
# ---------------------------------------------------------
time_a <- system.time({
  pseudo_outcomes <- estimate_pseudo_outcomes_old(Y = Y, 
                                              X = X, 
                                              C = C, 
                                              out_model = out_mod, 
                                              gps_model = gps_mod)
  
})
time_b <- system.time({
  pseudo_outcomes <- estimate_pseudo_outcomes_new(Y = Y, 
                                              X = X, 
                                              C = C, 
                                              out_model = out_mod, 
                                              gps_model = gps_mod)
  
})

message("Without data.table: ", time_a["elapsed"])
message("With data.table: ", time_b["elapsed"])

# > message("Without data.table: ", time_a["elapsed"])
# Without data.table: 44.45
# > message("With data.table: ", time_b["elapsed"])
# With data.table: 18.108
