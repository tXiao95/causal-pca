# =====================================================================
# 1. KERNEL FUNCTION
# Reference: Section 4 (Monte Carlo Simulations)
# =====================================================================

# The paper specifies a 4th-order kernel function, which is necessary to 
# achieve the required convergence rates for structural dimensions d_0 <= 6.
K4 <- function(u) {
  # Equation: K_4(u) = (105/64)(1 - 3u^2)(1 - u^2)^2 * 1(|u| <= 1)
  val <- (105/64) * (1 - 3 * u^2) * (1 - u^2)^2
  
  # Apply the indicator function 1(|u| <= 1)
  val[abs(u) > 1] <- 0
  return(val)
}

# =====================================================================
# 2. CROSS-VALIDATION OBJECTIVE FUNCTION (CV_M)
# Reference: Remark 3, Equations (2.11) and (2.12)
# =====================================================================

cv_objective_CMS <- function(params, X, Y, d, p, n) {
  # The optimizer passes a flattened vector of parameters.
  # We extract the (p-d)*d elements for C_Md and the d elements for bandwidths h_d.
  len_C <- (p - d) * d
  C_vec <- params[1:len_C]
  h <- params[(len_C + 1):length(params)]
  
  # Enforce strictly positive bandwidths
  if (any(h <= 0)) return(Inf)
  
  # Parameterize the basis matrix B_Md = (I_d, C_Md^T)^T to avoid identifiability 
  # problems using the Grassmann manifold local coordinate system.
  I_d <- diag(d)
  C_Md <- matrix(C_vec, nrow = p - d, ncol = d)
  B_Md <- rbind(I_d, C_Md)
  
  # Project covariates onto the subspace: U = X * B_Md
  U <- X %*% B_Md
  
  cv_error <- 0
  
  # Calculate the leave-one-out kernel estimator for each observation (Equation 2.11)
  for (i in 1:n) {
    u_i <- U[i, ]
    U_minus_i <- U[-i, , drop = FALSE]
    Y_minus_i <- Y[-i]
    
    # Calculate multivariate product kernel weights: \prod K_q((u_jk - u_ik)/h_k) / h_k
    weights <- rep(1, n - 1)
    for (k in 1:d) {
      u_diff <- (U_minus_i[, k] - u_i[k]) / h[k]
      weights <- weights * (K4(u_diff) / h[k])
    }
    
    sum_weights <- sum(weights)
    
    # Equation 2.11: \hat{g}_{C_Md}^{-i}(u)
    if (sum_weights == 0) {
      pred_i <- mean(Y_minus_i) # Fallback for isolated points
    } else {
      pred_i <- sum(weights * Y_minus_i) / sum_weights
    }
    
    # Accumulate squared error for Equation 2.12
    cv_error <- cv_error + (Y[i] - pred_i)^2
  }
  
  # Equation 2.12: CV_M(d_M, C_Md, h_d) = (1/n) * \sum (Y_i - pred_i)^2
  return(cv_error / n)
}

# =====================================================================
# 3. COMPUTATIONAL ALGORITHM
# Reference: Section 2.3 (Computational Algorithm)
# =====================================================================

estimate_CMS <- function(X, Y, max_iter = 1000) {
  n <- nrow(X)
  p <- ncol(X)
  
  # ----- Step 1 -----
  # Set CV(0) to the leave-one-out variance of Y (adapted from N_iy for CMS).
  cv_0 <- 0
  for (i in 1:n) {
    cv_0 <- cv_0 + (Y[i] - mean(Y[-i]))^2
  }
  CV_vals <- numeric(p)
  CV_vals[1] <- cv_0 / n # R is 1-indexed, so index 1 represents d=0
  
  best_models <- list()
  
  # Forward algorithm iterating through possible dimensions d = 1, ..., p-1
  for (d in 1:(p - 1)) {
    cat(sprintf("Testing structural dimension d = %d...\n", d))
    
    # ----- Step 2.1 -----
    # Initialize parameters. 
    # C_Md is initialized near 0.
    # Bandwidths are initialized proportionally to n^{-1/(2q+d)} where q=4.
    init_C <- rep(0.1, (p - d) * d)
    init_h <- rep(n^(-1 / (2*4 + d)), d) 
    init_params <- c(init_C, init_h)
    
    # ----- Step 2.2 -----
    # Perform nonlinear optimization. 
    # The paper notes the conjugate gradient algorithm of Fletcher and Reeves 
    # is a practical alternative to their custom Newton-CG line search to avoid
    # matrix operations. We use optim() with method = "CG".
    opt_res <- optim(
      par = init_params,
      fn = cv_objective_CMS,
      X = X, Y = Y, d = d, p = p, n = n,
      method = "CG",
      control = list(maxit = max_iter)
    )
    
    current_cv <- opt_res$value
    CV_vals[d + 1] <- current_cv
    
    # Reconstruct the optimal matrices to store them
    len_C <- (p - d) * d
    C_Md_opt <- matrix(opt_res$par[1:len_C], nrow = p - d, ncol = d)
    B_Md_opt <- rbind(diag(d), C_Md_opt)
    h_opt <- opt_res$par[(len_C + 1):length(opt_res$par)]
    
    best_models[[d]] <- list(
      d = d,
      C_M = C_Md_opt,
      B_M = B_Md_opt,
      h = h_opt,
      CV = current_cv
    )
    
    cat(sprintf("  CV(%d) = %f | CV(%d) = %f\n", d-1, CV_vals[d], d, current_cv))
    
    # ----- Step 3 -----
    # Implement Step 2 until CV(d) >= CV(d-1).
    if (current_cv >= CV_vals[d]) {
      cat(sprintf("Stopping rule triggered. Selected dimension: d_0 = %d\n", d - 1))
      
      if (d == 1) {
        return(list(d_0 = 0, message = "No predictive linear combinations found. d=0."))
      } else {
        return(best_models[[d - 1]])
      }
    }
  }
  
  # If the loop finishes without triggering the stopping rule, d = p-1 was optimal
  cat(sprintf("Reached maximum dimension. Selected dimension: d_0 = %d\n", p - 1))
  return(best_models[[p - 1]])
}

# --- Example Usage ---
set.seed(23)
n <- 500
p <- 4
X <- matrix(rnorm(n * p), n, p)
B_true <- matrix(c(1, 0, 0, 0,  0, 1, 1, 1), nrow = 4, ncol = 2) # True d=2
Y <- (X %*% B_true[,1])^2 + (X %*% B_true[,2]) + rnorm(n, 0, 0.5)
# 
result <- estimate_SDR_CMS(X, Y)
print(qr.Q(qr(result$B_M)))
      