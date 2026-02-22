library(MASS) # For ginv
#library(SuperLearner)
library(Rcpp)
library(RcppArmadillo) # Ensure this is loaded

# Helper ------------------------------------------------------------------

# Computes the projection matrix Pi(beta)
Pi <- function(beta) {
  beta %*% solve(t(beta) %*% beta) %*% t(beta)
}

Delta <- function(beta1, beta2, type = "F"){
  Pi1 <- Pi(beta1)
  Pi2 <- Pi(beta2)
  
  return( norm(Pi1 - Pi2, type = type) )
}

# Nuisances ---------------------------------------------------------------

# estimate_m_gradient_cpp
# Careful! This function only takes a matrix object X (no data.frame), and it cannot have column names
cppFunction('
  Rcpp::List estimate_m_gradient_cpp(arma::mat X, arma::colvec Y, arma::mat beta, double b) {
    int n = X.n_rows;
    int d = beta.n_cols;
    arma::mat betaX = X * beta;
    
    // Calculate size of Z matrix
    int p_z = 1 + d + (d * (d + 1)) / 2;
    
    // Pre-allocate output structures
    arma::colvec m_est(n);
    arma::mat m_prime_est(n, d, arma::fill::zeros);
    arma::mat Z(n, p_z, arma::fill::ones);
    
    // Pre-allocate the ridge penalty
    arma::mat ridge = arma::eye(p_z, p_z) * 1e-8;
    
    // Pre-calculate constants for the Gaussian kernel
    // FIX 1: Explicitly cast the integer d to double to avoid ambiguous pow() overloads
   // The C++ Loop 
    for (int i = 0; i < n; i++) {
      arma::rowvec beta_x_i = betaX.row(i);
      
      arma::mat diff_betaX = betaX;
      diff_betaX.each_row() -= beta_x_i;
      
      arma::colvec dists = arma::sqrt(arma::sum(arma::square(diff_betaX), 1));
      arma::colvec u = dists / b;
      
      // NUMERICAL FIX: Drop the scaling constants. 
      // The weights are now strictly between 0 and 1, preventing underflow 
      // and stabilizing the ridge penalty!
      arma::colvec weights = arma::exp(-0.5 * arma::square(u));
      
      Z.cols(1, d) = diff_betaX;
      int col_idx = d + 1;
      for (int k = 0; k < d; k++) {
        for (int l = k; l < d; l++) {
          Z.col(col_idx) = diff_betaX.col(k) % diff_betaX.col(l);
          col_idx++;
        }
      }
      
      arma::mat ZW = Z;
      ZW.each_col() %= weights;
      
      arma::mat ZTWZ = Z.t() * ZW;
      arma::colvec ZTWY = ZW.t() * Y;
      
      arma::colvec coeffs = arma::solve(ZTWZ + ridge, ZTWY);
      
      m_est(i) = coeffs(0);
      
      for(int k = 0; k < d; k++) {
         m_prime_est(i, k) = coeffs(k + 1);
      }
    }
    
    return Rcpp::List::create(Rcpp::Named("m_est") = m_est,
                              Rcpp::Named("m_prime_est") = m_prime_est);
  }
', depends = "RcppArmadillo")

# 1. Fast Vectorized Variance Estimator (Eliminates the 500ms estimate_sigma2 bottleneck)
estimate_sigma2 <- function(X, residuals, h) {
  p <- ncol(X)
  
  # Calculate all pairwise distances at once in C
  dist_mat <- as.matrix(dist(X, method = "euclidean"))
  
  # Apply kernel to the entire matrix
  #weight_mat <- (1 / h^p) * dnorm(dist_mat / h)
  weight_mat <- dnorm(dist_mat / h)
  
  num <- weight_mat %*% (residuals^2)
  den <- rowSums(weight_mat) #+ 1e-10
  
  as.vector(pmax(num / den, 1e-6))
}

# # Replaces the 'h' bandwidth parameter with your machine learning library
# estimate_sigma2_SL <- function(X, residuals, SL.library = c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger")) {
#   
#   # SuperLearner rigidly requires the predictor matrix to be a data.frame
#   X_df <- as.data.frame(X)
#   
#   # The target outcome for a variance estimator is the squared residuals
#   Y_var <- residuals^2
#   
#   # Fit the ensemble model
#   # suppressWarnings keeps the console clean from algorithm-specific convergence chatter
#   suppressWarnings({
#     sl_fit <- SuperLearner(
#       Y = Y_var,
#       X = X_df,
#       family = gaussian(),
#       SL.library = SL.library,
#       cvControl = list(V = 5) # 5-fold CV internally for the ensemble weights
#     )
#   })
#   
#   # Extract the final predicted variances for every observation
#   sigma2_est <- sl_fit$SL.predict
#   
#   # CRITICAL SAFETY: Unconstrained ML models (like GLM or Random Forests) 
#   # can occasionally predict zero or negative values for specific outliers. 
#   # We rigorously floor the variance to prevent Newton-Raphson explosions.
#   as.vector(pmax(sigma2_est, 1e-6))
# }

# Estimate the ratio E(X / sigma2 | beta'X) / E(1 / sigma2 | beta'X)
estimate_Eq_betaX <- function(X, beta, sigma2, b){
  betaX <- X %*% beta
  
  # Calculate pairwise distances of projected X for all points at once
  dist_mat_betaX <- as.matrix(dist(betaX, method = "euclidean"))
  
  # Calculate the full n x n weight matrix simultaneously (Gaussian kernel)
  #W <- (1 / b^d) * dnorm(dist_mat_betaX / b)
  W <-  dnorm(dist_mat_betaX / b)
  
  # Scale columns by sigma2 using an optimized C-level sweep
  W_combined <- sweep(W, 2, sigma2, "/")
  
  # The entire ratio_term loop condenses into one BLAS-optimized matrix multiplication!
  ratio_num  <- W_combined %*% X
  ratio_den  <- rowSums(W_combined) #+ 1e-10
  ratio_term <- ratio_num / ratio_den
  
  return(ratio_term)
}

# Efficient score ---------------------------------------------------------
# 2. Fully Vectorized Main Update (Eliminates the 2400ms ratio_term loop bottlenecks)
compute_efficient_score_and_update <- function(X, Y, beta, b, h, SL = FALSE, sigma2 = NULL, alpha = 1) {
  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(beta)
  
  # Nuisance estimation
  m_results  <- estimate_m_gradient_cpp(X, matrix(Y, ncol=1), beta, b); residuals <- Y - m_results$m_est
  sigma2     <- if (!is.null(sigma2)) sigma2 else (if (SL) estimate_sigma2_SL(X, residuals) else estimate_sigma2(X, residuals, h))
  ratio_term <- estimate_Eq_betaX(X, beta, sigma2, b)
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
  vec_beta_next <- vec_beta_k + alpha * MASS::ginv(avg_outer_S + diag(1e-8, p*d)) %*% mean_vec_S
  #vec_beta_next <- vec_beta_k + MASS::ginv(avg_outer_S ) %*% mean_vec_S
  
  matrix(vec_beta_next, nrow = p, ncol = d)
}

# Main iterative loop (No QR decomposition inside the loop)
run_efficient_estimator <- function(X, Y, beta_init, b=NULL, h=NULL, max_iters = 100, SL = FALSE, sigma2 = NULL, alpha = 1, threshold = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(beta_init)
  
  beta_current <- beta_init
  if(is.null(threshold)){threshold <- p / n }
  
  # Bandwidths
  b <- if(is.null(b)) n^(-1 / (d + 4)) else b
  h <- if(is.null(h)) n^(-1 / (4 * p)) else b
  
  for (k in 1:max_iters) {
    cat(sprintf("Starting iteration %d...\n", k))
    
    beta_next <- compute_efficient_score_and_update(X, Y, beta_current, b, h, SL = SL, sigma2 = sigma2,alpha = alpha)
    #beta_next <- qr.Q(qr(beta_next))
    
    dist <- Delta(beta_current, beta_next)
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