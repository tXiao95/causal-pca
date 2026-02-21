library(Rcpp)
library(RcppArmadillo) # Ensure this is loaded

# Use the depends argument directly in R
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
    double h_pow = std::pow(b, (double)d); 
    double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * arma::datum::pi);
    
    for (int i = 0; i < n; i++) {
      arma::rowvec beta_x_i = betaX.row(i);
      
      arma::mat diff_betaX = betaX;
      diff_betaX.each_row() -= beta_x_i;
      
      arma::colvec dists = arma::sqrt(arma::sum(arma::square(diff_betaX), 1));
      arma::colvec u = dists / b;
      arma::colvec weights = (1.0 / h_pow) * inv_sqrt_2pi * arma::exp(-0.5 * arma::square(u));
      
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
      
      // FIX 2: Manually map the coefficients to the matrix row to avoid transpose template errors
      for(int k = 0; k < d; k++) {
         m_prime_est(i, k) = coeffs(k + 1);
      }
    }
    
    return Rcpp::List::create(Rcpp::Named("m_est") = m_est,
                              Rcpp::Named("m_prime_est") = m_prime_est);
  }
', depends = "RcppArmadillo")