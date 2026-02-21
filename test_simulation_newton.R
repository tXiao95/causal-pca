# Data Generation: Model IV from Luo and Cai (2016)
generate_model_IV <- function(n = 200, p = 10) {
    # 1. Generate the covariance matrix Omega where entry (i,j) is 0.5^|i - j|
    Omega <- outer(1:p, 1:p, FUN = function(i, j) 0.5^abs(i - j))
    
    # 2. Generate the predictors X_N from N(0, Omega)
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Omega)
    
    # 3. Generate the independent error term epsilon_U from Uniform(-1, 1)
    epsilon <- runif(n, min = -1, max = 1)
    #epsilon <- rnorm(n, sd=0.6)
    # 4. Calculate the conditional mean function m_2(X_N)
    #m1 <- sin(0.3 * sum(colSums(X[,1:4])))
    m2 <-  10*X[, 1] / (0.5 + (1.5 + X[, 2])^2)
    
    # 5. Calculate the squared L2 norm for each row of X: ||X_N||_2^2
    X_norm_sq <- rowSums(X^2)
    
    # 6. Calculate the volatility multiplier based on the indicator functions
    # I(||X_N||_2^2 < p)/3 + 3I(||X_N||_2^2 >= p)
    volatility <- ifelse(X_norm_sq < p, 1/3, 3)
    
    # 7. Calculate the final Response Y
    Y <- m2 +  volatility*epsilon
    
    list(X = X, Y = Y)
}

set.seed(123)
B <- 20
distances <- vector(mode = "numeric", length = B)
distances2 <- vector(mode = "numeric", length = B)
for(i in 1:B){
  message(i)
# Generate the dataset
d <- 2  # dimension of central mean subspace
p <- 10 # number of predictors
n <- 200
#set.seed(1234)
sim_data_IV <- generate_model_IV(n = n, p = p)
X_test_IV <- sim_data_IV$X
Y_test_IV <- sim_data_IV$Y

# Set up an initial estimate (Step 1) to test the Newton-Raphson loop
# For Model IV, the central mean subspace is spanned by the 1st and 2nd predictors.
beta_init_IV <- MAVE::mave(Y_test_IV ~ X_test_IV, method = "meanMAVE")$dir[[2]]

# Calculate recommended bandwidths for testing
b_val <- n^(-1 / (d + 4))
h_val <- n^(-1 / (4 * p)) # Corresponds to the "EE 4" estimator

# Run the estimator (assuming run_efficient_estimator is loaded in your environment)
final_beta_IV <- run_efficient_estimator(X_test_IV, Y_test_IV, beta_init_IV, b_val, h_val)

beta0 <- matrix(0, nrow=p,ncol=d); beta0[1,1] <- 1; beta0[2,2] <- 1

distances[i]  <- frob_norm(projection_matrix(beta0), projection_matrix(beta_init_IV))
distances2[i] <- frob_norm(projection_matrix(beta0), projection_matrix(final_beta_IV))

}

mean(distances); sd(distances)
mean(distances2); sd(distances2)