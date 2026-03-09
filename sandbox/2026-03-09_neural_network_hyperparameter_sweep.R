library(testthat)
library(torch)
library(ggplot2)

# Prevent SLURM thread thrashing
torch_set_num_threads(1)
#torch_set_num_interop_threads(1)

# =========================================================================
# 1. Setup Additive Simulation Data (Train and Test Sets)
# =========================================================================
# We generate an additive DGP: Non-linearities exist, but X and C do not interact.
generate_additive_data <- function(n, seed) {
  set.seed(seed)
  p <- 3 # 3 Exposures/Treatments
  q <- 2 # 2 Confounders
  
  X <- matrix(runif(n * p, -2, 2), n, p)
  C <- matrix(rnorm(n * q, 0, 1), n, q)
  
  # Additive DGP: X1 is quadratic, X2 is sinusoidal, X3 is linear, C are linear
  # No interaction terms (e.g., no X1 * C1)
  Y_true <- (0.5 * X[,1]^2) + sin(2 * X[,2]) + (0.75 * X[,3]) + C[,1] - (0.5 * C[,2])
  Y_obs <- Y_true + rnorm(n, 0, 0.5)
  
  XC_df <- as.data.frame(cbind(X, C))
  colnames(XC_df) <- c(paste0("X", 1:p), paste0("C", 1:q))
  
  list(XC_df = XC_df, Y = Y_obs)
}

# 1000 observations for training, 500 for strict out-of-sample validation
train_data <- generate_additive_data(n = 2500, seed = 123)
test_data  <- generate_additive_data(n = 500,  seed = 456)

# =========================================================================
# 2. Define the Neural Network C_fitter
# =========================================================================
C_fitter <- function(Y, XC_df, epochs = 150, lr = 0.005, seed = NULL, ...) {
  if (!is.null(seed)) torch_manual_seed(seed)
  
  X_mat <- as.matrix(XC_df)
  Y_mat <- matrix(Y, ncol = 1)
  
  x_means <- colMeans(X_mat)
  x_sds   <- apply(X_mat, 2, sd)
  x_sds[x_sds == 0] <- 1 
  X_scaled <- scale(X_mat, center = x_means, scale = x_sds)
  
  x_tensor <- torch_tensor(X_scaled, dtype = torch_float())
  y_tensor <- torch_tensor(Y_mat, dtype = torch_float())
  
  model <- nn_sequential(
    nn_linear(ncol(X_mat), 100),
    nn_silu(),
    nn_linear(100, 50),
    nn_silu(),
    nn_linear(50, 1)
  )
  
  optimizer <- optim_adam(model$parameters, lr = lr)
  scheduler <- lr_step(optimizer, step_size = 50, gamma = 0.5)
  criterion <- nn_mse_loss()
  
  model$train()
  for (epoch in 1:epochs) {
    optimizer$zero_grad()
    output <- model(x_tensor)
    loss <- criterion(output, y_tensor)
    loss$backward()
    optimizer$step()
    scheduler$step()
  }
  
  res <- list(model = model, x_means = x_means, x_sds = x_sds)
  class(res) <- "nn_fit"
  return(res)
}

predict.nn_fit <- function(object, newdata, ...) {
  object$model$eval()
  X_mat <- as.matrix(newdata)
  X_scaled <- scale(X_mat, center = object$x_means, scale = object$x_sds)
  x_tensor <- torch_tensor(X_scaled, dtype = torch_float())
  
  with_no_grad({
    preds <- object$model(x_tensor)
    out <- as.numeric(preds)
  })
  
  rm(x_tensor, preds)
  return(out)
}

# =========================================================================
# 3. Conduct the Learning Rate Sweep
# =========================================================================
learning_rates <- c(1, 0.5, 0.1,0.05, 0.01, 0.005, 0.001, 0.0005)
results_list <- list()

message("Sweeping learning rates for the C_fitter...")

for (lr in learning_rates) {
  # Train on the training set
  fit <- C_fitter(Y = train_data$Y, XC_df = train_data$XC_df, lr = lr, seed = 99)
  
  # Predict on the unseen validation set
  preds <- predict(fit, test_data$XC_df)
  mse <- mean((test_data$Y - preds)^2)
  
  results_list[[as.character(lr)]] <- data.frame(LearningRate = lr, OutOfSample_MSE = mse)
  
  # Clean up memory between fits
  rm(fit, preds)
  gc(verbose = FALSE)
}

results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

print(results_df[order(results_df$OutOfSample_MSE), ])

# Identify the absolute best rate
best_lr <- results_df$LearningRate[which.min(results_df$OutOfSample_MSE)]
best_mse <- min(results_df$OutOfSample_MSE)

# =========================================================================
# 4. Formal Evaluation & Baseline Comparison
# =========================================================================
# Fit a simple linear baseline to ensure the NN is actually learning the additive curves
lm_fit <- lm(Y ~ ., data = cbind(Y = train_data$Y, train_data$XC_df))
lm_preds <- predict(lm_fit, test_data$XC_df)
lm_mse <- mean((test_data$Y - lm_preds)^2)

test_that("Optimal NN Learning Rate out-performs linear baseline on Additive DGP", {
  
  cat(sprintf("\nBaseline GLM MSE:    %.4f\n", lm_mse))
  cat(sprintf("Best C_fitter MSE: %.4f (at LR = %s)\n", best_mse, best_lr))
  
  # The optimally tuned neural network should strictly beat the linear baseline
  expect_lt(best_mse, lm_mse)
})
