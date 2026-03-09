# =========================================================================
# 1. Setup and Thread Management
# =========================================================================
library(torch)

# CRITICAL FOR SLURM: Force torch to use exactly 1 thread. 
# This prevents thread thrashing and CPU lockups.
torch_set_num_threads(1)
#torch_set_num_interop_threads(1)

# =========================================================================
# 2. Define the Neural Network Fitter and Predictor
# =========================================================================
nn_outcome_fitter <- function(Y, XC_df, epochs = 150, lr = 0.005, ...) {
  X_mat <- as.matrix(XC_df)
  Y_mat <- matrix(Y, ncol = 1)
  
  # Scaling parameters
  x_means <- colMeans(X_mat)
  x_sds   <- apply(X_mat, 2, sd)
  x_sds[x_sds == 0] <- 1 
  X_scaled <- scale(X_mat, center = x_means, scale = x_sds)
  
  x_tensor <- torch_tensor(X_scaled, dtype = torch_float())
  y_tensor <- torch_tensor(Y_mat, dtype = torch_float())
  
  # Architecture: 100x50 with SiLU
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
  
  # Training Loop
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
  })
  
  return(as.numeric(preds))
}

# =========================================================================
# 3. Your Specific `outcome_model` Wrapper
# =========================================================================
outcome_model <- function(Y, X, C, mu_fitter, ...) {
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL
  
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  p <- ncol(X_df)
  q <- ncol(C_df)
  
  if (is.null(orig_X_names)) colnames(X_df) <- paste0("X", 1:p) else colnames(X_df) <- make.names(orig_X_names, unique = TRUE)
  if (is.null(orig_C_names)) colnames(C_df) <- paste0("C", 1:q) else colnames(C_df) <- make.names(orig_C_names, unique = TRUE)
  
  df_train <- cbind(X_df, C_df)
  inner_fit <- mu_fitter(Y, df_train, ...)
  
  res <- list(inner_fit = inner_fit, X_names = colnames(X_df), C_names = colnames(C_df), p = p, q = q)
  class(res) <- "outcome_model"
  return(res)
}

predict.outcome_model <- function(object, newdata, ...) {
  newdata <- as.data.frame(newdata)
  req_cols <- c(object$X_names, object$C_names)
  newdata <- newdata[, req_cols, drop = FALSE]
  
  preds <- predict(object$inner_fit, newdata = newdata, ...)
  if (is.list(preds) && "pred" %in% names(preds)) return(as.numeric(preds$pred))
  return(as.numeric(preds))
}

# =========================================================================
# 4. Execute the Local Test
# =========================================================================
message("Generating synthetic causal data...")
set.seed(999)
n_train <- 500
p_dim <- 3
q_dim <- 2

# Dummy Data
X_train <- matrix(rnorm(n_train * p_dim), n_train, p_dim)
C_train <- matrix(rnorm(n_train * q_dim), n_train, q_dim)
# Non-linear outcome to give the neural net something to learn
Y_train <- 2 * sin(X_train[,1]) + X_train[,2]^2 + C_train[,1] + rnorm(n_train, 0, 0.5)

# Wrap the fitter so it matches the expected signature
my_nn_wrapper <- function(Y, XC_df, ...) nn_outcome_fitter(Y, XC_df, epochs = 150, lr = 0.01, ...)

message("Training 100x50 SiLU Neural Network...")
start_time <- Sys.time()
fit_obj <- outcome_model(Y = Y_train, X = X_train, C = C_train, mu_fitter = my_nn_wrapper)
end_time <- Sys.time()
print(end_time - start_time)

message("Predicting on target evaluation points...")
# Create a dummy evaluation grid (first 5 rows)
eval_grid <- cbind(X_train[1:5, ], C_train[1:5, ])
colnames(eval_grid) <- c(paste0("X", 1:p_dim), paste0("C", 1:q_dim))

predictions <- predict(fit_obj, newdata = eval_grid)
print(predictions)

message("Executing Garbage Collection cleanup...")
# Prove memory clears safely
rm(fit_obj, predictions, x_tensor, y_tensor) # x_tensor/y_tensor are inside the function scope, but good practice
gc()
message("Test Complete. Ready for the cluster!")