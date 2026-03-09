library(SuperLearner)
library(torch)

#' Fit a Global Outcome Regression Object E[Y | X, C]
#' 
#' @param Y Numeric vector of outcomes.
#' @param X Numeric matrix or data frame of observed treatments.
#' @param C Numeric matrix or data frame of observed confounders.
#' @param mu_fitter Function(Y, XC_df) that trains and returns a model.
#' @return An S3 object of class "outcome_model".

outcome_model <- function(Y, X, C, mu_fitter, ...) {
  # 1. Capture original names BEFORE any coercion or processing
  orig_X_names <- if(is.matrix(X) || is.data.frame(X)) colnames(X) else NULL
  orig_C_names <- if(is.matrix(C) || is.data.frame(C)) colnames(C) else NULL
  
  # 2. Convert to data frames
  X_df <- as.data.frame(X)
  C_df <- as.data.frame(C)
  p <- ncol(X_df)
  q <- ncol(C_df)
  
  # 3. Apply the Contract: Use original names if they exist, otherwise generate ours
  # We use 'tx_' and 'conf_' prefixes to guarantee zero collision if we generate them
  if (is.null(orig_X_names)) {
    colnames(X_df) <- paste0("X", 1:p)
  } else {
    colnames(X_df) <- make.names(orig_X_names, unique = TRUE)
  }
  
  if (is.null(orig_C_names)) {
    colnames(C_df) <- paste0("C", 1:q)
  } else {
    colnames(C_df) <- make.names(orig_C_names, unique = TRUE)
  }
  
  # 4. Final safety check for overlaps
  overlapping <- intersect(colnames(X_df), colnames(C_df))
  if (length(overlapping) > 0) {
    stop("Overlapping column names detected: ", paste(overlapping, collapse = ", "))
  }
  
  # Fit using the consistent names
  df_train <- cbind(X_df, C_df)
  inner_fit <- mu_fitter(Y, df_train, ...)
  
  res <- list(
    inner_fit = inner_fit,
    X_names = colnames(X_df),
    C_names = colnames(C_df),
    p = p,
    q = q
  )
  class(res) <- "outcome_model"
  return(res)
}

#' Predict Method for Outcome Model
predict.outcome_model <- function(object, newdata, ...) {
  
  # Ensure newdata is a data frame
  newdata <- as.data.frame(newdata)
  
  # Combine required column names
  req_cols <- c(object$X_names, object$C_names)
  
  # 1. Check for missing columns
  missing_cols <- setdiff(req_cols, colnames(newdata))
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from 'newdata': ", 
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # 2. Subset and FORCE column order to match the training data exactly
  newdata <- newdata[, req_cols, drop = FALSE]
  
  # Predict using the inner model
  preds <- predict(object$inner_fit, newdata = newdata, ...)
  
  # Extract numeric vector (SuperLearner returns a list with a $pred matrix)
  if (is.list(preds) && "pred" %in% names(preds)) {
    return(as.numeric(preds$pred))
  }
  
  return(as.numeric(preds))
}

# Define the wrapper
SL_outcome_fitter <- function(Y, XC_df, SL.lib = c("SL.glm", "SL.glmnet", "SL.xgboost", "SL.earth"), ...) {
  SuperLearner::SuperLearner(Y = Y, X = XC_df, family = gaussian(), SL.lib = SL.lib, ...)
}

# Neural Network ----------------------------------------------------------

torch_set_num_threads(1)

nn_outcome_fitter <- function(Y, XC_df, epochs = 150, lr = 0.005, seed = NULL, ...) {
  
  # ENFORCE REPRODUCIBILITY: Use the simulation iteration seed if provided
  if (!is.null(seed)) {
    torch_manual_seed(seed)
  }
  
  X_mat <- as.matrix(XC_df)
  Y_mat <- matrix(Y, ncol = 1)
  
  # Scaling parameters
  x_means <- colMeans(X_mat)
  x_sds   <- apply(X_mat, 2, sd)
  x_sds[x_sds == 0] <- 1 # Safety catch
  
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
    out <- as.numeric(preds)
  })
  
  # EXPLICIT CLEANUP: Destroy the C++ tensors to prevent memory leaks
  rm(x_tensor, preds)
  
  return(out)
}

NN_outcome_fitter <- function(Y, XC_df, ...) {
  # Passes any ... arguments (like seed) down to nn_outcome_fitter
  nn_outcome_fitter(Y, XC_df, epochs = 150, lr = 0.005, ...)
}
