library(testthat)
library(torch)
library(ggplot2)

# Ensure the environment is set up for testing
torch_set_num_threads(1)

source("R/nuisance_outcome_regression.R")
# =========================================================================
# 1. Setup Synthetic Non-Linear Data
# =========================================================================
set.seed(123)
n <- 100

# Create a 1D treatment (X) and 1D confounder (C)
X <- matrix(runif(n, -3, 3), ncol = 1)
C <- matrix(rnorm(n, 0, 1), ncol = 1)

# The Ground Truth: A highly non-linear sine wave + quadratic curve
# A linear model will fail completely here.
true_function <- function(x, c) { sin(3 * x) + (0.5 * x^2) + c }
Y <- true_function(X, C) #+ rnorm(n, 0, 0.1)

XC_df <- as.data.frame(cbind(X, C))
colnames(XC_df) <- c("X1", "C1")

# =========================================================================
# 2. Define a Simple Linear Benchmark
# =========================================================================
# This acts as our "dumb" linear model that cannot learn curves
glm_outcome_fitter <- function(Y, XC_df, ...) {
  fit <- lm(Y ~ ., data = XC_df)
  class(fit) <- c("glm_fit", class(fit))
  return(fit)
}
predict.glm_fit <- function(object, newdata, ...) {
  as.numeric(predict.lm(object, newdata))
}

# =========================================================================
# 3. Train Both Models
# =========================================================================
message("Training Linear Benchmark...")
linear_fit <- glm_outcome_fitter(Y, XC_df)

message("Training Neural Network...")
# Using 150 epochs to give it time to learn the sine wave
nn_fit <- nn_outcome_fitter(Y, XC_df, epochs = 1000, lr = 0.001, seed = 42)

# Generate Predictions
linear_preds <- predict(linear_fit, XC_df)
nn_preds <- predict(nn_fit, XC_df)

# =========================================================================
# 4. Formal Evaluation (The 'test_that' block)
# =========================================================================
test_that("NN captures non-linearity better than a linear baseline", {
  
  linear_mse <- mean((Y - linear_preds)^2)
  nn_mse <- mean((Y - nn_preds)^2)
  
  cat(sprintf("\nLinear Model MSE: %.4f\n", linear_mse))
  cat(sprintf("Neural Net MSE:   %.4f\n", nn_mse))
  
  # The NN MSE should be strictly and significantly lower
  expect_lt(nn_mse, linear_mse)
  
  # For this specific DGP, the NN should beat it by a massive margin 
  # (e.g., at least 50% better)
  expect_true(nn_mse < (0.5 * linear_mse))
})

# =========================================================================
# 5. Visual Evaluation
# =========================================================================
# To plot a clean 2D curve, we hold the confounder 'C' strictly at its mean (0)
plot_grid <- data.frame(
  X1 = seq(-3, 3, length.out = 300),
  C1 = 0 
)

# Predict across the grid
plot_grid$True_Y <- true_function(plot_grid$X1, plot_grid$C1)
plot_grid$Linear_Pred <- predict(linear_fit, plot_grid)

# FIX: Explicitly pass only the predictor columns to the neural network
plot_grid$NN_Pred <- predict(nn_fit, plot_grid[, c("X1", "C1")])

# Plotting the results
ggplot() +
  # Actual noisy data points (subsetting so the plot isn't too crowded)
  geom_point(data = data.frame(X1 = X, Y = Y)[1:1000,], 
             aes(x = X1, y = Y), color = "gray70", alpha = 0.5) +
  # The true underlying mathematical function
  geom_line(data = plot_grid, aes(x = X1, y = True_Y, color = "True Function"), 
            linewidth = 1.2, linetype = "dashed") +
  # The Linear Model's attempt (will just be a straight line)
  geom_line(data = plot_grid, aes(x = X1, y = Linear_Pred, color = "Linear Model"), 
            linewidth = 1) +
  # The Neural Network's attempt (should wrap tightly around the true function)
  geom_line(data = plot_grid, aes(x = X1, y = NN_Pred, color = "Neural Network"), 
            linewidth = 1.2) +
  scale_color_manual(values = c(
    "True Function" = "black", 
    "Linear Model" = "red", 
    "Neural Network" = "blue"
  )) +
  labs(
    title = "Neural Network vs. Linear Model on Non-Linear Data",
    subtitle = "Visualizing the Nuissance Fitter's ability to learn a sine wave",
    x = "Treatment (X)",
    y = "Outcome (Y)",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
