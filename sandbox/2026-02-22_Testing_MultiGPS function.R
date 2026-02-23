library(here)

source(here("R", "multiGPS.R"))

# ---------------------------------------------------------
# 1. Simulate Test Data with Known Conditional Means
# ---------------------------------------------------------
set.seed(123)
n <- 1500

# Single confounder for easy verification
C_matrix <- cbind(Temp = rnorm(n, mean = 20, sd = 5))

# Exposures depend linearly on Temp
# If Temp = 20, E[PM2.5] = 10 + 0.5(20) = 20
# If Temp = 20, E[NO2]   = 20 - 0.2(20) = 16
PM2.5 <- 10 + 0.5 * C_matrix[, "Temp"] + rnorm(n, sd = 2)
NO2   <- 20 - 0.2 * C_matrix[, "Temp"] + rnorm(n, sd = 2)

X_matrix <- cbind(PM2.5 = PM2.5, NO2 = NO2)

# Ensure ols_fitter is in the environment
ols_fitter <- function(W, C) {
  train_data <- data.frame(W = W, C)
  lm(W ~ ., data = train_data)
}

# ---------------------------------------------------------
# 2. Define the Evaluation Grid and Fixed Confounders
# ---------------------------------------------------------
# We will evaluate a 20x20 grid around the expected conditional means
pm_seq <- seq(12, 28, length.out = 5)
no2_seq <- seq(8, 24, length.out = 5)
eval_grid <- expand.grid(PM2.5 = pm_seq, NO2 = no2_seq)

# We want the density evaluated specifically at Temp = 20
C_new <- data.frame(Temp = 20)

# Pre-calculate the rule-of-thumb bandwidth [cite: 427] to keep it constant across grid points
bw_fixed <- 1.25 * apply(X_matrix, 2, sd) * (n^(-0.2))

# ---------------------------------------------------------
# 3. Iterate over the Grid to Estimate Densities
# ---------------------------------------------------------
predicted_densities <- numeric(nrow(eval_grid))

# Define the libraries you want to ensemble
# e.g., standard GLM, Lasso, and Random Forest
my_libraries <- c("SL.glm")

# Set up custom cross-validation controls if needed
my_cv_control <- list(V = 5) # 5-fold CV for the SuperLearner meta-weights

# Note: This loop trains a new OLS model for every point in the grid
for (i in 1:nrow(eval_grid)) {
  target_t <- as.numeric(eval_grid[i, ])
  
  fit_gps <- multigps(
    x_eval    = target_t,
    X         = X_matrix,
    C         = C_matrix,
    h1        = bw_fixed,
    mu_fitter = function(W, C) SL_fitter(W, C, 
                                         SL.library = my_libraries, 
                                         cvControl = my_cv_control)
  )
  
  predicted_densities[i] <- predict(fit_gps, newdata = C_new)
}

# ---------------------------------------------------------
# 4. Plot the Density Surface and Verify Integration
# ---------------------------------------------------------
# Reshape the vector of densities back into a 20x20 matrix
density_matrix <- matrix(predicted_densities, nrow = length(pm_seq), ncol = length(no2_seq))

# Plot the contour map
contour(pm_seq, no2_seq, density_matrix,
        main = "Estimated Joint Density f(PM2.5, NO2 | Temp = 20)",
        xlab = "PM2.5", ylab = "NO2",
        nlevels = 10, col = terrain.colors(10))

# Add a point at the theoretical true mean (20, 16)
points(20, 16, col = "red", pch = 19, cex = 1.5)

# Calculate the numerical integral (Volume under the surface)
# Volume = Sum of (Density * Area of grid cell)
cell_area <- (pm_seq[2] - pm_seq[1]) * (no2_seq[2] - no2_seq[1])
integral_approx <- sum(density_matrix) * cell_area

cat("Approximate integral over the grid (should be ~1):", integral_approx, "\n")
