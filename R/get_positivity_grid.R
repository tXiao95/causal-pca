library(geometry)
library(ks)
library(MASS)

get_convex_positivity_grid <- function(data, percentile = 0, n_points = 10) {
  # Ensure data is a matrix for mathematical operations
  data <- as.matrix(data)
  dims <- ncol(data)
  
  # 1. Determine the bounds for the rectangular grid
  ranges <- list()
  for (i in 1:dims) {
    # Calculate percentiles (R's quantile uses probs between 0 and 1)
    lower <- quantile(data[, i], probs = percentile / 100, names = FALSE)
    upper <- quantile(data[, i], probs = (100 - percentile) / 100, names = FALSE)
    
    # Create n_points linearly spaced between the bounds
    ranges[[i]] <- seq(lower, upper, length.out = n_points)
  }
  
  # Set list names to preserve column names if they exist
  if (!is.null(colnames(data))) {
    names(ranges) <- colnames(data)
  }
  
  # 2. Generate the Cartesian Product (Rectangular Grid)
  # expand.grid generates all combinations; convert to matrix for linear algebra
  grid_points <- as.matrix(expand.grid(ranges))
  
  # 3. Compute the Convex Hull of the source data
  # The "n" option tells Qhull to return the normals and offsets of the facets
  hull <- geometry::convhulln(data, options = "n")
  
  # Extract A (normals) and b (offsets)
  # hull$normals is an M x (D+1) matrix, where M is the number of facets
  A <- hull$normals[, 1:dims, drop = FALSE]
  b <- hull$normals[, dims + 1]
  
  # 4. Filter points: Check which are inside the Convex Hull
  # Equation form: Ax + b <= 0
  # Calculate A * x.T for all points (grid_points %*% t(A))
  prod_matrix <- grid_points %*% t(A)
  
  # Add the offset 'b' to each corresponding column
  # sweep() applies the operation "+" across margin 2 (columns)
  eq_vals <- sweep(prod_matrix, 2, b, "+")
  
  # A point is inside if all equations evaluate to <= 0 (allowing for tiny float error)
  epsilon <- 1e-9
  is_inside <- apply(eq_vals, 1, function(row) all(row <= epsilon))
  
  # Return only the points satisfying the hull constraints
  return(grid_points[is_inside, , drop = FALSE])
}

get_kde_positivity_grid <- function(data, percentile = 0, n_points = 20, threshold = 0.05, bw_selection = c("silverman", "plugin")) {
  data <- as.matrix(data)
  
  # 1. Create the base rectangular grid using percentiles
  ranges <- lapply(1:ncol(data), function(i) {
    lower <- quantile(data[, i], probs = percentile / 100)
    upper <- quantile(data[, i], probs = (100 - percentile) / 100)
    seq(lower, upper, length.out = n_points)
  })
  
  # Create Cartesian product
  grid_points <- as.matrix(expand.grid(ranges))
  colnames(grid_points) <- colnames(data)
  
  # 2. Perform Kernel Density Estimation (KDE)
  # Hpi is a plug-in bandwidth selector—it determines how "tight" the fit is.
  #H_mat <- ks::Hscv(x = data)
  # Replace ks::Hpi(x = data) with:
  if(bw_selection == "silverman"){
    H_mat <- diag(apply(data, 2, sd) * (4/(ncol(data)+2)/nrow(data))^(1/(ncol(data)+4)))
  } else if(bw_selection == "plugin"){
    H_mat <- ks::Hpi(x = data)
  }else{
    stop("Please select a valid 'bw_selection' argument from 'silverman' or 'plugin'")
  }
  fhat <- ks::kde(x = data, H = H_mat)
  
  # 3. Calculate density for our grid points vs. actual data
  grid_density <- predict(fhat, x = grid_points)
  actual_density <- predict(fhat, x = data)
  
  # 4. Set a threshold
  # Any grid point with density lower than the bottom 5% of our 
  # observed data density is discarded.
  threshold <- quantile(actual_density, threshold)
  
  # Filter the grid
  keep_idx <- which(grid_density >= threshold)
  return(grid_points[keep_idx, , drop = FALSE])
}
