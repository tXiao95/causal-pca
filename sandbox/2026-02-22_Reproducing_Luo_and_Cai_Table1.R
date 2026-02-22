### RUNTIME: 13min 20 seconds
## Slurm Job_id=35636637 Name=single_thread_job Ended, Run time 00:13:20, COMPLETED, ExitCode 0
library(here)

source(here("R", "sims_from_papers.R"))
source(here("R", "estimate_Seff.R"))

reproduce_table_1 <- function(n_reps = 1, n = 200, p = 10) {
  
  # Store the data generating functions in a list for easy iteration
  models <- list(
    sim_luo_cai_mod1,
    sim_luo_cai_mod2,
    sim_luo_cai_mod3,
    sim_luo_cai_mod4,
    sim_luo_cai_mod5
  )
  
  model_names <- c("I", "II", "III", "IV", "V")
  
  # 1. Pre-allocate a single dataframe for all aggregated results
  raw_results <- data.frame(
    Model = model_names,
    MAVE_mean   = numeric(5), MAVE_sd   = numeric(5),
    EE4_mean    = numeric(5), EE4_sd    = numeric(5),
    MAVE_mean_S = numeric(5), MAVE_sd_S = numeric(5),
    EE4_mean_S  = numeric(5), EE4_sd_S  = numeric(5)
  )
  
  cat("Starting simulation suite (", n_reps, " replications per model)...\n", sep = "")
  
  for (m in 1:5) {
    cat(sprintf("\nRunning Model %s...\n", model_names[m]))
    
    # 2. Use a matrix to store the 4 distance metrics for the current model
    # Columns: 1=MAVE(F), 2=EE4(F), 3=MAVE(S), 4=EE4(S)
    dist_mat <- matrix(NA_real_, nrow = n_reps, ncol = 4)
    
    for (i in seq_len(n_reps)) {
      # Use a dynamic seed to ensure distinct datasets per replicate 
      # but identical datasets if the whole suite is run again
      current_seed <- i * 100 + m
      sim_data <- models[[m]](n = n, p = p, seed = current_seed)
      
      X <- sim_data$X
      Y <- sim_data$Y
      beta_true <- sim_data$beta_true
      d <- sim_data$d
      
      # Fit MAVE initialization
      suppressWarnings({
        fit_mave <- MAVE::mave(Y ~ X, method = "meanMAVE")
      })
      beta_mave <- fit_mave$dir[[d]]
      
      # Calculate Theoretical Bandwidths for EE4
      b_val <- n^(-1 / (d + 4))
      h_val <- n^(-1 / (4 * p))
      
      # Fit EE4 
      beta_ee4 <- run_efficient_estimator(
        X = X, 
        Y = Y, 
        beta_init = beta_mave, 
        b = b_val, 
        h = h_val,
        max_iters = 100, 
        SL = FALSE, threshold = 0.01
      )
      
      # Calculate Subspace Distances
      dist_mat[i, 1] <- Delta(beta_true, beta_mave, "F")
      dist_mat[i, 2] <- Delta(beta_true, beta_ee4,  "F")
      dist_mat[i, 3] <- Delta(beta_true, beta_mave, "2")
      dist_mat[i, 4] <- Delta(beta_true, beta_ee4,  "2")
      
      # Print a tiny progress tracker every 10 iterations
      if (i %% 10 == 0) cat(i, " ")
    }
    
    # 3. Aggregate results seamlessly using matrix column operations
    means <- colMeans(dist_mat)
    sds   <- apply(dist_mat, 2, sd)
    
    # Store in the pre-allocated dataframe
    raw_results[m, c("MAVE_mean",   "EE4_mean",   "MAVE_mean_S", "EE4_mean_S")] <- means
    raw_results[m, c("MAVE_sd",     "EE4_sd",     "MAVE_sd_S",   "EE4_sd_S")]   <- sds
    
    cat("\n  MAVE Mean (F):", round(means[1], 3), " | EE4 Mean (F):", round(means[2], 3))
    cat("\n  MAVE Mean (S):", round(means[3], 3), " | EE4 Mean (S):", round(means[4], 3))
  }
  
  # Format the final table to match the paper's style (Mean over SD)
  formatted_table <- data.frame(
    Model  = model_names,
    MAVE   = sprintf("%.3f\n(%.3f)", raw_results$MAVE_mean,   raw_results$MAVE_sd),
    EE4    = sprintf("%.3f\n(%.3f)", raw_results$EE4_mean,    raw_results$EE4_sd),
    MAVE_S = sprintf("%.3f\n(%.3f)", raw_results$MAVE_mean_S, raw_results$MAVE_sd_S),
    EE4_S  = sprintf("%.3f\n(%.3f)", raw_results$EE4_mean_S,  raw_results$EE4_sd_S)
  )
  
  cat("\n\n--- Final Results (", n_reps, " Replications) ---\n", sep = "")
  print(formatted_table, row.names = FALSE)
  
  # Return the standardized output
  invisible(list(
    formatted = formatted_table,
    raw = raw_results
  ))
}

# Execute the simulation
# Note: This will take some time depending on your C++ and BLAS optimizations
results <- reproduce_table_1(n_reps = 5, n = 200, p = 10)

resultspath <- here("outputs", "experiments", "reproducing_luo_and_cai_Table1.rds")
saveRDS(results, resultspath)

tab <- readRDS(resultspath)


