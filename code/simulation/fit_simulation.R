library(data.table)
library(here)
library(parallel)
library(SuperLearner)
library(MAVE)

source(here("R/pCCA.R"))
source(here("R/csMAVE.R"))
source(here("R/Seff.R"))
source(here("R/crossfit_ERS.R"))
source(here("R/compute_new_response_and_exposure.R"))

source(here("R/estimate_ERS.R"))
source(here("R/estimate_pseudo_outcomes.R"))
source(here("R/estimate_residualized_pair.R"))

source(here("R/nuisance_outcome_regression.R"))
source(here("R/nuisance_outcome_C.R"))
source(here("R/nuisance_gps.R"))

source(here("R/misc.R"))
source(here("R/simulate_old_papers.R"))
source(here("R/simulate_data.R"))
source(here("R/simulate_pfas_birthweight.R"))

# -------------------------------------------------------------------------
# Well-Specified Fitters for simulate_causal_sdr DGP
# -------------------------------------------------------------------------

# 1. Outcome Regression: E[Y | X, C]
SL_fitter_simple <- function(Y, XC_df, ...) {
  SL_outcome_fitter(Y, XC_df, SL.lib = c("SL.glm", "SL.glmnet", "SL.earth"), ...)
}

SL_fitter <- function(Y, XC_df, ...) {
  SL_outcome_fitter(Y, XC_df, SL.lib = c("SL.glm", "SL.glmnet", 
                                         "SL.earth", "SL.gam",
                                         "SL.ranger", "SL.xgboost"), ...)
}

# DR_outcome_fitter <- function(Y, XC_df, ...) {
#    nn_outcome_fitter(Y, XC_df, arch = "100_50", epochs = 150, lr = 0.05, ...)
# }

DR_outcome_fitter <- function(Y, XC_df, ...) {
   nn_outcome_fitter(Y, XC_df, arch = "128_64_32", epochs = 1000, lr = 0.01, ...)
}

# 2. GPS Model: f(X | C)
well_specified_gps_fitter <- function(X, C, ...) {
  mvn_fitter(X, C, method_gps = "linear", ...)
}

# -------------------------------------------------------------------------
# Evaluation Module for Parallel Processing
# -------------------------------------------------------------------------

evaluate_method_group <- function(group, sim, n, task_id, x_eval_fixed, mu_true_fixed) {
  Y <- sim$Y; X <- sim$X; C <- sim$C
  beta0 <- sim$beta_true
  d0 <- sim$d
  p_dims <- ncol(X)
  n_eval <- nrow(x_eval_fixed)
  
  # ========================================================================
  # HELPER FUNCTIONS
  # ========================================================================
  
  # 1. Format Error DT
  format_err <- function(method, frob_d0, frob_dhat, time, dhat_val) {
    data.table(ID = task_id, n = n, method = method, time = time, dhat = dhat_val,
               error_type = c("d0", "dhat"), frob_norm = c(frob_d0, frob_dhat))
  }
  
  # 2. Format Diagonal DT (Long Format)
  format_diag <- function(method, b_d0) {
    if (!is.matrix(b_d0) || any(is.na(b_d0))) {
      diag_vals <- rep(NA_real_, p_dims)
    } else {
      Q <- tryCatch(qr.Q(qr(b_d0)), error = function(e) b_d0)
      diag_vals <- diag(tcrossprod(Q))
    }
    data.table(ID = task_id, n = n, method = method, dimension = 1:p_dims, P_diag_val = diag_vals)
  }
  
  # 3. Evaluate ERS DT (Wide Format for Z, est, CIs)
  evaluate_ers_dt <- function(method, b_hat) {
    dt <- data.table(ID = task_id, n = n, method = method, eval_id = 1:n_eval)
    for(j in 1:p_dims) dt[[paste0("X_", j)]] <- x_eval_fixed[, j]
    dt$mu_true <- mu_true_fixed
    for(j in 1:p_dims) dt[[paste0("Z_", j)]] <- NA_real_
    dt[, c("est", "ci_lower", "ci_upper")] <- NA_real_
    
    if (!is.matrix(b_hat) || any(is.na(b_hat))) return(dt)
    
    Z_train <- X %*% b_hat
    z_eval  <- x_eval_fixed %*% b_hat
    
    for(j in 1:ncol(b_hat)) dt[[paste0("Z_", j)]] <- z_eval[, j]
    
    res <- tryCatch({
      suppressWarnings({
        crossfit_ERS(Y = Y, X = Z_train, C = C, x_eval = z_eval, 
                     estimator = "DR", L = 5, optimize_bw = TRUE,
                     outcome_fitter = DR_outcome_fitter,
                     #outcome_fitter = SL_fitter_simple,
                     #outcome_fitter = SL_fitter,
                     gps_fitter = well_specified_gps_fitter)
      })
    }, error = function(e) NULL)
    
    if(!is.null(res) && !is.null(res$results)) {
      dt$est <- res$results$estimate
      dt$ci_lower <- res$results$ci_lower
      dt$ci_upper <- res$results$ci_upper
    }
    return(dt)
  }
  
  # 4. Unified Recorder (Bundles all 3 steps seamlessly)
  record_method <- function(name, b_d0, b_dhat, time, dhat_val, do_ers = TRUE) {
    err_dt  <- format_err(name, Delta(beta0, b_d0, "F"), 
                          if(is.matrix(b_dhat)) Delta(beta0, b_dhat, "F") else NA, 
                          time, dhat_val)
    diag_dt <- format_diag(name, b_d0)
    ers_dt  <- if(do_ers) evaluate_ers_dt(name, b_dhat) else data.table()
    
    list(err = err_dt, diag = diag_dt, ers = ers_dt)
  }
  
  # ========================================================================
  # METHOD EVALUATION
  # ========================================================================
  res_lists <- list()
  
  if (group == "Base") {
    message("Running baseline methods")
    
    t_pca <- system.time({ pca <- prcomp(X, center = TRUE, scale. = TRUE) })
    res_lists[[length(res_lists) + 1]] <- record_method("PCA", pca$rotation[, 1:d0, drop = FALSE], NA, t_pca["elapsed"], NA, do_ers = FALSE)
    
    t_pcca <- system.time({ b_pcca <- pCCA(Y, X, C)[, 1:d0, drop = FALSE] })
    res_lists[[length(res_lists) + 1]] <- record_method("pCCA", b_pcca, b_pcca, t_pcca["elapsed"], d0, do_ers = TRUE)
    
    t_mave <- system.time({ reg_MAVE <- MAVE::mave(Y ~ X, method = "meanMAVE") })
    dhat_m <- MAVE::mave.dim(reg_MAVE)$dim.min
    b_mave_d0 <- reg_MAVE$dir[[d0]]
    b_mave_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) reg_MAVE$dir[[dhat_m]] else NA
    res_lists[[length(res_lists) + 1]] <- record_method("MAVE", b_mave_d0, b_mave_dhat, t_mave["elapsed"], dhat_m, do_ers = TRUE)
    
    t_ee <- system.time({ b_ee_d0 <- run_efficient_estimator(X, Y, beta_init = b_mave_d0) })
    b_ee_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) run_efficient_estimator(X, Y, beta_init = b_mave_dhat) else NA
    res_lists[[length(res_lists) + 1]] <- record_method("EE", b_ee_d0, b_ee_dhat, t_ee["elapsed"], NA, do_ers = FALSE)
    
    # Full X Benchmark (Only ERS)
    res_lists[[length(res_lists) + 1]] <- list(err = data.table(), diag = data.table(), ers = evaluate_ers_dt("Full_X", diag(p_dims)))
    
  } else if (group == "RA_DR_PO") {
    message("Running Causal Pipeline (RA, DR, PO)...")
    causal_methods <- c("RA", "DR", "PO")
    
    t_cre <- system.time({
      cre_obj <- compute_new_response_and_exposure(Y = Y, X = X, C = C, method = causal_methods, L = 5,
                                                   #outcome_fitter = SL_fitter_simple,
                                                   #outcome_fitter = SL_fitter,
                                                   outcome_fitter = DR_outcome_fitter,
                                                   gps_fitter = well_specified_gps_fitter)
    })
    amortized_time <- t_cre["elapsed"] / length(causal_methods)
    
    for (m in causal_methods) {
      new_X <- cre_obj$new_X[[m]]; new_Y <- cre_obj$new_Y[[m]]
      df <- data.frame(newY = new_Y, new_X)
      
      t_mave <- system.time({ reg_MAVE <- MAVE::mave(newY ~ ., data = df, method = "meanMAVE") })
      dhat_m <- MAVE::mave.dim(reg_MAVE)$dim.min
      b_mave_d0 <- reg_MAVE$dir[[d0]]
      b_mave_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) reg_MAVE$dir[[dhat_m]] else NA
      
      res_lists[[length(res_lists) + 1]] <- record_method(paste0(m, "-MAVE"), b_mave_d0, b_mave_dhat, 
                                                          amortized_time + t_mave["elapsed"], dhat_m, do_ers = TRUE)
      
      t_ee <- system.time({ b_ee_d0 <- run_efficient_estimator(new_X, new_Y, beta_init = b_mave_d0) })
      b_ee_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) run_efficient_estimator(new_X, new_Y, beta_init = b_mave_dhat) else NA
      res_lists[[length(res_lists) + 1]] <- record_method(paste0(m, "-EE"), b_ee_d0, b_ee_dhat, t_ee["elapsed"], NA, do_ers = FALSE)
    }
    
  } else if (group == "RP") {
    message("Running Residualized Pair (RP)...")
    t_cre <- system.time({
      cre_obj <- compute_new_response_and_exposure(Y = Y, X = X, C = C, method = "RP", L = 5,
                                                   #C_fitter = well_specified_C_fitter)
                                                   #C_fitter = DR_outcome_fitter)
                                                   C_fitter = SL_fitter_simple)
    })
    
    new_X <- cre_obj$new_X; new_Y <- cre_obj$new_Y
    df <- data.frame(newY = new_Y, new_X)
    
    t_mave <- system.time({ reg_MAVE <- MAVE::mave(newY ~ ., data = df, method = "meanMAVE") })
    dhat_m <- MAVE::mave.dim(reg_MAVE)$dim.min
    b_mave_d0 <- reg_MAVE$dir[[d0]]
    b_mave_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) reg_MAVE$dir[[dhat_m]] else NA
    
    res_lists[[length(res_lists) + 1]] <- record_method("RP-MAVE", b_mave_d0, b_mave_dhat, t_cre["elapsed"] + t_mave["elapsed"], dhat_m, do_ers = TRUE)
    
    t_ee <- system.time({ b_ee_d0 <- run_efficient_estimator(new_X, new_Y, beta_init = b_mave_d0) })
    b_ee_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) run_efficient_estimator(new_X, new_Y, beta_init = b_mave_dhat) else NA
    res_lists[[length(res_lists) + 1]] <- record_method("RP-EE", b_ee_d0, b_ee_dhat, t_ee["elapsed"], NA, do_ers = FALSE)
    
  } else if (group == "Oracle") {
    message("Running Oracle")
    t_mave <- system.time({ obj <- MAVE::mave(sim$mu_X ~ X, method = "meanMAVE") })
    dhat_m <- MAVE::mave.dim(obj)$dim.min
    b_mave_d0 <- obj$dir[[d0]]
    b_mave_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) obj$dir[[dhat_m]] else NA
    
    res_lists[[length(res_lists) + 1]] <- record_method("Oracle-MAVE", b_mave_d0, b_mave_dhat, t_mave["elapsed"], dhat_m, do_ers = TRUE)
    
    t_ee <- system.time({ b_ee_d0 <- run_efficient_estimator(X, Y, beta_init = b_mave_d0) })
    b_ee_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) run_efficient_estimator(X, Y, beta_init = b_mave_dhat) else NA
    res_lists[[length(res_lists) + 1]] <- record_method("Oracle-EE", b_ee_d0, b_ee_dhat, t_ee["elapsed"], NA, do_ers = FALSE)
  }
  
  # ========================================================================
  # COMBINE AND RETURN
  # ========================================================================
  return(list(
    err  = rbindlist(lapply(res_lists, `[[`, "err"), use.names = TRUE, fill = TRUE),
    diag = rbindlist(lapply(res_lists, `[[`, "diag"), use.names = TRUE, fill = TRUE),
    ers  = rbindlist(lapply(res_lists, `[[`, "ers"), use.names = TRUE, fill = TRUE)
  ))
}

# -------------------------------------------------------------------------
# Main Execution Loop
# -------------------------------------------------------------------------

main <- function() {
  message("TASK_ID: ", TASK_ID)
  message("Number of cores: ", N_CORES)
  message("Experiment: ", EXPERIMENT)
  
  N_vector <- c(100, 500, 1000, 2500, 5000)
  #N_vector <- c(303, 500, 1000)
  #groups   <- c("Base", "RA_DR_PO", "RP", "Oracle")
  groups   <- c("Base", "RA_DR_PO", "Oracle")
  
  # -------------------------------------------------------------------
  # Generate a Fixed, 100-Point Evaluation Grid
  # -------------------------------------------------------------------
  message("Generating fixed evaluation grid for ERS...")
  
  # HARDCODED SEED: Ensures the grid is mathematically identical 
  # across ALL tasks, ALL sample sizes, and ALL methods.
  set.seed(99999)
  
  is_weak   <- (EXPERIMENT == "weak_dim")
  int_coef  <- if(EXPERIMENT == "additive") 0.0 else 5.0  #(MAIN EXPERIMENT)
  is_sparse <- if(NAME == "main_paper_final_results_nnet") TRUE else FALSE
  
  #int_coef <- if(EXPERIMENT == "additive") 0.0 else 10.0   #(PFAS EXPERIMENT)
  #int_coef <- if(EXPERIMENT == "additive") 0.0 else 200.0   #(PFAS EXPERIMENT)
  
  # MAIN EXPERIMENT
  sim_grid <- simulate_causal_sdr_simple(n = 100, p = 10, q = 5, noise_sd = 0.5, 
                                         rho_X = 0.8, 
                                         interaction_coef = int_coef, 
                                         weak_dim_signal = is_weak,
                                         sparse = is_sparse)
  
  #sim_grid <- simulate_pfas_birthweight(n = 100, interaction_coef = int_coef)
  
  x_eval_fixed  <- sim_grid$X
  mu_true_fixed <- sim_grid$mu_X
  
  message("Running with ", N_CORES, " cores.")
  
  # Run simulations across N_vector
  tables_list <- lapply(N_vector, function(n) {
    message("\nSimulation experiment with sample size ", n)
    
    # Reset the seed using the TASK_ID so the training data varies properly
    # across the SLURM array, while retaining reproducibility.
    set.seed(TASK_ID)
    
    # MAIN EXPERIMENT
    sim <- simulate_causal_sdr_simple(n = n, p = 10, q = 5, noise_sd = 0.5, rho_X = 0.8,
                                      interaction_coef = int_coef, 
                                      weak_dim_signal = is_weak,
                                      sparse = is_sparse) 
    
    # PFAS EXPERIMENT
    #sim <- simulate_pfas_birthweight(n = n, interaction_coef = int_coef)
    
    # Run the families in parallel
    res_list <- mclapply(groups, function(g) {
      evaluate_method_group(group = g, sim = sim, n = n, task_id = TASK_ID,
                            x_eval_fixed = x_eval_fixed, mu_true_fixed = mu_true_fixed)
    }, mc.cores = N_CORES)
    
    # Combine outputs for this sample size N
    list(
      err  = rbindlist(lapply(res_list, `[[`, "err"), use.names = TRUE, fill = TRUE),
      diag = rbindlist(lapply(res_list, `[[`, "diag"), use.names = TRUE, fill = TRUE),
      ers  = rbindlist(lapply(res_list, `[[`, "ers"), use.names = TRUE, fill = TRUE)
    )
  })
  
  # Aggregate all sample sizes
  final_err  <- rbindlist(lapply(tables_list, `[[`, "err"), use.names = TRUE, fill = TRUE)
  final_diag <- rbindlist(lapply(tables_list, `[[`, "diag"), use.names = TRUE, fill = TRUE)
  final_ers  <- rbindlist(lapply(tables_list, `[[`, "ers"), use.names = TRUE, fill = TRUE)
  
  
  out_dir <- here("outputs/simulation/jasa-initial-submission", NAME, EXPERIMENT)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  base_filename <- paste0("sim-", sprintf("%03d", TASK_ID))
  
  data.table::fwrite(final_err,  file.path(out_dir, paste0(base_filename, "-error.csv")))
  data.table::fwrite(final_diag, file.path(out_dir, paste0(base_filename, "-diag.csv")))
  data.table::fwrite(final_ers,  file.path(out_dir, paste0(base_filename, "-ers.csv")))
  
  message("Results saved to ", out_dir)
}

# Arguments for main ------------------------------------------------------
if(interactive()){
  TASK_ID    <- 2
  N_CORES    <- 1
  NAME       <- "pfas"
  EXPERIMENT <- "interaction"
} else{
  TASK_ID    <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  N_CORES    <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
  NAME       <- as.character( commandArgs(trailingOnly=TRUE)[1] )
  EXPERIMENT <- as.character( commandArgs(trailingOnly=TRUE)[2] )
}

main()