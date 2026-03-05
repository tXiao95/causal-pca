library(data.table)
library(here)
library(parallel)
library(SuperLearner)
library(MAVE)

source(here("R/pCCA.R"))
source(here("R/csMAVE.R"))
source(here("R/Seff.R"))
source(here("R/compute_new_response_and_exposure.R"))

source(here("R/estimate_ERS.R"))
source(here("R/estimate_pseudo_outcomes.R"))
source(here("R/estimate_residualized_pair.R"))

source(here("R/nuisance_outcome_regression.R"))
source(here("R/nuisance_outcome_C.R"))
source(here("R/nuisance_gps.R"))

source(here("R/misc.R"))
source(here("R/sims_from_papers.R"))
source(here("R/simulate_data.R"))

# -------------------------------------------------------------------------
# Well-Specified Fitters for simulate_causal_sdr DGP
# -------------------------------------------------------------------------

# 1. Outcome Regression: E[Y | X, C]
# The DGP features non-linearities (tanh, pnorm) and interactions.
# SL.earth (MARS) and SL.gam are excellent for smooth, bounded non-linear surfaces.
well_specified_outcome_fitter <- function(Y, XC_df, ...) {
  SL_outcome_fitter(Y, XC_df, SL.lib = c("SL.glm", "SL.glmnet", "SL.gam", 
                                         "SL.earth", "SL.xgboost"), ...)
}

# 2. GPS Model: f(X | C)
# The DGP strictly generates X | C as a multivariate normal with a linear mean.
# A linear MVN model is perfectly specified and lightning fast.
well_specified_gps_fitter <- function(X, C, ...) {
  mvn_fitter(X, C, method_gps = "linear", ...)
}

# 3. C_fitter: E[Y | C] and E[X | C]
# E[X|C] is linear, but E[Y|C] is non-linear after marginalizing over X.
# Including SL.glm and SL.earth allows the ensemble to adapt appropriately.
well_specified_C_fitter <- function(target, C_df, ...) {
  SL_nuisance_fitter(target, C_df, SL.lib = c("SL.glm", "SL.glmnet", "SL.earth"), ...)
}

# -------------------------------------------------------------------------
# Evaluation Module for Parallel Processing
# -------------------------------------------------------------------------

evaluate_method_group <- function(group, sim, n, task_id) {
  
  Y <- sim$Y; X <- sim$X; C <- sim$C
  beta0 <- sim$beta_true
  d0 <- sim$d
  
  # Helper to structure the output data.table
  format_res <- function(method_name, frob_d0, frob_dhat, time_elapsed, dhat_val) {
    data.table(
      ID = task_id, n = n, method = method_name, time = time_elapsed, dhat = dhat_val,
      error_type = c("d0", "dhat"), frob_norm = c(frob_d0, frob_dhat)
    )
  }
  
  results <- list()
  
  if (group == "Base") {
    message("Running baseline methods")
    # PCA
    t_pca <- system.time({ pca <- prcomp(X, center = TRUE, scale. = FALSE) })
    b_pca <- pca$rotation[, 1:d0, drop = FALSE]
    results[[1]] <- format_res("PCA", Delta(beta0, b_pca, "F"), NA, t_pca["elapsed"], NA)
    
    # pCCA
    t_pcca <- system.time({ b_pcca <- pCCA(Y, X, C)[, 1:d0, drop = FALSE] })
    results[[2]] <- format_res("pCCA", Delta(beta0, b_pcca, "F"), NA, t_pcca["elapsed"], NA)
    
    # MAVE
    t_mave <- system.time({ reg_MAVE <- MAVE::mave(Y ~ X, method = "meanMAVE") })
    dhat_m <- MAVE::mave.dim(reg_MAVE)$dim.min
    err_mave_d0 <- Delta(beta0, reg_MAVE$dir[[d0]], "F")
    err_mave_dhat <- if(dhat_m >= 1) Delta(beta0, reg_MAVE$dir[[dhat_m]], "F") else NA
    results[[3]] <- format_res("MAVE", err_mave_d0, err_mave_dhat, t_mave["elapsed"], dhat_m)
    
    # EE
    t_ee <- system.time({ b_ee_d0 <- run_efficient_estimator(X, Y, beta_init = reg_MAVE$dir[[d0]]) })
    b_ee_dhat <- if(dhat_m >= 1) run_efficient_estimator(X, Y, beta_init = reg_MAVE$dir[[dhat_m]]) else NA
    err_ee_d0 <- Delta(beta0, b_ee_d0, "F")
    err_ee_dhat <- if(dhat_m >= 1) Delta(beta0, b_ee_dhat, "F") else NA
    results[[4]] <- format_res("EE", err_ee_d0, err_ee_dhat, t_ee["elapsed"], NA)
    
  } else if (group %in% c("RA", "DR", "PO", "RP")) {
    message("Running method ", group)
    # 1. Run the unified csMAVE pipeline once for the family
    t_mave <- system.time({
      obj <- csMAVE(Y = Y, X = X, C = C, method = group,
                    args_compute_new_response = list(
                      outcome_fitter = well_specified_outcome_fitter,
                      gps_fitter     = well_specified_gps_fitter,
                      C_fitter       = well_specified_C_fitter,
                      L              = 5 # 5-fold cross-fitting
                    ),
                    args_MAVE = list(method = "meanMAVE"))
    })
    
    dhat_m <- obj$d_hat
    err_mave_d0 <- Delta(beta0, obj$mave_fit$dir[[d0]], "F")
    err_mave_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) Delta(beta0, obj$mave_fit$dir[[dhat_m]], "F") else NA
    results[[1]] <- format_res(paste0(group, "-MAVE"), err_mave_d0, err_mave_dhat, t_mave["elapsed"], dhat_m)
    
    # 2. Run EE on the generated pseudo-data
    new_Y <- obj$new_data$new_Y; new_X <- obj$new_data$new_X
    t_ee <- system.time({ b_ee_d0 <- run_efficient_estimator(new_X, new_Y, beta_init = obj$mave_fit$dir[[d0]]) })
    b_ee_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) run_efficient_estimator(new_X, new_Y, beta_init = obj$mave_fit$dir[[dhat_m]]) else NA
    
    err_ee_d0 <- Delta(beta0, b_ee_d0, "F")
    err_ee_dhat <- if(!is.na(dhat_m) && dhat_m >= 1) Delta(beta0, b_ee_dhat, "F") else NA
    results[[2]] <- format_res(paste0(group, "-EE"), err_ee_d0, err_ee_dhat, t_ee["elapsed"], NA)
    
  } else if (group == "Oracle") {
    message("Running Oracle")
    t_mave <- system.time({ obj <- MAVE::mave(sim$mu_X ~ X, method = "meanMAVE") })
    dhat_m <- MAVE::mave.dim(obj)$dim.min
    err_mave_d0 <- Delta(beta0, obj$dir[[d0]], "F")
    err_mave_dhat <- if(dhat_m >= 1) Delta(beta0, obj$dir[[dhat_m]], "F") else NA
    results[[1]] <- format_res("Oracle-MAVE", err_mave_d0, err_mave_dhat, t_mave["elapsed"], dhat_m)
    
    t_ee <- system.time({ b_ee_d0 <- run_efficient_estimator(X, Y, beta_init = obj$dir[[d0]]) })
    b_ee_dhat <- if(dhat_m >= 1) run_efficient_estimator(X, Y, beta_init = obj$dir[[dhat_m]]) else NA
    err_ee_d0 <- Delta(beta0, b_ee_d0, "F")
    err_ee_dhat <- if(dhat_m >= 1) Delta(beta0, b_ee_dhat, "F") else NA
    results[[2]] <- format_res("Oracle-EE", err_ee_d0, err_ee_dhat, t_ee["elapsed"], NA)
  }
  
  return(rbindlist(results))
}

# -------------------------------------------------------------------------
# Main Execution Loop
# -------------------------------------------------------------------------

main <- function() {
  message("TASK_ID: ", TASK_ID)
  message("Number of cores: ", N_CORES)
  message("Experiment: ", EXPERIMENT)
  
  # Actual simulation
  N_vector <- c(100, 500, 1000, 2500, 5000)
  groups   <- c("Base", "RA", "DR", "PO", "RP", "Oracle")
  
  # For testing
  #N_vector <- c(100, 500, 1000)
  #N_vector <- c(100)
  #groups   <- c("Base", "RA")
  
  # Retrieve number of cores from SLURM environment for parallelization
  message("Running with ", N_CORES, " cores.")
  
  tables <- lapply(N_vector, function(n) {
    message("\nSimulation experiment with sample size ", n)
    set.seed(TASK_ID)
    
    # Baseline experiment parameters
    if(EXPERIMENT == "baseline"){
      sim <- simulate_causal_sdr(n = n, 
                                 p = 10,
                                 q = 5,
                                 rho_X = 0.7,
                                 causal_conf_strength = 1.0,   # Keeps Z bounded in [-3, 3] grid
                                 spurious_strength = 2.0,      # Keeps the MAVE trap strong
                                 var_scale = 5,                # Keeps the PCA trap strong
                                 signal_multiplier = 2.0,      # Knob 1: ERS Signal strength
                                 noise_sd = 0.5,               # Knob 2: Noise strength (SNR control)
                                 interaction_coef = 10,        # Huge misleading interaction
                                 heteroskedastic = FALSE) 
    }
    
    if(EXPERIMENT == "additive_confounding"){
      sim <- simulate_causal_sdr(n = n, 
                                 p = 10,
                                 q = 5,
                                 rho_X = 0.7,
                                 causal_conf_strength = 1.0,   
                                 spurious_strength = 2.0,      
                                 var_scale = 5,                
                                 signal_multiplier = 2.0,      
                                 noise_sd = 0.5,               
                                 interaction_coef = 0,         # Turning off interaction between X and C
                                 heteroskedastic = FALSE) 
    }
    
    if(EXPERIMENT == "weak_dim_signal"){
      sim <- simulate_weak_dim_signal (n = n,
                                       p = 10,
                                       q = 5,
                                       rho_X = 0.3,
                                       confounding_strength = 3,
                                       spurious_strength = 6,
                                       var_scale = 5,            
                                       interaction_coef = 1,
                                       sigma2 = 0.25,
                                       heteroskedastic = FALSE)
    }
    
    # Run the families in parallel
    res_list <- mclapply(groups, function(g) {
      evaluate_method_group(group = g, sim = sim, n = n, task_id = TASK_ID)
    }, mc.cores = N_CORES)
    
    return(rbindlist(res_list))
  }) |> rbindlist()
  
  # Save Output
  out_dir <- here("outputs/simulation/jasa-initial-submission", EXPERIMENT)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  filename <- paste0("sim-", sprintf("%03d", TASK_ID), ".csv")
  data.table::fwrite(tables, file = file.path(out_dir, filename))
  message("Results saved to ", filename)
}

# Arguments for main ------------------------------------------------------
if(interactive()){
  TASK_ID    <- 100
  N_CORES    <- 1
  EXPERIMENT <- "baseline"
} else{
  TASK_ID    <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  N_CORES    <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
  EXPERIMENT <- as.character( commandArgs(trailingOnly=TRUE)[1] )
}

main()