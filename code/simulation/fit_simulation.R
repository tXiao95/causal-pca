library(data.table)
library(here)
library(parallel)
library(SuperLearner)
library(tmle)
library(MAVE)

source(here("R/simulate_data.R"))
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

# SL library
SL.lib <- c("SL.glm", 
            "SL.glmnet", 
            "SL.earth", 
            "SL.xgboost", 
            "SL.ksvm",
            "SL.nnet")
SL.lib <- c("SL.glm", "SL.gam", "SL.xgboost")

N_vector <- c(100, 500, 1000, 5000, 10000)
N_vector <- c(100, 500, 1000)
N_vector <- c(100)

sdr_methods <- c("PCA", "PCCA", "MAVE", "EE", 
                 "RA-MAVE", "DR-MAVE", "PO-MAVE", "RP-MAVE",
                 "RA-EE", "DR-EE", "PO-EE", "RP-EE")

mean_methods <- c("mu_X", "mu_Z")

# Simulate data
# Fit all SDR methods to get different dhat and beta
# 1. Given  a) runtime b) beta c) dhat

# Main code ---------------------------------------------------------------
main <- function(){
  tables <- lapply(N_vector, function(n){
    message("Simulation experiment with sample size ", n)
    set.seed(TASK_ID)
    
    # Simulate Data -----------------------------------------------------------
     # sim <- simulate_data(n = n, p=10, q=5, 
     #                      h_Z_coef = 0.01,
     #                      g_C_coef = .1, 
     #                      interaction_coef = 0, 
     #                      Z1_coef = 10, 
     #                      Z2_coef = 10, 
     #                      Z12_coef = 0,
     #                      rho = 0.5, var_scale = 5)
    
    #sim <- simulate_data(n = n)
    sim <- simulate_causal_sdr(n = n)
    
    #sim <- sim_nabi_case1_p6(n = n)
    
    # Data
    Y <- sim$Y
    X <- sim$X; p  <- ncol(X)
    Z <- sim$Z; d0 <- ncol(Z)
    C <- sim$C; q  <- ncol(C)
    
    # Desired beta, projection 
    beta0   <- sim$beta_true
    P_beta0 <- Pi(beta0)
    d0      <- sim$d
    
    # PCA ---------------------------------------------------------------------
    message("PCA")
    time_pca <- system.time({
      pca <- prcomp(sim$X, center = TRUE, scale. = FALSE)
    })
    beta_pca  <- pca$rotation[, 1:d0, drop = FALSE]
    error_pca <- Delta(beta0, beta_pca, type = "F")

    # pCCA --------------------------------------------------------------------
    message("pCCA")
    time_pcca <- system.time({
      pcca <- pCCA(Y, X, C)
    })
    beta_pcca  <- pcca[, 1:d0, drop = FALSE]
    error_pcca <- Delta(beta0, beta_pcca, type = "F")

    # Regular MAVE ------------------------------------------------------------
    message("MAVE")
    time_MAVE <- system.time({
      reg_MAVE <- MAVE::mave(Y ~ X, method = "meanMAVE")
    })
    beta_MAVE_d0    <- reg_MAVE$dir[[d0]]
    error_MAVE_d0   <- Delta(beta0, beta_MAVE_d0, type = "F")
    dhat_MAVE       <- MAVE::mave.dim( reg_MAVE )$dim.min
    beta_MAVE_dhat  <- reg_MAVE$dir[[dhat_MAVE]]
    error_MAVE_dhat <- Delta(beta0, beta_MAVE_dhat)

    # Efficient Estimation: EE ------------------------------------------------
    message("EE")
    time_EE <- system.time({
      beta_EE_d0 <- run_efficient_estimator(X, Y, beta_init = beta_MAVE_d0)
    })
    beta_EE_dhat <- run_efficient_estimator(X, Y, beta_init = beta_MAVE_dhat)
    error_EE_d0   <- Delta(beta0, beta_EE_d0, type = "F")
    error_EE_dhat <- Delta(beta0, beta_EE_dhat, type = "F")

    # RA-MAVE -----------------------------------------------------------------
    message("RA-MAVE")
    time_RA_MAVE <- system.time({
      obj_RA_MAVE <- csMAVE(Y = Y, X = X, C = C, method = "RA", 
                            args_outcome = list(SL.lib = SL.lib))
    })
    beta_RA_MAVE_d0    <- obj_RA_MAVE$mave_fit$dir[[d0]]
    error_RA_MAVE_d0   <- Delta(beta0, beta_RA_MAVE_d0)
    dhat_RA_MAVE       <- obj_RA_MAVE$d_hat
    beta_RA_MAVE_dhat  <- obj_RA_MAVE$mave_fit$dir[[dhat_RA_MAVE]]
    error_RA_MAVE_dhat <- Delta(beta0, beta_RA_MAVE_dhat)
    
    # RA-EE -------------------------------------------------------------------
    message("RA-EE")
    time_RA_EE <- system.time({
      beta_RA_EE_d0 <- run_efficient_estimator(X = X, Y = obj_RA_MAVE$new_data$new_Y, 
                                            beta_init = beta_RA_MAVE_d0)
    })
    beta_RA_EE_dhat <- run_efficient_estimator(X = X, Y = obj_RA_MAVE$new_data$new_Y, 
                                            beta_init = beta_RA_MAVE_dhat)
    error_RA_EE_d0   <- Delta(beta0, beta_RA_EE_d0, type = "F")
    error_RA_EE_dhat <- Delta(beta0, beta_RA_EE_dhat, type = "F")

    # DR-MAVE -----------------------------------------------------------------
    message("DR-MAVE")
    time_DR_MAVE <- system.time({
      obj_DR_MAVE <- csMAVE(Y = Y, X = X, C = C, method = "DR", 
                            args_outcome = list(SL.lib = SL.lib),
                            args_gps = list(SL.lib = SL.lib))
    })
    beta_DR_MAVE_d0    <- obj_DR_MAVE$mave_fit$dir[[d0]]
    error_DR_MAVE_d0   <- Delta(beta0, beta_DR_MAVE_d0)
    dhat_DR_MAVE       <- obj_DR_MAVE$d_hat
    beta_DR_MAVE_dhat  <- obj_DR_MAVE$mave_fit$dir[[dhat_DR_MAVE]]
    error_DR_MAVE_dhat <- Delta(beta0, beta_DR_MAVE_dhat)
    
    # DR-EE -------------------------------------------------------------------
    message("DR-EE")
    time_DR_EE <- system.time({
      beta_DR_EE_d0   <- run_efficient_estimator(X = X, Y = obj_DR_MAVE$new_data$new_Y, 
                                            beta_init = beta_DR_MAVE_d0)
    })
    beta_DR_EE_dhat <- run_efficient_estimator(X = X, Y = obj_DR_MAVE$new_data$new_Y, 
                                            beta_init = beta_DR_MAVE_dhat)
    error_DR_EE_d0   <- Delta(beta0, beta_DR_EE_d0, type = "F")
    error_DR_EE_dhat <- Delta(beta0, beta_DR_EE_dhat, type = "F")
    
    # PO-MAVE -----------------------------------------------------------------
    message("PO-MAVE")
    time_PO_MAVE <- system.time({
      obj_PO_MAVE <- csMAVE(Y = Y, X = X, C = C, method = "PO",
                            args_outcome = list(SL.lib = SL.lib),
                            args_gps = list(SL.lib = SL.lib))
    })
    beta_PO_MAVE_d0    <- obj_PO_MAVE$mave_fit$dir[[d0]]
    error_PO_MAVE_d0   <- Delta(beta0, beta_PO_MAVE_d0)
    dhat_PO_MAVE       <- obj_PO_MAVE$d_hat
    if (dhat_PO_MAVE >= 1L) {
      beta_PO_MAVE_dhat  <- obj_PO_MAVE$mave_fit$dir[[dhat_PO_MAVE]]
      error_PO_MAVE_dhat <- Delta(beta0, beta_PO_MAVE_dhat)
    } else {
      beta_PO_MAVE_dhat  <- NULL
      error_PO_MAVE_dhat <- NA_real_   # or a penalty value; see below
    }
    
    # PO-EE -------------------------------------------------------------------
    message("PO-EE")
    time_PO_EE <- system.time({
      beta_PO_EE_d0 <- run_efficient_estimator(X = X, Y = obj_PO_MAVE$new_data$new_Y, 
                                            beta_init = beta_PO_MAVE_d0)
    })
    error_PO_EE_d0   <- Delta(beta0, beta_PO_EE_d0, type = "F")
    if (dhat_PO_MAVE >= 1L) {
      beta_PO_EE_dhat <- run_efficient_estimator(X = X, Y = obj_PO_MAVE$new_data$new_Y, 
                                            beta_init = beta_PO_MAVE_dhat)
      error_PO_EE_dhat <- Delta(beta0, beta_PO_EE_dhat, type = "F")
    } else {
      beta_PO_EE_dhat  <- NULL
      error_PO_EE_dhat <- NA_real_   # or a penalty value; see below
    }
    
    # RP-MAVE -----------------------------------------------------------------
    message("RP-MAVE")
    time_RP_MAVE <- system.time({
      obj_RP_MAVE <- csMAVE(Y = Y, X = X, C = C, method = "RP", 
                            args_C = list(SL.lib = SL.lib))
    })
    beta_RP_MAVE_d0    <- obj_RP_MAVE$mave_fit$dir[[d0]]
    error_RP_MAVE_d0   <- Delta(beta0, beta_RP_MAVE_d0)
    dhat_RP_MAVE       <- obj_RP_MAVE$d_hat
    beta_RP_MAVE_dhat  <- obj_RP_MAVE$mave_fit$dir[[dhat_RP_MAVE]]
    error_RP_MAVE_dhat <- Delta(beta0, beta_RP_MAVE_dhat)
    
    # RP-EE -------------------------------------------------------------------
    message("RP-EE")
    time_RP_EE <- system.time({
      beta_RP_EE_d0 <- run_efficient_estimator(X = X, Y = obj_RP_MAVE$new_data$new_Y, 
                                            beta_init = beta_RP_MAVE_d0)
    })
    beta_RP_EE_dhat <- run_efficient_estimator(X = X, Y = obj_RP_MAVE$new_data$new_Y, 
                                            beta_init = beta_RP_MAVE_dhat)
    error_RP_EE_d0 <- Delta(beta0, beta_RP_EE_d0, type = "F")
    error_RP_EE_dhat <- Delta(beta0, beta_RP_EE_dhat, type = "F")

    # Oracle MAVE ------------------------------------------------------------------
    message("Oracle")
    time_oracle_MAVE <- system.time({
      obj_oracle_MAVE <- MAVE::mave(sim$mu_X ~ X, method = "meanOPG")
    })
    beta_oracle_MAVE_d0    <- obj_oracle_MAVE$dir[[d0]]
    error_oracle_MAVE_d0   <- Delta(beta0, beta_oracle_MAVE_d0, type = "F")
    dhat_oracle_MAVE       <- MAVE::mave.dim( obj_oracle_MAVE )$dim.min
    beta_oracle_MAVE_dhat  <- obj_oracle_MAVE$dir[[dhat_oracle_MAVE]]
    error_oracle_MAVE_dhat <- Delta(beta0, beta_oracle_MAVE_dhat)

    # Oracle-EE ------------------------------------------------
    message("EE")
    time_oracle_EE <- system.time({
      beta_oracle_EE_d0 <- run_efficient_estimator(X, Y, beta_init = beta_oracle_MAVE_d0)
    })
    beta_oracle_EE_dhat <- run_efficient_estimator(X, Y, beta_init = beta_oracle_MAVE_dhat)
    error_oracle_EE_d0 <- Delta(beta0, beta_oracle_EE_d0, type = "F")
    error_oracle_EE_dhat <- Delta(beta0, beta_oracle_EE_dhat, type = "F")
    
    # Organize estimates into errors data.table
    errors <- data.table::CJ(ID = TASK_ID, 
                             n = n,
                             error_type = c("d0", "dhat"),
                             method = c("PCA","pCCA" ,"MAVE", "EE",
                                        "RA-MAVE", "RA-EE", 
                                        "DR-MAVE", "DR-EE",
                                        "PO-MAVE", "PO-EE",
                                        "RP-MAVE", "RP-EE",
                                        "Oracle-MAVE", "Oracle-EE"),
                              sorted = FALSE)
    
    errors$frob_norm <- c(error_pca, error_pcca, error_MAVE_d0, error_EE_d0, 
                          error_RA_MAVE_d0, error_RA_EE_d0, 
                          error_DR_MAVE_d0, error_DR_EE_d0,
                          error_PO_MAVE_d0, error_PO_EE_d0,
                          error_RP_MAVE_d0, error_RP_EE_d0,
                          error_oracle_MAVE_d0, error_oracle_EE_d0, 
                          error_pca, error_pcca, error_MAVE_dhat, error_EE_dhat,
                          error_RA_MAVE_dhat, error_RA_EE_dhat, 
                          error_DR_MAVE_dhat, error_DR_EE_dhat,
                          error_PO_MAVE_dhat, error_PO_EE_dhat,
                          error_RP_MAVE_dhat, error_RP_EE_dhat,
                          error_oracle_MAVE_dhat, error_oracle_EE_dhat)
    
    errors$time <- rep(c(
                         time_pca["elapsed"], 
                         time_pcca["elapsed"],
                         time_MAVE["elapsed"], 
                         time_EE["elapsed"], 
                         time_RA_MAVE["elapsed"],
                         time_RA_EE["elapsed"], 
                         time_DR_MAVE["elapsed"],
                         time_DR_EE["elapsed"],
                         time_PO_MAVE["elapsed"],
                         time_PO_EE["elapsed"],
                         time_RP_MAVE["elapsed"],
                         time_RP_EE["elapsed"],
                         time_oracle_MAVE["elapsed"],
                         time_oracle_EE["elapsed"]), 2 
    )
    
    errors$dhat <- rep(c(
      NA_integer_, NA_integer_, dhat_MAVE, NA_integer_, 
      dhat_RA_MAVE, NA_integer_,
      dhat_DR_MAVE, NA_integer_,
      dhat_PO_MAVE, NA_integer_,
      dhat_RP_MAVE, NA_integer_,
      dhat_oracle_MAVE, NA_integer_
    ), 2)
    
    return(errors)
  }) |> rbindlist()
  
  filename <- paste0("sim-", sprintf("%03d", TASK_ID))
  data.table::fwrite(tables, file = here(paste0("outputs/simulation/enar-student-paper2_highZ1/", filename, ".csv")) )
}

# Arguments for main ------------------------------------------------------
if(interactive()){
  # Also the random seed
  TASK_ID  <- 4
} else{
  TASK_ID    <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

main()