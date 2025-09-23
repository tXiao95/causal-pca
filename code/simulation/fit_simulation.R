library(data.table)
library(here)
library(SuperLearner)
library(tmle)
library(hal9001)
library(MAVE)

source(here("R/simulate_data.R"))
source(here("R/pCCA.R"))
source(here("R/csPCA.R"))
source(here("R/gcomp.R"))

#' Projection into the d-dim subspace of R^p
project <- function(P){
  P %*% solve(t(P) %*% P) %*% t(P)
}

predict_mu <- function(beta, Y, X, C){
  Z      <- X %*% beta
  mu_hat <- gcomp(Y, Z, C, SL.library = c("SL.glm", 
                                          "tmle.SL.dbarts2", 
                                          "SL.glmnet"))
  return(mu_hat)
}

# Main code ---------------------------------------------------------------
main <- function(){
  #tables <- lapply(c(100, 200, 400, 800, 1600, 3200, 6400), function(n){
  tables <- lapply(c(100, 200, 400, 800, 1600), function(n){
  #tables <- lapply(c(100), function(n){
    # For reproducibility - set seed as Task ID
    message("Simulation experiment with sample size ", n)
    set.seed(TASK_ID)
    
    # Generate data
    sim <- simulate_data(n = n, p=20, q=10, 
                         h_Z_coef = 0.1,
                         g_C_coef = 0.1, 
                         interaction_coef = 0.1, 
                         rho = 0.8, var_scale = 3)
    
    # Data
    Y <- sim$Y
    X <- sim$X
    C <- sim$C
    
    # Desired beta, projection 
    beta   <- sim$beta
    P_beta <- sim$P_beta
    
    # True mu_X with MAVE. 
    time.truth <- system.time({
      truth <- MAVE::mave(sim$mu_X ~ X, method = "meanMAVE")
    })
    P_truth         <- project( truth$dir[[2]] )
    dhat_truth      <- MAVE::mave.dim(truth)$dim.min
    muhat_truth     <- predict_mu(truth$dir[[2]], Y, X, C)
    mu_Z_mse_truth  <- mean( (muhat_truth - sim$mu_X)^2 )
    
    
    # csPCA using MAVE 
    time.cs <- system.time({
      cs   <- csPCA(Y, X, C, 
                    L = 1,
                    mu_fun = gcomp, 
                    mu_args = list(SL.library = c("SL.glm", 
                                                  "tmle.SL.dbarts2",
                                                  "SL.glmnet" 
                                                )
                                   ))
    })
    P_cs        <- project( cs$mave$dir[[2]] )
    dhat_cs     <- MAVE::mave.dim( cs$mave )$dim.min
    muhat_cs    <- predict_mu(cs$mave$dir[[2]], Y, X, C)
    mu_Z_mse_cs <- mean( (muhat_cs - sim$mu_X)^2 )
    mu_X_mse_cs <- mean( (cs$mu_X - sim$mu_X )^2 )
    
    # csPCA (L = 5) cross-fitting with five folds using MAVE 
    time.cs_cf <- system.time({
      cs_cf   <- csPCA(Y, X, C, 
                    L = 5,
                    mu_fun = gcomp, 
                    mu_args = list(SL.library = c("SL.glm",
                                                  "tmle.SL.dbarts2",
                                                  "SL.glmnet" 
                                                )
                                   ))
    })
    P_cs_cf        <- project( cs_cf$mave$dir[[2]] )
    dhat_cs_cf     <- MAVE::mave.dim( cs_cf$mave )$dim.min
    muhat_cs_cf    <- predict_mu(cs_cf$mave$dir[[2]], Y, X, C)
    mu_Z_mse_cs_cf <- mean( (muhat_cs_cf - sim$mu_X)^2 )
    mu_X_mse_cs_cf <- mean( (cs_cf$mu_X - sim$mu_X )^2 )
    
    # pCCA
    time.pcca <- system.time({
      pcca <- pCCA(Y, X, C)
    })
    P_pcca <- pcca[, 1:2, drop = FALSE] |> project()
    
    # Regular MAVE (as SDR benchmark)
    time.sdr <- system.time({
      reg   <- MAVE::mave(Y ~ X, method = "meanMAVE")
    })
    P_reg        <- project( reg$dir[[2]] )
    dhat_reg     <- MAVE::mave.dim( reg )$dim.min
    muhat_reg    <- predict_mu( reg$dir[[2]], Y, X, C)
    mu_Z_mse_reg <- mean( (muhat_reg - sim$mu_X)^2 )
    mu_X_mse_reg <- mean( (Y - sim$mu_X)^2 )
    
    # PCA
    time.pca <- system.time({
      pca <- prcomp(sim$X, center = TRUE, scale. = FALSE)
    })
    P_pca <- pca$rotation[, 1:2, drop = FALSE] |> project()
    
    # Frobenius norm
    error.truth.f <- norm(P_truth - P_beta, "F")
    error.cs.f    <- norm(P_cs - P_beta, "F")
    error.cs_cf.f <- norm(P_cs_cf - P_beta, "F")
    error.sdr.f   <- norm(P_reg - P_beta, "F")
    error.pca.f   <- norm(P_pca - P_beta, "F")
    error.pcca.f  <- norm(P_pcca - P_beta, "F")
    
    # Spectral 2-norm
    error.truth.2 <- norm(P_truth - P_beta, "2")
    error.cs.2    <- norm(P_cs - P_beta, "2")
    error.cs_cf.2 <- norm(P_cs_cf - P_beta, "2")
    error.sdr.2   <- norm(P_reg - P_beta, "2")
    error.pca.2   <- norm(P_pca - P_beta, "2")
    error.pcca.2  <- norm(P_pcca - P_beta, "2")
    
    # Organize estimates into errors data.table
    errors <- data.table::CJ(ID = TASK_ID, 
                             n = n,
                             error_type = c("2-norm", "Frobenius"), 
                             method = c("Truth", "CS","CS-CF" ,"PCA", "SDR", "pCCA"))
    
    errors$value <- c(error.cs.2, error.cs_cf.2, error.pca.2, error.sdr.2, error.truth.2, error.pcca.2,
                      error.cs.f, error.cs_cf.f, error.pca.f, error.sdr.f, error.truth.f, error.pcca.f)
    
    errors$time <- rep(c(
                         time.cs["elapsed"], 
                         time.cs_cf["elapsed"],
                         time.pca["elapsed"], 
                         time.sdr["elapsed"], 
                         time.truth["elapsed"],
                         time.pcca["elapsed"]), 2 
    )
    
    errors[, mu_X_mse := ifelse(method == "CS", mu_X_mse_cs, 
                                ifelse(method == "CS-CF", mu_X_mse_cs_cf, 
                                       ifelse(method == "Truth", 0, 
                                              ifelse(method == "SDR", mu_X_mse_reg, NA))))]
    errors[, dhat := ifelse(method == "CS", dhat_cs, 
                            ifelse(method == "CS-CF", dhat_cs_cf, 
                                   ifelse(method == "Truth", dhat_truth, 
                                          ifelse(method == "SDR", dhat_reg, NA))))]
    errors[, mu_Z_mse := ifelse(method == "CS", mu_Z_mse_cs, 
                            ifelse(method == "CS-CF", mu_Z_mse_cs_cf, 
                                   ifelse(method == "Truth", mu_Z_mse_truth, 
                                          ifelse(method == "SDR", mu_Z_mse_reg, NA))))]
    
    return(errors)
  }) |> rbindlist()
  
  filename <- paste0("sim-", sprintf("%03d", TASK_ID))
  #data.table::fwrite(tables, file = here(paste0("outputs/simulation/enar-student-paper2_low-rho_p6/", filename, ".csv")) )
  data.table::fwrite(tables, file = here(paste0("outputs/simulation/enar-student-paper2_high-rho_p20/", filename, ".csv")) )
  #data.table::fwrite(tables, file = here(paste0("outputs/simulation/enar-student-paper2_high-rho_p6/", filename, ".csv")) )
  #data.table::fwrite(tables, file = here(paste0("outputs/simulation/enar-student-paper2_low-rho_p20/", filename, ".csv")) )
}

# Arguments for main ------------------------------------------------------
if(interactive()){
  # Also the random seed
  TASK_ID  <- 4
} else{
  TASK_ID    <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

main()