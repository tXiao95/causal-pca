library(here)

source(here("R", "gcomp.R"))
source(here("R", "kernel_smooth.R"))
source(here("R", "csPCA.R"))
source(here("R", "doPCA.R"))
source(here("R", "misc.R"))
source(here("R", "sim_additive.R"))
source(here("R", "simdata_nabi2022.R"))
source(here("R", "population_obj.R"))
source(here("R", "est_conditional_moments.R"))

# Main code ---------------------------------------------------------------
main <- function(){
  set.seed(TASK_ID)
  # sim <- sim_additive(n = 500, seed = TASK_ID, f = function(x) x^2, g = function(x) sin(x), 
  #                     confounding_type = "multiplicative" )
  
  # Baseline
  if(scenario == 1){
    sim <- sim_additive(n = n, 
                        seed = TASK_ID, 
                        f = function(x) x, 
                        g = function(x) x, 
                        e.dist = "normal",
                        confounding_type = "additive")
  }
  
  # Nonlinear X
  if(scenario == 2){
    sim <- sim_additive(n = n, 
                        seed = TASK_ID, 
                        f = function(x) x^2, 
                        g = function(x) x, 
                        e.dist = "normal",
                        confounding_type = "additive")
  }
  
  # Nonlinear X, nonlinear C
  if(scenario == 3){
    sim <- sim_additive(n = n, 
                        seed = TASK_ID, 
                        f = function(x) x^2, 
                        g = function(x) sin(x), 
                        e.dist = "normal",
                        confounding_type = "additive")
  }
  
  # Nonlinear X, nonlinear C, Laplace Error
  if(scenario == 4){
    sim <- sim_additive(n = n, 
                        seed = TASK_ID, 
                        f = function(x) x^2, 
                        g = function(x) sin(x), 
                        e.dist = "laplace",
                        confounding_type = "additive")
  }
  
  # Nonlinear X, nonlinear C, Laplace Error, Multiplicative confounding error
  if(scenario == 5){
    sim <- sim_additive(n = n, 
                        seed = TASK_ID, 
                        f = function(x) x^2, 
                        g = function(x) sin(x), 
                        e.dist = "laplace",
                        confounding_type = "multiplicative")
  }
  
  # sim <- simdata_nabi(n = 500, seed = TASK_ID, case = 2)
  
  # # Empirical covariance matrices. Cov(X,Y) is p x 1 and Sigma_XX is p x p.
  Sigma_XY <- cov(sim$X, sim$Y) |> as.numeric()
  Sigma_XX <- cov(sim$X)
  
  # # Associational PCA
  omega.pca <- prcomp(sim$X, center = TRUE, scale. = FALSE)$rotation[,1]
  omega.pls <- unitvec(Sigma_XY)                             # PLS
  omega.cca <- unitvec(solve(Sigma_XX, Sigma_XY))            # CCA
  
  # Associational PCA w/ confounders
  partial_moments <- est_conditional_moments(sim$Y, sim$X, sim$C)
  omega.ppls      <- unitvec(partial_moments$Cov_XY_given_C) |> as.numeric()
  omega.pcca      <- unitvec(solve(partial_moments$Var_X_given_C, partial_moments$Cov_XY_given_C)) |>as.numeric()
  
  # Causal PCA
  # Timing and storing causal PCA methods only
  time.cs <- system.time({
    cs <- csPCA(sim$Y, sim$X, sim$C, maxit = 500, omega0 = omega.pcca)
  })
  
  time.do <- system.time({
    do <- doPCA(sim$Y, sim$X, sim$C, maxit = 500, omega0 = omega.pcca)
  })
  
  time.do.s <- system.time({
    do.s <- doPCA(sim$Y, sim$X, sim$C, maxit = 500, omega0 = omega.pcca, scaled = TRUE)
  })
  
  # Organize estimates
  est <- rbind(sim$beta |> unitvec(),
               cs$omega, do$omega, do.s$omega,
               omega.pca, omega.pls, omega.cca,
               omega.ppls, omega.pcca) |> data.frame()
  
  rownames(est) <- NULL
  colnames(est) <- paste0("X", 1:ncol(sim$X))
  est[est$X1 < 0, ] <- -est[est$X1 < 0, ]
  
  est$method    <- c("Truth", 
                     "cs", "do", "do_scaled", 
                     "PCA", "PLS", "CCA", "pPLS", "pCCA")
  est$ID        <- TASK_ID
  est$i <- i
  est$n <- n
  est$scenario <- scenario
  
  # Convergence
  cs.convergence   <- cs$opt$status
  do.convergence   <- do$opt$status
  do.s.convergence <- do.s$opt$status
  
  convergence <- c(NA, 
                   cs.convergence, do.convergence, do.s.convergence, 
                   NA, NA, NA, NA, NA)
  est$convergence <- convergence
  
  # Runtime
  runtime <- c(NA, 
               time.cs[["elapsed"]], time.do[["elapsed"]], time.do.s[["elapsed"]], 
               NA, NA, NA, NA, NA) 
  est$runtime <- runtime
  
  # Estimation error (by Euclidean norm for the vector)
  Xcols    <- paste0("X", 1:ncol(sim$X))
  X_matrix <- est[, Xcols] |> as.matrix()
  truth    <- X_matrix[1,]
  
  error <- apply(X_matrix, 1, function(omega){
    norm(t(truth - omega), "F")
  })
  est$error <- error
  
  # Final object
  obj <- list(data = sim,
              est = est, 
              cs = cs, do = do, do.s = do.s)
  
  filename <- paste0("sim-", TASK_ID, "_", scenario, "-", i, "-", n)
  
  saveRDS(obj, file = paste0("results/sim_obj/", filename, ".rds"))
  readr::write_csv(est, file = paste0("results/sim_results/", filename, ".csv"))
}

if(interactive()){
  TASK_ID  <- 1
  scenario <- 1
  i        <- 2
  n        <- 200
  # Each experiment should be for the data model. Different sample sizes within each.  
} else{
  TASK_ID    <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) + 1000
  taskmap <- read.csv(here("config/taskmap.csv"))
  scenario <- taskmap[TASK_ID, "scenario"]
  i        <- taskmap[TASK_ID, "i"]
  n        <- taskmap[TASK_ID, "n"]
}

main()