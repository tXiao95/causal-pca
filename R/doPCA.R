library(nloptr)
library(ManifoldOptim)

doPCA <- function(Y, X, C, d = 2, mu_est = "gcomp", scaled = FALSE, maxit = 5000, verbose = FALSE, omega0 = NULL){
  # n <- nrow(X); p <- ncol(X); q <- ncol(C); d <- 2
  
  ## 0) Ensure numeric matrices and reasonable scale
  X <- as.matrix(X); storage.mode(X) <- "double"
  C <- as.matrix(C); storage.mode(C) <- "double"
  Y <- as.numeric(Y)
  
  # scale to improve conditioning
  Xs <- scale(X); Cs <- scale(C)
  n  <- nrow(Xs); p <- ncol(Xs); d <- 2
  
  ## 1) Fixed CV folds for determinism
  set.seed(1)
  V <- 5
  fold_id <- sample(rep_len(1:V, n))
  validRows <- split(seq_len(n), fold_id)
  
  ## 2) Libraries: fast (phase 1) vs fuller (phase 2)
  SL_fast <- c("SL.glm", "SL.gam", "SL.glmnet")
  SL_full <- c("SL.glm", "SL.gam", "SL.glmnet", "SL.ranger")  # add more later if you want
  
  ## 3) Optional: subsample index for phase 1 objective
  n_sub <- min(n, 2000L)
  sub_id <- if (n > n_sub) sample.int(n, n_sub) else seq_len(n)
  
  ## 4) Objective factory so we can switch settings easily
  gcomp_obj_factory <- function(SL.library, use_subsample = TRUE) {
    function(beta) {
      B <- matrix(beta, nrow = p, ncol = d)       # p×d
      Z <- Xs %*% B                                # n×d
      
      if (use_subsample) {
        i <- sub_id
        mu_hat <- gcomp(Y[i], Z[i, , drop = FALSE], Cs[i, , drop = FALSE],
                        SL.library = SL.library,
                        cvControl  = list(V = V, validRows = lapply(validRows, intersect, y = i)))
      } else {
        mu_hat <- gcomp(Y, Z, Cs,
                        SL.library = SL.library,
                        cvControl  = list(V = V, validRows = validRows))
      }
      
      # scale objective so magnitudes are O(1)
      - var(mu_hat, na.rm = TRUE) / var(Y, na.rm = TRUE)
    }
  }
  
  f_fast <- gcomp_obj_factory(SL_fast, use_subsample = TRUE)
  f_full <- gcomp_obj_factory(SL_full, use_subsample = FALSE)
  
  ## 5) Problem + manifold
  mod  <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
  prob_fast <- new(mod$RProblem, f_fast)
  
  # Good initial value (p×d) and orthonormalize
  mfit <- MAVE::mave(Y ~ Xs, method = "meanMAVE")
  x0   <- as.matrix(mfit$dir)[[d]]
  x0   <- qr.Q(qr(x0))  # ensure Stiefel
  
  mani.defn <- get.stiefel.defn(p, d)
  
  ## 6) Solver params: coarse then refine
  mani.params_fast   <- get.manifold.params(IsCheckParams = FALSE)
  solver.params_fast <- get.solver.params(
    Tolerance     = 1e-3,
    Max_Iteration = 200,
    OutputGap     = 10,
    DEBUG         = 0
  )
  
  mani.params_ref   <- get.manifold.params(IsCheckParams = FALSE)
  solver.params_ref <- get.solver.params(
    Tolerance     = 5e-4,
    Max_Iteration = 600,
    OutputGap     = 10,
    DEBUG         = 0
  )
  
  ## 7) Phase 1: fast, robust method (RBFGS/RCG) with subsampling
  res1 <- manifold.optim(
    prob_fast, mani.defn,
    method        = "RBFGS",             # or "RCG"
    mani.params   = mani.params_fast,
    solver.params = solver.params_fast,
    x0            = x0
  )
  
  ## 8) Phase 2: refine with fuller library & trust-region
  prob_full <- new(mod$RProblem, f_full)
  
  res2 <- manifold.optim(
    prob_full, mani.defn,
    method        = "RTRSR1",
    mani.params   = mani.params_ref,
    solver.params = solver.params_ref,
    x0            = res1$xopt
  )
  
  res <- res2
  
  
}
