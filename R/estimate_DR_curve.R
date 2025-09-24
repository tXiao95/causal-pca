library(KernSmooth)
library(SuperLearner)
library(tmle)

#' @description Implementation of Kennedy et al. (2017) for estimate dose-response curve for 1D continuous treatment.
#' @param Y (n x 1) vector of outcome/response
#' @param Z (n x 1) vector of treatment
#' @param C (n x q) matrix of confounders
#' @param Z.new (m x 1) vector of treatment values to evaluate E(Y(z)) on
#' @param SL.library SuperLearner Library for all nuisance functions
#' @param bw.range
#' @param tol
#' 
#' @importFrom SuperLearner SuperLearner 
#' @importFrom KernSmooth locpoly
#' @export

estimate_DR_curve <- function(Y, Z, C,
                              Z.new = NULL,
                              SL.library = c(
                                             "SL.earth",
                                             "tmle.SL.dbarts2",
                                             "SL.glmnet",
                                             "SL.glm"),
                              bw.range   = c(0.01, 50),
                              tol        = 0.01)
{
  ## ------------------------------------------------------------------ ##
  ## 0.  Basic checks & helpers
  ## ------------------------------------------------------------------ ##
  
  C <- data.frame(C)
  n <- length(Z)
  if(is.null(Z.new)){
    Z.new <- Z
  }
  m <- length(Z.new)
  if (nrow(C) != n || length(Y) != n){
    stop("Lengths of Y, Z, and rows of C must be equal")
  }
  
  ## simple spline interpolator used several times
  approx_fn <- function(x, y, x_out){
    predict(stats::smooth.spline(x, y), x = x_out)$y
  }
  
  ## ------------------------------------------------------------------ ##
  ## 1.  Build prediction designs:  (C_i , Z_j)  for all i,j
  ## ------------------------------------------------------------------ ##
  
  # Original covariates and treatment
  CZ          <- cbind(C, Z = Z)                               # n rows
  # Replicate the covariates and treatment 'm' times for evaluation points.
  CZ.pred     <- cbind( C[rep(seq_len(n), m), , drop = FALSE],
                        Z = rep(Z.new, rep(n, m)) )                  # n*m rows
  
  # Combine original and new data into one
  CZ.new      <- rbind(CZ, CZ.pred)                            # n + n^2
  # Isolate the covariates pieces
  C.new       <- CZ.new[, -dim(CZ.new)[2]] 
  
  ## ------------------------------------------------------------------ ##
  ## 2.  Nuisance functions for Z | C
  ## ------------------------------------------------------------------ ##
  ## 2.1  π̂(C) = E[Z|C]
  pimod   <- SuperLearner::SuperLearner(Y = Z, 
                                        X = C,
                                        SL.library = SL.library,
                                        newX = C.new)
  pimod.vals <- pimod$SL.predict                                # n + n^2
  
  ## 2.2  σ̂²(C) = Var[Z|C]a
  sq.res     <- (Z - pimod.vals[seq_len(n)])^2
  log.sq.res <- log( pmax(sq.res, 1e-8) )
  pi2mod     <- SuperLearner::SuperLearner(Y = log.sq.res, 
                                       X = C,
                                       SL.library = SL.library,
                                       newX = C.new)
  pi2mod.vals <- exp( pi2mod$SL.predict )                               # n + n^2 )
  
  ## ------------------------------------------------------------------ ##
  ## 3.  Conditional density  f̂_{Z|C}
  ## ------------------------------------------------------------------ ##
  # (Z - E(Z|C)) / sqrt( Var(Z|C) ) Standardize to mean 0, var 1.  
  Z.std    <- (CZ.new[, "Z"] - pimod.vals) / sqrt(pi2mod.vals)
  
  # Taking density over the Z (standardized)
  dens     <- stats::density(Z.std[seq_len(n)])                # kernel on sample resid.
  
  # Approximating the density function and evaluating at Z (standardized)
  pihat.vals <- approx_fn(dens$x, dens$y, Z.std)
  pihat      <- pihat.vals[seq_len(n)]
  pihat.mat  <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = m)
  
  # Variance
  varpihat     <- approx_fn(Z.new, apply(pihat.mat, 2, mean), Z)
  varpihat.mat <- matrix(rep(apply(pihat.mat, 2, mean), n), byrow = TRUE, nrow = n)
  
  ## ------------------------------------------------------------------ ##
  ## 4.  Outcome regression  μ̂(C,Z)
  ## ------------------------------------------------------------------ ##
  mumod   <- SuperLearner::SuperLearner(Y = Y,
                                         X = CZ,
                                         SL.library = SL.library,
                                         newX = CZ.new)
  muhat.vals <- mumod$SL.predict                                # n + n^2
  muhat      <- muhat.vals[seq_len(n)]
  muhat.mat  <- matrix( muhat.vals[-seq_len(n)], nrow = n, ncol = m )
  mhat       <- approx_fn(Z.new, apply(muhat.mat, 2, mean), Z)
  mhat.mat   <- matrix( rep(apply(muhat.mat, 2, mean), n), byrow = TRUE, nrow = n)  # m̂(Z_j) = E_C μ̂(C,Z_j)
  
  ## ------------------------------------------------------------------ ##
  ## 5.  Pseudo-outcome ξ̂_i
  ## ------------------------------------------------------------------ ##
  pseudo.out <- (Y - muhat) / (pihat/varpihat) + mhat
  
  ## ------------------------------------------------------------------ ##
  ## 6.  Bandwidth selection by leave-one-out CV (local-linear)
  ## ------------------------------------------------------------------ ##
  kern <- function(x) stats::dnorm(x)
  
  ## hat-value function w(a; bw) for *all* a in Z
  w.fn <- function(bw) {
    w.vec <- numeric(0)
    for (a.val in Z) {
      z.std   <- (Z - a.val) / bw
      k.std   <- kern(z.std) / bw
      m1      <- mean(k.std)
      m2      <- mean(z.std * k.std)
      m3      <- mean(z.std^2 * k.std)
      w.a     <- m3 * kern(0) / bw /
        (m1 * m3 - m2^2)
      w.vec   <- c(w.vec, w.a / n)   # divide by n: hat-matrix diagonal element
    }
    w.vec
  }
  
  ## function returning diagonal of hat matrix at sample points
  hatvals <- function(bw){approx(Z, w.fn(bw), xout=Z)$y }
  
  ## local-linear smoother of pseudo.out at sample points
  cts.eff <- function(out, bw) {
    approx(KernSmooth::locpoly(x=Z, y=out, bandwidth = bw), xout=Z)$y
  }
  
  ## note: choice of bandwidth range depends on problem
  h.opt <- optimize( f = function(h){
    hats <- hatvals(h)
    mean( ((pseudo.out - cts.eff(out = pseudo.out, bw=h)) / (1-hats))^2)
  }, 
  interval = bw.range, 
  tol = tol)$minimum
  
  ## ------------------------------------------------------------------ ##
  ## 7.  Final dose–response estimate at original Z values
  ## ------------------------------------------------------------------ ##
  est <- approx(locpoly(Z, pseudo.out, bandwidth = h.opt), xout = Z.new)$y
  

  # Estimate sandwich-style pointwise confidence band -----------------------
  
  ## ------------------------------------------------------------------ ##
  ## 8.  Return
  ## ------------------------------------------------------------------ ##
  list(bandwidth  = h.opt,
       Z          = Z,
       DR_curve   = est,
       nuisance   = list(pihat   = pihat,
                         varpihat = varpihat,
                         muhat   = muhat,
                         mhat    = mhat,
                         pseudo.out   = pseudo.out))
}

