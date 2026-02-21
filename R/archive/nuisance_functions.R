library(SuperLearner)

# Main function --------------------------------------------------------------------
#' Regress each column of phi(X, C) on C using SuperLearner
#'
#' @param phi_XC n x (p*d) numeric matrix of feature/basis evaluations.
#'   Column j corresponds to one element of a (p-by-d) basis expansion (flattened).
#' @param C n x q confounder matrix/data.frame.
#' @param ... passed to SuperLearner::SuperLearner (e.g., SL.library, family, cvControl).
#'
#' @return n x (p*d) matrix of out-of-fold predictions \hat E[phi(X,C) | C].
fit_phi_given_C <- function(phi_XC, C, ...) {
  phi_XC <- as.matrix(phi_XC)
  C <- as.data.frame(C)
  
  stopifnot(
    nrow(phi_XC) == nrow(C),
    is.numeric(phi_XC)
  )
  
  n  <- nrow(phi_XC)
  pd <- ncol(phi_XC)
  
  phi_hat_given_C <- matrix(NA_real_, nrow = n, ncol = pd)
  
  for (j in seq_len(pd)) {
    sl_fit <- SuperLearner::SuperLearner(Y = phi_XC[, j], X = C, ...)
    phi_hat_given_C[, j] <- sl_fit$SL.predict
  }
  
  phi_hat_given_C
}

#' @param Y Numeric vector of length \eqn{n}. Outcome.
#' @param X Numeric matrix of dimension \eqn{n \times p}. Exposure/covariate matrix.
#' @param C Data frame or matrix of dimension \eqn{n \times q}. Confounders used to fit
#'   \eqn{E(\phi \mid C)}.
#' @param EY_given_XC Numeric vector of length \eqn{n}. Estimated outcome regression
#'   \eqn{\hat E(Y \mid X, C)} evaluated at each observation.
#' @param pX Numeric vector of length \eqn{n}. Stabilization weight numerator (e.g.,
#'   \eqn{\hat p(X)} or analogous quantity depending on your weighting scheme).
#' @param pXC Numeric vector of length \eqn{n}. Stabilization weight denominator / propensity
#'   component (e.g., \eqn{\hat p(X \mid C)} or propensity score term). Must be nonzero.
#' @param m_hat Numeric vector of length \eqn{n}. Estimated regression function evaluated at
#'   the index, \eqn{\hat m(\beta^\top X_i)}.
#' @param mprime_hat Numeric matrix of dimension \eqn{n \times d}. Estimated derivative(s)
#'   of \eqn{m} with respect to the index, \eqn{\hat m'(\beta^\top X_i)}.
#' @param sigma2_hat Numeric vector of length \eqn{n}. Estimated conditional variance
#'   \eqn{\hat\sigma^2(X_i)}. Must be positive.
#' @param EX_over_sigma2 Numeric matrix of dimension \eqn{n \times p}. Estimated conditional
#'   expectation \eqn{\widehat{E\{X/\hat\sigma^2(X) \mid \beta^\top X_i\}}}.
#' @param E_inv_sigma2 Numeric vector of length \eqn{n}. Estimated conditional expectation
#'   \eqn{\widehat{E\{1/\hat\sigma^2(X) \mid \beta^\top X_i\}}}. Must be nonzero.
#'   
#' @return A list with:
#' \describe{
#'   \item{S_bar_vec}{Numeric vector of length \eqn{p d}. Sample mean
#'     \eqn{\bar S = n^{-1}\sum_i S_i}.}
#'   \item{S2_bar}{Numeric matrix of dimension \eqn{(p d) \times (p d)}. Empirical
#'     second moment \eqn{n^{-1}\sum_i S_i S_i^\top}, computed as
#'     \code{crossprod(S_i_vec) / n}.}
#' }
compute_score_eq7 <- function(
    Y,                   # n 
    X,                   # n x p
    C,                   # n x q
    EY_given_XC,         # n (outcome regression)
    pX,                  # n (stabilisation weights)
    pXC,                 # n (propensity score)
    m_hat,               # n : m̂(β^T X_i)
    mprime_hat,          # n x d : m̂'(β^T X_i)
    sigma2_hat,          # n : σ̂^2(X_i)
    EX_over_sigma2,      # n x p : Ê{X/σ̂^2(X) | β^T X_i}
    E_inv_sigma2         # n : Ê{1/σ̂^2(X) | β^T X_i}
) {
  stopifnot(is.matrix(X))
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(Y) == n,
            length(m_hat) == n,
            nrow(mprime_hat) == n,
            length(sigma2_hat) == n,
            nrow(EX_over_sigma2) == n,
            ncol(EX_over_sigma2) == p,
            length(E_inv_sigma2) == n)
  
  d <- ncol(mprime_hat)
  
  # Propensity weight that is n x 1 vector
  W <- pX / pXC
  
  # p-vector inside brackets for each i
  # g_i = (1/sigma2_hat_i) * bracket_i
  ratio_mat <- EX_over_sigma2 / E_inv_sigma2  # n x p (recycles E_inv_sigma2 by row)
  bracket   <- X - ratio_mat                  # n x p 
  g_mat     <- bracket / sigma2_hat               # n x p (recycles sigma2_hat by row, each row gets same sigma)
  
  # residual part
  r_Y  <- Y - m_hat                              
  r_EY <- EY_given_XC - m_hat
  
  # Pass 1: compute vec(Yterm) and vec(EYterm) which are n x (pd) matrices
  Yterm_vec  <- matrix(0, nrow = n, ncol = p*d)
  phi_vec <- matrix(0, nrow = n, ncol = p*d)
  
  for (i in 1:n) {
    temp            <- (g_mat[i, ] %o% mprime_hat[i, ]) * W[i]  # p x d
    Yterm_vec[i, ]  <- as.vector(temp * r_Y[i])      # vec(p x d)
    phi_vec[i, ]    <- as.vector(temp * r_EY[i])     # vec(p x d)
  }
  
  # Fit E(phi|C) ≈ E(EYterm | C) componentwise, also n x (pd)
  Ephi_given_C_vec <- fit_phi_given_C(phi_XC = phi_vec, C = C, 
                              SL.library = c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger"), 
                              cvControl = list(V = 2))
  
  # Compute score for each i, its sample average, and the cross product. 
  S_i_vec    <- Yterm_vec - phi_vec + Ephi_given_C_vec
  S_bar_vec  <- colMeans(S_i_vec)
  S2_bar     <- crossprod(S_i_vec) / n
  
  return(list(S_bar_vec = S_bar_vec, S2_bar = S2_bar))
}

# Test --------------------------------------------------------------------
set.seed(3)

# Dimensions
n <- 1000
p <- 20
q <- 10
d <- 1

# Basic data
Y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
C <- matrix(rnorm(n * q), n, q)

# Nuisance-type quantities
EY_XC <- rnorm(n)              # outcome regression prediction

pX  <- runif(n, 0.2, 0.8)      # stabilisation weights
pXC <- runif(n, 0.2, 0.8)      # propensity scores

m_hat      <- rnorm(n)                 # m̂(β^T X_i)
mprime_hat <- matrix(rnorm(n * d), n, d)

sigma2_hat <- abs(rnorm(n)) + 0.5      # ensure positive

EX_over_sigma2 <- matrix(rnorm(n * p), n, p)
E_inv_sigma2   <- abs(rnorm(n)) + 0.5  # positive

# Nuisance function groups
# 1. E(Y|X,C)
# 2. p(X) and p(X|C)
# 3. m and m'
# 4. sigma2(X)
# 5. E(X / sigma2(X) | beta'X) and E(1 / sigma2(X) | beta'X)

# Now test the function
score <- compute_score_eq7(
            Y = Y,
            X = X,
            C = C,
            EY_given_XC = EY_XC,
            pX = pX,
            pXC = pXC,
            m_hat = m_hat,
            mprime_hat = mprime_hat,
            sigma2_hat = sigma2_hat,
            EX_over_sigma2 = EX_over_sigma2,
            E_inv_sigma2 = E_inv_sigma2
)


# newton Raphson step
MASS::ginv(score$S2_bar) %*% score$S_bar_vec

#debugonce(compute_score_eq7)
#debugonce(fit_pd_list)
