set.seed(1)

library(MAVE)
library(SuperLearner)

source("R/gcomp.R")

## ---------------------------
## Simulation setup
## ---------------------------

n  <- 1000
p  <- 3
q  <- 3
beta0 <- c(0, 0, 1)
B0 <- matrix(beta0 / sqrt(sum(beta0^2)), p, 1)   # normalized true direction

# Confounders
C <- matrix(rnorm(n * q), n, q)

# Exposure model: X = m(C) + U
A <- matrix(c(1, -0.5, 0.3,
              0.5,  0.2, 0.4,
              -0.3,  0.6, 0.1), q, p, byrow = TRUE)
mC <- C %*% A
U  <- matrix(rnorm(n * p), n, p)
X  <- mC + U

# Outcome model: Y = g(beta0'X) + h(C) + error
g <- function(t) sin(t)
h <- function(c) 2 * c[,1] - c[,2] + 0.5 * c[,3]
eps <- rnorm(n)
Y <- g(X %*% beta0) + h(C) + eps

muX <- g(X %*% beta0)
muX_hat <- gcomp(Y, X, C, SL.library = c("SL.glm", "SL.gam", "SL.glmnet", "SL.xgboost", "SL.randomForest"))

# Residualization
# Ytilde <- residuals(lm(Y ~ C))
Ytilde <- residuals(loess(Y ~ C))
# Xtilde <- apply(X, 2, function(x) residuals(lm(x ~ C)))
Xtilde <- apply(X, 2, function(x) residuals(loess(x ~ C)))

## ---------------------------
## Helper functions
## ---------------------------

normalize_dir <- function(b) b / sqrt(sum(b^2))

proj_mat <- function(B) B %*% solve(t(B) %*% B) %*% t(B)

frobenius_norm <- function(P1, P2) sqrt(sum((P1 - P2)^2))
spectral_norm  <- function(P1, P2) max(svd(P1 - P2)$d)

compare_proj <- function(Bhat, B0) {
  P1 <- proj_mat(Bhat)
  P0 <- proj_mat(B0)
  list(
    fro = frobenius_norm(P1, P0),
    spec = spectral_norm(P1, P0)
  )
}

## ---------------------------
## Fit MAVE
## ---------------------------

fit_raw   <- mave(Y ~ X, method = "meanMAVE")
fit_mu    <- mave(muX ~ X, method = "meanMAVE")
fit_muhat    <- mave(muX_hat ~ X, method = "meanMAVE")
fit_resid <- mave(Ytilde ~ Xtilde, method = "meanMAVE")

# Extract directions (standardized)
B_raw   <- normalize_dir(fit_raw$dir[[1]])
B_mu    <- normalize_dir(fit_mu$dir[[1]])
B_muhat    <- normalize_dir(fit_muhat$dir[[1]])
B_resid <- normalize_dir(fit_resid$dir[[1]])

## ---------------------------
## Compare with true span(beta0)
## ---------------------------

errs <- rbind(
  observed   = unlist(compare_proj(B_raw,   B0)),
  causal_truth    = unlist(compare_proj(B_mu,    B0)),
  causal_est    = unlist(compare_proj(B_muhat,    B0)),
  residual = unlist(compare_proj(B_resid, B0))
)

print(round(errs, 4))
