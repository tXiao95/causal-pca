library(data.table)
library(here)

args <- commandArgs(trailingOnly = TRUE)
N <- as.numeric(args[1])

# Setup -------------------------------------------------------------------

# datpath <- here("data", "birthweight", "ATL20County_Linked.RData")
# load(datpath)
# 
# dat <- data.table(dat)
# Y_var <- "grams"
# old_X_vars <- names(dat)[grepl("Total.", names(dat))] |>
#         setdiff(c("Total.Temp", "Total.dewp")) 
# 
# C_vars <- c(
#         "m_race", 
#         "m_age", 
#         "edu", 
#         "married", 
#         "tobacco", 
#         "alcohol", 
#         "prev", 
#         "gestweeks", 
#         "county", 
#         "BG_POV_LVL"
# )
# 
# dat_subset <- dat[, c(old_X_vars, C_vars, Y_var), with = FALSE]
# 
# X_vars <- old_X_vars |> stringr::str_replace("Total.", "")
# 
# setnames(dat_subset, old_X_vars, X_vars)
# 
# yesno_vars  <- c("tobacco", "alcohol", "prev")
# factor_only <- "gestweeks"
# 
# # Recode to Yes/No factors
# dat_subset[, (yesno_vars) := lapply(.SD, function(x) {
#         factor(ifelse(x > 0, "Yes", "No"),
#                levels = c("No", "Yes"))
# }), .SDcols = yesno_vars]
# 
# # Convert gestweeks to factor (no relabeling)
# dat_subset[, (factor_only) := lapply(.SD, factor),
#     .SDcols = factor_only]
# 
# fwrite(dat_subset, here( "data", "birthweight", "ATL_bw_data.csv"))

# Analysis ----------------------------------------------------------------

dat_subset <- fread(here("data", "birthweight", "ATL_bw_data.csv"), 
                    stringsAsFactors = TRUE)

X_vars <- c("CO", "EC", "NH4", "NO2", "NO3", "NOx", 
            "OC", "O3", "PM10", "PM25", "SO2", "SO4")

C_vars <- c(
        "m_race",
        "m_age",
        "edu",
        "married",
        "tobacco",
        "alcohol",
        "prev",
        "gestweeks",
        "county",
        "BG_POV_LVL"
)

X <- dat_subset[, ..X_vars]
C <- dat_subset[, ..C_vars]
Y <- dat_subset$grams

n <- nrow(dat_subset)
p <- ncol(X)
q <- ncol(C)

source(here("R", "csMAVE.R"))
source(here("R", "estimate_Seff.R"))
source(here("R", "outcome_regression.R"))

set.seed(12)
sampled_idx <- sample(1:n, N, replace = FALSE)
dat_sampled <- dat_subset[sampled_idx,]

SL.lib <- c( 
            "SL.glmnet",
            "SL.speedglm",
            "SL.xgboost", 
            "SL.earth")

# 7 minutes for training
start_time <- proc.time()
out_model_X <- outcome_model(Y = dat_sampled$grams, 
                             X = dat_sampled[, ..X_vars], 
                             C = dat_sampled[, ..C_vars], 
                             mu_fitter = SL_outcome_fitter, 
                             SL.lib = SL.lib, 
                             cvControl = list(V = ifelse(N >= 10000, 2, 5)))
out_model_X_time <- proc.time() - start_time

# I think this gets slower the bigger the model is
start_time <- proc.time()
mu_X <- gcomp(out_model_X, C= dat_sampled[, ..C_vars], X_new = dat_sampled[, ..X_vars])
gcomp_time <- proc.time() - start_time

dat_sampled$mu_X <- mu_X

form     <- reformulate(X_vars, response = "mu_X")
start_time <- proc.time()
mave_fit <- MAVE::mave(form, data = dat_sampled, method = "meanOPG", max.dim = p)
mave_time <- proc.time() - start_time

start_time <- proc.time()
dhat     <- MAVE::mave.dim(mave_fit)
dhat_time <- proc.time() - start_time

beta <- mave_fit$dir[[dhat$dim.min]] |>
    qr() |>
    qr.Q()

d <- ncol(beta)

rownames(beta) <- X_vars

P <- Pi(beta)
lambdas <- sort((diag(P)))

# Transform the X vars
Z <- as.matrix(dat_sampled[, ..X_vars]) %*% beta
Z_vars <- paste0("Z", 1:d)

dat_sampled[, (Z_vars) := as.data.table(Z)]

start_time <- proc.time()
out_model_Z <- outcome_model(Y = dat_sampled$grams, 
                             X = dat_sampled[, ..Z_vars], 
                             C = dat_sampled[, ..C_vars], 
                             mu_fitter = SL_outcome_fitter, 
                             SL.lib = SL.lib, 
                             cvControl = list(V = ifelse(N >= 10000, 2, 5)))
out_model_Z_time <- proc.time() - start_time

obj <- list(dat_sampled = dat_sampled, # Dataset
            out_model_X = out_model_X, # Outcome regression 
            mu_X = mu_X, # pseudo-outcomes 
            mave_fit = mave_fit, # MAVE fit for 1:p
            dhat = dhat, # CV dimension selection object
            beta = beta, # Selected beta
            P = P, # projection matrix
            lambdas = lambdas, # sorted eigenvalues of P
            out_model_Z = out_model_Z, 
            out_model_X_time = out_model_X_time,
            gcomp_time = gcomp_time,
            mave_time = mave_time,
            dhat_time = dhat_time,
            out_model_Z_time = out_model_Z_time
            )

saveRDS(obj, file = here("results", paste0("bw_obj_N=", N, ".rds")))
