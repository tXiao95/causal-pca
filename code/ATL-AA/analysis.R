library(data.table)
library(ggplot2)
library(here)
library(readxl)
library(qgcomp)
library(bkmr)
library(pcpr)
library(SuperLearner)
library(MAVE)
library(torch) # Required for the Neural Network

# Ensure you have all wrappers sourced
source(here("R/nuisance_outcome_regression.R")) 
source(here("R/nuisance_gps.R")) # Contains mvn_fitter
source(here("R/crossfit_ERS.R")) # Contains crossfit_ERS
source(here("R/estimate_ERS.R"))

X_path <- here("data/ATL-AA/ATL_AA_PFAS_N=532_two_visits.csv")
C_path <- here("data/ATL-AA/MPTB_AA Cohort Data_N766, clean, DEidentified, data dictionary_12.5.2025.xlsx")

# Data path 
X_dat <- fread(X_path)[Visit == 1, .(subjectid, PFHXS, PFOS, PFOA, PFNA, PFDA, PFUNDA, PFDODA)]
setnames(X_dat, "subjectid", "Subjectid")

C_data <- read_excel(C_path, sheet = 1) |> data.table()
dict   <- read_excel(C_path, sheet = 2) |> data.table()

# Pick covariates
C_dat <- C_data[AllFullTerm == 1, .(Subjectid, 
                    age_enrollment, 
                    Education_4.level,
                #    income_5cat, 
                   # Parity_3cat,
                    FirstPrenatalBMI,
                  #  MarriedCohab_Not,
                    TobaccoUse_MRorSR,
                    #AlcoholUse_MRorSR,
                    MarijuanaUse_MRorSR,
                    Sex,
                    #birthga,
                    birth_weight)]

C_dat[, birth_weight := ifelse(birth_weight == "NA", NA, as.numeric(birth_weight) )]

dat <- merge(X_dat, C_dat, by = "Subjectid")


# Apply log2 transformation to all PFAS concentrations to reduce influenc eof outliers
dat[, `:=`(PFHXS = log2(PFHXS), 
           PFOS = log2(PFOS),
           PFOA = log2(PFOA),
           PFNA = log2(PFNA))]
           #PFDA = log2(PFDA),
           #PFUNDA = log2(PFUNDA),
           #PFDODA = log2(PFDODA))]

# 
dat[, `:=`(edu = factor(Education_4.level, levels = 1:4, labels = c("Less than HS", 
                                                                    "HS or GED", 
                                                                    "Some college or tech school",
                                                                    "4-yr college or more")),
           bmi = FirstPrenatalBMI,
           age = age_enrollment,
           tobacco = factor(TobaccoUse_MRorSR, levels = 0:1, labels = c("No", "Yes")),
           marijuana = factor(MarijuanaUse_MRorSR, levels = 0:1, labels = c("No", "Yes")),
           sex = as.factor(Sex))]

dat <- dat[, .(PFOS, PFOA, PFNA, PFHXS, 
               #PFDA, PFUNDA, PFDODA, 
               age, bmi, edu, tobacco, marijuana, sex, birth_weight)]

# Setup & Libraries -------------------------------------------------------

X_vars <- c("PFOS", "PFOA", "PFNA", "PFHXS")
C_vars <- c("age", "bmi", "edu", "tobacco", "marijuana", "sex")
Y_var  <- "birth_weight"

p <- length(X_vars)
n <- nrow(dat)

# Other methods -----------------------------------------------------------
# 2. Matrix Preparation for BKMR and PCPR ---------------------------------

# Convert to standard data.frame for safer formula/matrix parsing
dat_df <- as.data.frame(dat)

# Exposures matrix (Z)
Z_mat <- as.matrix(dat_df[, X_vars])

# Outcome vector (Y)
Y_vec <- dat_df[[Y_var]]

# Covariates matrix (X_mat)
# model.matrix creates the necessary dummy variables for your factors (edu, tobacco, etc.)
# We use [,-1] to drop the intercept column since the models fit their own intercept
formula_C <- as.formula(paste("~", paste(C_vars, collapse = " + ")))
C_mat <- model.matrix(formula_C, data = dat_df)[, -1]


# 3. Quantile G-Computation (qgcomp) --------------------------------------
cat("\nRunning qgcomp...\n")

formula_full <- as.formula(
  paste("birth_weight ~", paste(c(X_vars, C_vars), collapse = " + "))
)

set.seed(123)
qgcomp_fit <- qgcomp.noboot(
  f = formula_full,
  expnms = X_vars,
  data = dat_df,
  family = gaussian(),
  q = 4 # Quartiles
)

print(summary(qgcomp_fit))
plot(qgcomp_fit)


set.seed(123)
qgcomp_curve_fit <- qgcomp.boot(
  f = formula_full,
  expnms = X_vars,
  data = dat_df,
  family = gaussian(),
  q = 4,
  B = 500,
  degree = 3   # degree = 2 allows for a quadratic (curved) dose-response shape
)

# Plot the curved dose-response line
plot(qgcomp_curve_fit)

ggsave(here("results/jasa-initial-submission/ATL-AA/qgcomp.pdf"),
  width = 7.2, 
  height = 3.2, 
  dpi = 300
)

# 4. Bayesian Kernel Machine Regression (bkmr) ----------------------------
cat("\nRunning bkmr (this may take a while)...\n")

set.seed(123)
bkmr_fit <- kmbayes(
  y = Y_vec,
  Z = Z_mat,
  X = C_mat,
  iter = 10000, 
  family = "gaussian",
  varsel = TRUE # Enables calculation of Posterior Inclusion Probabilities (PIPs)
)

# BKMR Diagnostics and Outputs
TracePlot(fit = bkmr_fit, par = "beta")
TracePlot(fit = bkmr_fit, par = "r")
TracePlot(fit = bkmr_fit, par = "sigsq.eps")

pips <- ExtractPIPs(bkmr_fit)
print(pips)
# Example plot: Univariate predictor-response functions
# Your existing base plot
pred_resp <- PredictorResponseUnivar(fit = bkmr_fit)

ggplot(pred_resp, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable, scales = "free_x") + # 'free_x' allows each chemical to have its own x-axis range
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) + 
  theme_bw() +
  # Add your clean, descriptive labels here
  labs(
    x = "Log2-Transformed PFAS Concentration",
    y = "Estimated Change in Birth Weight (g)"
  )

ggsave(here("results/jasa-initial-submission/ATL-AA/bkmr_univariate.pdf"),
       width = 7.2, 
       height = 4.2, 
       dpi = 300
)

# Code to check overall mixture effect
risks.overall <- OverallRiskSummaries(fit = bkmr_fit, 
                                      y = Y_vec, 
                                      Z = Z_mat, 
                                      X = C_mat, 
                                      qs = seq(0, 1, by = 0.005), 
                                      q.fixed = 0.5)
risks.overall

ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange() + 
  theme_bw()

ggplot(risks.overall, aes(x = quantile, y = est)) + 
  # Add a reference line at zero for easy interpretation
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  
  # geom_ribbon creates the shaded uncertainty band. 
  # alpha = 0.2 makes it transparent so you can see the grid behind it.
  # Note: mapping ymin and ymax must be inside an aes() call.
  geom_ribbon(aes(ymin = est - 1.96*sd, ymax = est + 1.96*sd), fill = "black", alpha = 0.2) + 
  
  # geom_line draws the main estimate curve
  geom_line(size = 1, color = "black") + 
  
  theme_bw() +
  labs(
    x = "Joint PFAS Quantile",
    y = "Expected Difference in Birth Weight (vs. Median)"
  )

ggsave(here("results/jasa-initial-submission/ATL-AA/bkmr_diff.pdf"),
       width = 7.2, 
       height = 4.2, 
       dpi = 300
)

# CSDR --------------------------------------------------------------------


# 2. Fit the Global Outcome Model E[Y | X, C]
# We use SuperLearner to flexibly model the outcome surface
SL.lib <- c("SL.glmnet", "SL.glm", "SL.earth", "SL.xgboost", "SL.ranger")

message("Fitting global outcome regression...")
set.seed(123)
out_model_X <- outcome_model(
  Y = dat[[Y_var]], 
  X = dat[, ..X_vars], 
  C = dat[, ..C_vars], 
  mu_fitter = SL_outcome_fitter, 
  SL.lib = SL.lib,
  cvControl = list(V = 10) 
)

# 3. Estimate Pseudo-Outcomes
# Evaluate the outcome model to extract the pseudo-outcomes (mu_X) using Regression Adjustment
message("Estimating pseudo-outcomes...")
mu_X_obj <- estimate_ERS(
  Y = dat[[Y_var]],
  X = dat[, ..X_vars],
  C = dat[, ..C_vars],
  x_eval = dat[, ..X_vars],
  out_model = out_model_X,
  estimator = "RA",
  return_vector = TRUE # Returns just the vector of predictions
)

dat$mu_X <- mu_X_obj

# 4. Dimension Reduction via MAVE
message("Performing dimension reduction...")
form <- reformulate(X_vars, response = "mu_X")

# Fit MAVE and use cross-validation to select the optimal dimension
mave_fit <- MAVE::mave(form, data = dat, method = "meanMAVE", max.dim = p)
dhat <- MAVE::mave.dim(mave_fit)


# The cross-validation is run on dimensions of 0 1 2 3 4 
# Dimension	0 	1 	2 	3 	4 	
# CV-value	3296.32 	5.25 	5.58 	6.3 	7.45 

# Extract the orthonormal projection matrix (beta) for the selected dimension
beta <- mave_fit$dir[[dhat$dim.min]] |>
  qr() |>
  qr.Q()

d <- ncol(beta)
rownames(beta) <- X_vars

# 5. Transform Exposures into Low-Dimensional Index (Z)
message(sprintf("Selected optimal dimension: d = %d", d))
Z <- as.matrix(dat[, ..X_vars]) %*% beta
Z_vars <- paste0("Z", 1:d)
dat[, (Z_vars) := as.data.table(Z)]

# 6. Final Causal Evaluation on Z using Cross-Fitting and DR
message("Estimating final causal dose-response curve with cross-fitting (L=5)...")

if (d == 1) {
  # Standard 1D grid
  z_grid_vals <- seq(quantile(dat$Z1, 0.001), quantile(dat$Z1, 0.95), length.out = 70)
  z_eval_df <- data.frame(Z1 = z_grid_vals)
  
} else if (d == 2) {
  # 2D Mesh Grid for a Surface/Contour Plot
  z1_vals <- seq(quantile(dat$Z_assoc1, 0.05), quantile(dat$Z_assoc1, 0.95), length.out = 30)
  z2_vals <- seq(quantile(dat$Z_assoc2, 0.05), quantile(dat$Z_assoc2, 0.95), length.out = 30)
  
  # expand.grid creates every possible combination of Z1 and Z2 (900 rows)
  z_assoc_eval_df <- expand.grid(Z_assoc1 = z1_vals, Z_assoc2 = z2_vals)
}

# Estimate the final curve over the grid using 5-fold cross-fitting
# This automatically trains NN_outcome_fitter and mvn_fitter inside each fold
final_ERS <- crossfit_ERS(
  Y = dat[[Y_var]],
  X = dat[, ..Z_vars],
  C = dat[, ..C_vars],
  x_eval = z_eval_df,
  estimator = "DR",
  L = 5,
  outcome_fitter = SL_outcome_fitter,  # Neural net for smooth outcome estimation
  gps_fitter = mvn_fitter,             # Multivariate normal GPS
  optimize_bw = TRUE,                  # Dynamically calculates AMSE-optimal bandwidth
  seed = 42                            # For reproducible CV folds
)

ers_plt <- ggplot(final_ERS$results, aes(x = -Z1, y = estimate)) + 
  # 1. Subtle density marks (rug plot) using the observed data
  #geom_point() + 
  geom_rug(data = dat, aes(x = Z1), inherit.aes = FALSE, 
           alpha = 0.2, sides = "b", length = unit(0.05, "npc")) +
  
  # 2. Main curve and CI ribbon
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, fill = "steelblue") + 
  geom_line(linewidth = 1.2, color = "darkblue") + 
  
  # 3. Slide-ready theme and text sizes
  theme_bw(base_size = 11) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  ) +
  
  # 4. Clean mathematical x-axis label
  xlab(expression(Z == 0.77*PFOS + 0.14*PFOA + 0.26*PFNA - 0.56*PFHxS )) + 
  ylab("Birthweight (g)") + 
  xlim(c(-1.8, 2.3)) +  
  ggtitle(NULL)

ers_plt
ggsave(filename = "results/jasa-initial-submission/ATL-AA/causal_ers_no_title.pdf",
       plot = ers_plt, width = 7.2, height = 4)


# MAVE --------------------------------------------------------------------

message("Performing Association Dimension Reduction (Standard MAVE)...")

# 1. Direct formula: Regress Y directly on X (skip pseudo-outcomes)
form_assoc <- reformulate(X_vars, response = Y_var)

# 2. Fit standard MAVE and select dimension
mave_fit_assoc <- MAVE::mave(form_assoc, data = dat, method = "meanMAVE", max.dim = p)
dhat_assoc <- MAVE::mave.dim(mave_fit_assoc)

# 3. Extract the naive association projection matrix
#beta_assoc <- mave_fit_assoc$dir[[dhat_assoc$dim.min]] |>
beta_assoc <- mave_fit_assoc$dir[[dhat$dim.min]] |>
  qr() |>
  qr.Q()

d_assoc <- ncol(beta_assoc)
rownames(beta_assoc) <- X_vars

# 4. Transform Exposures into Association Index (Z_assoc)
message(sprintf("Selected optimal association dimension: d = %d", d_assoc))
Z_assoc_mat <- as.matrix(dat[, ..X_vars]) %*% beta_assoc
Z_assoc_vars <- paste0("Z_assoc", 1:d_assoc)
dat[, (Z_assoc_vars) := as.data.table(Z_assoc_mat)]

# Association
z_assoc_grid_vals <- seq(quantile(dat$Z_assoc1, 0.01), quantile(dat$Z_assoc1, 0.95), length.out = 100)
z_assoc_eval_df <- data.frame(Z_assoc1 = z_assoc_grid_vals)

# Handle multi-dimensional grids if MAVE selected d > 1
if (d_assoc > 1) {
  for (dim_idx in 2:d_assoc) {
    z_assoc_eval_df[[paste0("Z_assoc", dim_idx)]] <- median(dat[[paste0("Z_assoc", dim_idx)]])
  }
}

# Cross-fit the curve along the Association dimension
assoc_ERS <- crossfit_ERS(
  Y = dat[[Y_var]],
  X = dat[, ..Z_assoc_vars],
  C = dat[, ..C_vars],
  x_eval = z_assoc_eval_df,
  estimator = "DR",
  L = 5,
  outcome_fitter = SL_outcome_fitter,
  gps_fitter = mvn_fitter,
  optimize_bw = TRUE,
  seed = 42
)

assoc_plt <- ggplot(assoc_ERS$results, aes(x = Z_assoc1, y = estimate)) + 
  # 1. Subtle density marks (rug plot) using the observed data
  geom_rug(data = dat, aes(x = Z_assoc1), inherit.aes = FALSE, 
           alpha = 0.2, sides = "b", length = unit(0.05, "npc")) +
  
  # 2. Main curve and CI ribbon
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, fill = "steelblue") + 
  geom_line(linewidth = 1.2, color = "darkblue") + 
  
  # 3. Slide-ready theme and text sizes
  theme_bw(base_size = 18) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  ) +
  
  # 4. Clean mathematical x-axis label
  xlab(expression(Z == 0.27*PFHxS - 0.41*PFOS - 0.83*PFOA - 0.27*PFNA)) + 
  ylab("Birthweight (g)") + 
  xlim(c(-1.95, 2.1)) +  
  #ylim(c(3200, 3651)) + 
  ggtitle("Association-based Exposure Response Surface")
assoc_plt

# Causal ERS
ggsave(filename = "results/jasa-initial-submission/ATL-AA/assocation_ers.pdf",
       plot = assoc_plt, width = 8, height = 6)