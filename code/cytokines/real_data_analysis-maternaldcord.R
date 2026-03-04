library(data.table)
library(here)
library(MAVE)
library(SuperLearner)
library(tmle)
library(magrittr)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggplot2)
library(stringr)

# Source functions --------------------------------------------------------
source(here("R/csMAVE.R"))
source(here("R/gcomp.R"))
source(here("R/estimate_Seff.R"))
#source(here("R/estimate_DR_curve.R"))
source(here("R/plot_factor_loadings.R"))

#--- csPCA estimation ----------------------------------------------------------------
set.seed(123)
DT     <- fread(here("data/cytokines/processed/analysis_maternal-to-cord.csv"),
                stringsAsFactors = TRUE, na.strings = c("", NA))
y_name <- "IL-6_c" 
x_cols <- 2:10
c_cols <- 24:32

x_cols <- which( DT[, lapply(.SD, function(x) uniqueN(x) / 100), .SDcols = x_cols][1,] |>
  as.numeric() > .5) + 1

Y <- DT[[y_name]]
X <- DT[, x_cols, with = FALSE]
C <- DT[, c_cols, with = FALSE] 

SL.library <- c("SL.glm", 
                "SL.glmnet", 
                "SL.earth",
                "SL.ranger",
                "tmle.SL.dbarts2")

fit <- csMAVE(
  Y = Y,
  X = X,
  C = C, 
  mu_args   = list(SL.library = SL.library, method = "method.NNLS2")
)

cs.dim <- mave.dim(fit$mave)
X <- as.matrix(DT[, x_cols, with = FALSE])
# No column names allowed for Rcpp
colnames(X) <- NULL
beta <- run_efficient_estimator(X = X, Y = fit$mu_X, beta_init = fit$mave$dir[[1]])


fit.mave <- MAVE::mave(Y ~ ., data = data.frame(Y, X), method = "meanMAVE")

mave.dim <- mave.dim(fit.mave)

# Plot of csPCA and MAVE loadings -----------------------------------------------------------------
cv_table <- data.frame(d=1:10, cs.dim$cv)
plt1 <- plot_factor_loadings(fit$mave, 1) + ggtitle("d=1") + theme(plot.title = element_text(size = 10))
plt2 <- plot_factor_loadings(fit$mave, 2) + ggtitle("d=2") + theme(plot.title = element_text(size = 10))

plt <- plt1 / plt2 
# ggsave(here("results/figures/csPCA Loadings for cytokines (Maternal to Cord Blood IL-6.png"), 
#        width = 11, height=4, units = "in")
# ggsave(here("results/figures/csPCA Loadings for cytokines (Maternal to Cord Blood IL-6.pdf"), 
#        width = 11, height=4, units = "in")
ggsave(here("results/figures/csPCA_IL6.pdf"), 
       width = 11, height=4, units = "in")
ggsave(here("results/figures/csPCA_IL6.png"), 
       width = 11, height=4, units = "in")

plt1.mave <- plot_factor_loadings(fit.mave, 1) 
plt2.mave <- plot_factor_loadings(fit.mave, 2)
plt3.mave <- plot_factor_loadings(fit.mave, 3)

plt.mave <- plt1.mave / plt2.mave / plt3.mave
ggsave(here("results/figures/MAVE Loadings for cytokines (Maternal to Cord Blood IL-6.png"),
       width = 11, height=7, units = "in")
ggsave(here("results/figures/MAVE Loadings for cytokines (Maternal to Cord Blood IL-6.pdf"), 
       width = 11, height=7, units = "in")
ggsave(here("results/figures/MAVE_IL6.png"), 
       width = 11, height=7, units = "in")

# Scatterplot of Z against mu(X) --------------------------------------
Z       <- as.matrix(X) %*% fit$mave$dir[[1]]
# we want TNF-alpha to align with positive. Try to make this more robust.
Z <- -1*Z

newdf <- cbind(muhat_X = fit$mu_X, 
               Y, Z, X)

plt.scatter <- ggplot(newdf, aes(dir1, muhat_X)) + 
                geom_point() + 
                theme_bw() + 
  xlab(expression(Z[1])) + 
  ylab(expression(mu[x]))

# ggsave(filename = here("results/figures/Real Data Analysis Scatter of Z1 against mu(X).png"), 
#        units = "in", width = 6, height = 4)
ggsave(filename = here("results/figures/Real Data Analysis Scatter of Z1 against mu(X).pdf"), 
       units = "in", width = 6, height = 4)
ggsave(filename = here("results/figures/RDA_Z1_muX.pdf"), 
       units = "in", width = 6, height = 4)
ggsave(filename = here("results/figures/RDA_Z1_muX.png"), 
       units = "in", width = 6, height = 4)

# Estimate DR causal curve with Z -----------------------------------------
# set.seed(123)
# muhat.Z <- gcomp(Y, Z, C, SL.library = SL.library)
# muhat.Z <- estimate_DR_curve(Y, Z, C, Z.new = 0)
# newdf <- cbind(muhat_Z = muhat.Z, 
#                muhat_X = fit$mu_X, 
#                Y, Z, X)
# 
# ggplot(newdf, aes(dir1, dir2)) + geom_point()
# 
