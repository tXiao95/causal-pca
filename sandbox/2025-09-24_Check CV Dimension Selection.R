library(data.table)
library(ggplot2)
library(here)
library(MAVE)
library(parallel)
library(patchwork)

source(here("R/simulate_data.R"))

# Intro -------------------------------------------------------------------

# In a previous simulation study, I found that the cross-validation dimension selection
# method in MAVE did a poor job of selecting the right dimension. I want to now check it that is true.

# Previous simulation setting ---------------------------------------------
set.seed(123)
B <- 100

results.old <- mclapply(1:B, function(i){
  data <- simulate_data(500, p=20, q=10, rho=0.8, h_Z_coef = 0.1, 
                        Z1_coef = 1, Z2_coef = 1, Z12_coef = 1)
  mave <- MAVE::mave(Y ~ ., data = data.frame(Y = data$mu_X, data$X), method = "meanMAVE")
  return(mave.dim(mave)$dim.min)
}, mc.cores = 4)

mean(results.old==2)
# 0.0

# Different Z -------------------------------------------------------------
set.seed(123)

results.highZ1 <- mclapply(1:B, function(i){
  data <- simulate_data(500, p=20, q=10, rho=0.8, h_Z_coef = 0.1, 
                        Z1_coef = 10, Z2_coef = 1, Z12_coef = 1)
  mave <- MAVE::mave(Y ~ ., data = data.frame(Y = data$mu_X, data$X), method = "meanMAVE")
  return(mave.dim(mave)$dim.min)
}, mc.cores = 4)

mean(results.highZ1==2)
# [1] 1

# Testing a new Z for simulation-------------------------------------------------------------
set.seed(123)
results.new <- mclapply(1:B, function(i){
  data <- simulate_data(500, p=20, q=10, rho=0.8, h_Z_coef = 0.01, 
                        Z1_coef = 15, Z2_coef = 5, Z12_coef = 5)
  mave <- MAVE::mave(Y ~ ., data = data.frame(Y = data$mu_X, data$X), method = "meanMAVE")
  return(mave.dim(mave)$dim.min)
}, mc.cores = 4)

mean(results.new==2)
# [1] 1

# Visualize difference in Z -------------------------------------------------------------
set.seed(123)
data.old    <- simulate_data(500, p=20, q=10, rho=0.8, h_Z_coef = 0.1,
                      Z12_coef = 1, Z1_coef = 1, Z2_coef = 1)
data.highZ1 <- simulate_data(500, p=20, q=10, rho=0.8, h_Z_coef = 0.01,
                      Z12_coef = 1, Z1_coef = 10, Z2_coef = 1)
data.new <- simulate_data(500, p=20, q=10, rho=0.8, h_Z_coef = 0.01,
                      Z12_coef = 5, Z1_coef = 15, Z2_coef = 5)

# Z variables
Z.old    <- data.old$X %*% data.old$beta
Z.highZ1 <- data.highZ1$X %*% data.highZ1$beta
Z.new    <- data.new$X %*% data.new$beta

n <- nrow(Z.new)

newdf <- data.table(rbind(Z.old, Z.highZ1, Z.new), 
                    type = rep(c("Old", "High Z1", "New"), each = n), 
                    mu_X = c(data.old$mu_X, data.highZ1$mu_X, data.new$mu_X))

plt1 <- ggplot(newdf, aes(V1, mu_X)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~type) + 
  geom_smooth() + 
  xlab("Z1") + 
  ggtitle("Z1 vs. mu(X)")

plt2 <- ggplot(newdf, aes(V2, mu_X)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~type) + 
  geom_smooth() + 
  xlab("Z2") + 
  ggtitle("Z2 vs. mu(X)")

plt <- plt1 / plt2

plt
ggsave(here("sandbox/Z1 and Z2 difference for dimension reduction selection.png"), 
       units = "in", height = 6, width = 11)