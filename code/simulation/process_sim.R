library(data.table)
library(ggplot2)
library(here)

# Hyperparameters ---------------------------------------------------------
corr <- "high"
nvar <- 20

# Process -----------------------------------------------------------------
dir      <- paste0("enar-student-paper2_", corr, "-rho", "_p", nvar)
sim_path <- here("outputs", "simulation", dir)
dt       <- lapply(list.files(sim_path, full.names = TRUE), fread) |>
                   rbindlist(fill = TRUE)

dt[, method := ifelse(method == "CS", "csPCA (L=1)",
                      ifelse(method == "CS-CF", "csPCA (L=5)", method))]

dt.summary <- dt[, .(mean = mean(value), sd = sd(value),
       time_avg = mean(time), time_sd = sd(time)), .(n, error_type, method)]

# Boxplot of errors -------------------------------------------------------
ggplot(dt, aes(method, value)) + 
  geom_boxplot(aes(col = method)) + 
  facet_grid(error_type~n, scales = "free") + 
  theme_bw()
ggsave("results/simulation-boxplot.png", width = 11, height = 8.5)

# Table of Errors ---------------------------------------------------------
dt.summary[, value := paste0(round(mean, 4), " (", round(sd, 4), ")")]
dt.wide <- dt.summary[error_type == "2-norm"] |> 
  dcast(n ~ method, value.var = "value")
dt.wide <- dt.summary[error_type == "Frobenius"] |> 
  dcast(n ~ method, value.var = "value")

dt.wide[, .(n, PCA, pCCA, SDR, `csPCA (L=1)`, `csPCA (L=5)`)]

# Sqrt(n) Error -----------------------------------------------------------
ggplot(dt.summary[method != "Truth"], aes(n^{-1/2}, mean)) + 
  geom_point() + 
  geom_line(aes(group = method, col = method)) +
  facet_wrap(~error_type, scales = "free") + 
  xlab("Sample size (n)") + 
  ylab("Error") + 
  theme_bw() + ggtitle("Error") + 
  #ylim(c(0, 2))  + 
  xlim(c(0, .1)) + 
  #geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)
ggsave("results/simulation-avg-error.png")

# MSE of mu(X) ------------------------------------------------------------
ggplot(a, aes(n^{-1/2}, V1)) + 
  geom_point() + 
  #ylim(c(0,1)) + 
  geom_hline(yintercept = 0) +
  xlim(c(0, .1)) +
  geom_abline(slope=1, intercept=0)

# MSE of mu(Z) ------------------------------------------------------------

dt[, .(times = mean(dhat == 2)), .(n, method)]
dt[, .(mean(mu_Z_mse, na.rm = TRUE), sd(mu_Z_mse, na.rm = TRUE)), .(n, method)]
a <- dt[method %in% c("CS", "CS-CF", "SDR", "Truth"), 
   .(mse=mean(mu_X_mse, na.rm = TRUE), sd=sd(mu_X_mse, na.rm = TRUE)), .(n, method)]

a <- dt[method %in% c("csPCA (L=1)", "csPCA (L=5)", "SDR", "Truth"), 
   .(mse=mean(mu_Z_mse, na.rm = TRUE), sd=sd(mu_Z_mse, na.rm = TRUE)), .(n, method)]

a[, value := paste0(round(mse, 4), " (", round(sd, 4), ")")]
a.wide <- a |> 
  dcast(n ~ method, value.var = "value")
a.wide[, .(n, Truth, SDR, `csPCA (L=1)`, `csPCA (L=5)`)]

a[, lower := mse - 1.96*sd]
a[, upper := mse + 1.96*sd]

ggplot(a, aes(n, mse)) + 
  geom_point()+ 
  #geom_errorbar(aes(ymin = lower, ymax = upper)) + 
  geom_line(aes(group = method, col = method)) + 
  #ylim(c(0,0.3)) + 
  geom_hline(yintercept = 0) + 
  theme_bw() + 
  xlim(c(0,1600)) + 
  geom_vline(xintercept = 0) + 
  ylab("MSE of mu(Z)") + 
  xlab("Sample size (n)")  
  #facet_wrap(~method)

