library(ggplot2)
library(here)
library(dplyr)

plot_zero_noise_variance <- function(n = 300, p = 10) {
  # 1. Generate STRICTLY ZERO NOISE data
  X <- matrix(runif(n * p, -2, 2), n, p)
  beta_true <- matrix(0, p, 1)
  beta_true[1:4, 1] <- 0.5
  
  # Deterministic non-linear relationship (no epsilon!)
  betaX <- as.vector(X %*% beta_true)
  Y <- sin(betaX) 
  
  # 2. Get MAVE estimates for the projected space
  fit_mave <- MAVE::mave(Y ~ X, method = "meanMAVE")
  beta_est <- fit_mave$dir[[1]]
  betaX_est <- as.vector(X %*% beta_est)
  
  # 3. Estimate the conditional variance using a simple Nadaraya-Watson / Local kernel
  # (Simulating what your loop does internally)
  b <- n^(-1/5) # standard bandwidth
  dist_mat <- as.matrix(dist(betaX_est))
  W <- dnorm(dist_mat / b)
  W_norm <- sweep(W, 1, rowSums(W), "/")
  
  # m_hat is the kernel smoothed estimate of Y
  m_hat <- as.vector(W_norm %*% Y)
  
  # sigma2_hat is the estimated variance of the residuals
  residuals <- Y - m_hat
  sigma2_hat <- as.vector(W_norm %*% (residuals^2))
  
  # 4. Plot the results
  df <- data.frame(
    Index = betaX_est,
    True_Y = Y,
    Estimated_m = m_hat,
    Estimated_Variance = sigma2_hat
  )
  
  p1 <- ggplot(df, aes(x = Index)) +
    geom_point(aes(y = True_Y, color = "True Deterministic Y"), alpha = 1.5) +
    geom_line(aes(y = Estimated_m, color = "Kernel Estimate m(x)"), size = 1) +
    theme_minimal() +
    labs(title = "The Deterministic Fit (epsilon = 0)", y = "Y") +
    scale_color_manual(values = c("black", "red")) +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  p2 <- ggplot(df, aes(x = Index, y = Estimated_Variance)) +
    geom_point(color = "blue", alpha = 0.6) +
    #geom_smooth(se = FALSE, color = "darkblue", linetype = "dashed") +
    theme_minimal() +
    labs(
      title = "What is sigma^2 actually capturing?",
      subtitle = "It spikes at the boundaries and high-curvature points where the kernel is biased, not true noise.",
      y = "Estimated Variance (Penalty Weight)"
    )
  
  # Combine plots using patchwork (if installed) or simply print them
  if(requireNamespace("patchwork", quietly = TRUE)) {
    print(p1 / p2)
  } else {
    print(p1)
    print(p2)
  }
}

# Run the plot
set.seed(42)
plt <- plot_zero_noise_variance()
ggsave(filename = here("outputs/experiments/2026-02-23_Understanding_Zero_Variance.pdf"))