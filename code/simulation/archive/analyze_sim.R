library(data.table)
library(dplyr)
library(ggplot2)
library(gt)
library(stringr)
library(tidyr)

files <- list.files(c(here("results/sim_results")), 
                      full.names = TRUE, pattern = "*.csv")
                    
# Extract the number after "sim-" and convert to numeric
sim_nums <- str_extract(basename(files), "(?<=sim-)[0-9]+") |> as.numeric()


files1 <- files[sim_nums %in% 1:1119]
files2 <- list.files(c(here("results/sim_results/temp")), 
                      full.names = TRUE, pattern = "*.csv")

files <- c(files1, files2)

dt    <- lapply(files, fread) |> rbindlist()

dt <- dt |>
  mutate(scenario_full = case_when(
    scenario == "1" ~ "1: Linear X, Linear C, Additive Confounding",
    scenario == "2" ~ "2: Quadratic X",
    scenario == "3" ~ "3: + sin(x) Confounding",
    scenario == "4" ~ "4: + Laplace Error in X",
    scenario == "5" ~ "5: + Multiplicative Error in X",
    TRUE ~ "Other"
  ))

dt[, method := ifelse(method == "do_scaled", "do.s", method)]

dt$method <- factor(dt$method, c("Truth", "pCCA", "CCA", 
                                 "cs", "do", "do.s",
                                 "PLS", "pPLS", "PCA"))

# Compute biases
# Get the Truth values by ID
truth_vals <- dt[method == "Truth", .SD, .SDcols = paste0("X", 1:6), by = ID]

# Merge the truth values onto the original table by ID
dt_bias <- merge(dt, truth_vals, by = "ID", suffixes = c("", "_truth"))

# Compute bias columns: difference from Truth
for (j in 1:6) {
  x_col <- paste0("X", j)
  bias_col <- paste0("bias_X", j)
  truth_col <- paste0(x_col, "_truth")
  dt_bias[, (bias_col) := get(x_col) - get(truth_col)]
}

# Summary of bias mean and sd by method, sample size, and scenario
bias_summary <- dt_bias[method != "Truth", .(
  mean_bias_X1 = mean(bias_X1), sd_bias_X1 = sd(bias_X1),
  mean_bias_X2 = mean(bias_X2), sd_bias_X2 = sd(bias_X2),
  mean_bias_X3 = mean(bias_X3), sd_bias_X3 = sd(bias_X3),
  mean_bias_X4 = mean(bias_X4), sd_bias_X4 = sd(bias_X4),
  mean_bias_X5 = mean(bias_X5), sd_bias_X5 = sd(bias_X5),
  mean_bias_X6 = mean(bias_X6), sd_bias_X6 = sd(bias_X6)
), by = .(method, n, scenario)]

# Step 1: Format and prepare bias data
bias_gt_data <- bias_summary %>%
  mutate(across(starts_with("mean_bias_X"), ~sprintf("%.3f", .), .names = "mean_fmt_{.col}"),
         across(starts_with("sd_bias_X"), ~sprintf("%.2f", .), .names = "sd_fmt_{.col}")) %>%
  rowwise() %>%
  mutate(
    X1 = paste0(mean_fmt_mean_bias_X1, " (", sd_fmt_sd_bias_X1, ")"),
    X2 = paste0(mean_fmt_mean_bias_X2, " (", sd_fmt_sd_bias_X2, ")"),
    X3 = paste0(mean_fmt_mean_bias_X3, " (", sd_fmt_sd_bias_X3, ")"),
    X4 = paste0(mean_fmt_mean_bias_X4, " (", sd_fmt_sd_bias_X4, ")"),
    X5 = paste0(mean_fmt_mean_bias_X5, " (", sd_fmt_sd_bias_X5, ")"),
    X6 = paste0(mean_fmt_mean_bias_X6, " (", sd_fmt_sd_bias_X6, ")")
  ) %>%
  ungroup() %>%
  mutate(
    n = factor(n, levels = c(200, 800, 1600)),
    scenario = factor(scenario, levels = sort(unique(scenario))),
    group_label = paste0("Scenario ", scenario, " — n = ", n)
  ) %>%
  arrange(scenario, n, method)

# Keep only relevant columns for gt
gt_data <- bias_gt_data %>%
  select(group_label, method, X1:X6)

# Step 2: Build initial table
bias_gt_table <- gt_data %>%
  gt(rowname_col = "method", groupname_col = "group_label") %>%
  tab_header(title = "Component-wise Bias: Mean (SD)") %>%
  cols_label(
    X1 = "X1", X2 = "X2", X3 = "X3",
    X4 = "X4", X5 = "X5", X6 = "X6"
  )

# Step 3: Bold minimum absolute bias per column within each group
for (col in paste0("X", 1:6)) {
  min_rows <- gt_data %>%
    mutate(abs_bias = as.numeric(sub(" .*", "", !!sym(col)))) %>%
    group_by(group_label) %>%
    filter(abs(abs_bias) == min(abs(abs_bias))) %>%
    ungroup() %>%
    select(group_label, method)
  
  for (i in seq_len(nrow(min_rows))) {
    bias_gt_table <- bias_gt_table %>%
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(
          columns = col,
          rows = method == min_rows$method[i] & group_label == min_rows$group_label[i]
        )
      )
  }
}


# Distribution of Error for each sample size and scenario
for(N in unique(dt$n)){
  plt <- ggplot(dt[method != "PCA" & method != "Truth" & n == N], 
         aes(as.factor(method), error)) + 
    geom_boxplot(size = 0.3, outlier.size = 1) + 
    facet_wrap(~scenario_full, scales = "free") + 
    expand_limits(y = 0) + 
    theme_bw() + 
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    xlab("Method") +
    ylab("Error (Frobenius norm)") + 
    ggtitle(paste0("Sample Size (N=", N, ")"))
  ggsave(paste0("results/figures/sim_new_", N, ".png"), width=11, height=8.5)
}

# Table of runtimes (outcome regression with SL.gam, glm, and glmnet only)
df <- dt[, .(mean = mean(runtime / 60), 
       sd = sd(runtime / 60)), .(method, scenario,n)][!is.na(mean)]
# Step 1: Combine mean and sd into a formatted string
df_fmt <- df %>%
  mutate(label = sprintf("%.2f (%.2f)", mean, sd)) %>%
  select(method, scenario, n, label)

# Step 2: Pivot to wide format
df_wide <- df_fmt %>%
  pivot_wider(
    names_from = n,
    values_from = label,
    names_prefix = "n"
  ) %>%
  mutate(scenario = paste("Scenario", scenario))  # Add "Scenario" prefix

# Step 3: Build gt table
df_final <- df_wide %>%
  arrange(scenario, method) %>%
  gt(rowname_col = "method", groupname_col = "scenario") %>%
  tab_spanner(label = "n = 200", columns = "n200") %>%
  tab_spanner(label = "n = 800", columns = "n800") %>%
  tab_spanner(label = "n = 1600", columns = "n1600") %>%
  cols_label(
    n200  = "",
    n800  = "",
    n1600 = ""
  ) %>%
  tab_header(
    title = "Simulation Runtimes (in minutes)"
  )

df_final

# Save Bias (SD) table
gtsave(bias_gt_table, "results/bias_table.html")
as_latex(bias_gt_table) %>% 
  as.character() %>% 
  writeLines("results/bias_table.tex")

# Save Runtime Table
gtsave(df_final, "results/runtime_table.html")
as_latex(df_final) %>% 
  as.character() %>% 
  writeLines("results/runtime_table.tex")
