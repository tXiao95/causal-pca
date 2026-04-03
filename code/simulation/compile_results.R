library(data.table)
library(gt)
library(ggplot2)
library(here)

# =========================================================================
# 1. Data Processing Function
# =========================================================================

# Helper to load and format specific file types (-error.csv, -diag.csv, -ers.csv)
prep_sim_data <- function(sim_path, file_pattern) {
  files <- list.files(sim_path, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  dt <- rbindlist(lapply(files, fread), fill = TRUE)
  
  # Define ALL methods we need for any table
  all_methods <- c(
    "PCA", "pCCA", "MAVE", "EE", "Oracle-MAVE", "Oracle-EE",
    "RA-MAVE", "RA-EE", "DR-MAVE", "DR-EE", 
    "PO-MAVE", "PO-EE", "RP-MAVE", "RP-EE", "Full_X"
  )
  
  # Subset and factorize
  dt_sub <- dt[method %in% all_methods]
  dt_sub[, method := factor(method, levels = all_methods)]
  
  # Dynamically extract sample sizes and create ordered factors
  available_n <- sort(unique(dt_sub$n))
  factor_levels_n <- paste0("n = ", available_n)
  dt_sub[, n_str := factor(paste0("n = ", n), levels = factor_levels_n)]
  
  return(dt_sub)
}

# Define the core methods for the main tables globally
main_methods <- c("PCA", "pCCA", "MAVE", "Oracle-MAVE", 
                  "RA-MAVE", "DR-MAVE", "PO-MAVE", "RP-MAVE", "Full_X")

# =========================================================================
# 2. Table Generators (Original Tables)
# =========================================================================

create_table_frob_d0 <- function(dt_sub) {
  if(is.null(dt_sub)) return(NULL)
  dt_main <- dt_sub[method %in% main_methods & method != "Full_X"]
  dt_main[, method := droplevels(method)]
  
  dt_error <- dt_main[error_type == "d0", 
                      .(mean_err = mean(frob_norm, na.rm = TRUE),
                        sd_err   = sd(frob_norm, na.rm = TRUE)), 
                      by = .(method, n_str)]
  
  dt_error[, val_formatted := sprintf("%.3f (%.3f)", mean_err, sd_err)]
  tab_wide <- dcast(dt_error, method ~ n_str, value.var = "val_formatted", drop = FALSE)
  
  tab_wide |>
    gt() |>
    tab_header(
      title = "Frobenius Norm Error (Evaluated at d0)",
      subtitle = "Mean (SD) between estimated and true subspace"
    ) |>
    cols_label(method = "Method") |>
    sub_missing(missing_text = "-") |>
    opt_align_table_header(align = "left")
}

create_table_frob_dhat <- function(dt_sub) {
  if(is.null(dt_sub)) return(NULL)
  dt_main <- dt_sub[method %in% main_methods & method != "Full_X"]
  dt_main[, method := droplevels(method)]
  
  dt_err_dhat <- dt_main[error_type == "dhat", 
                         .(mean_err = mean(frob_norm, na.rm = TRUE),
                           sd_err   = sd(frob_norm, na.rm = TRUE)), 
                         by = .(method, n_str)]
  
  dt_err_dhat[, val_formatted := sprintf("%.3f (%.3f)", mean_err, sd_err)]
  dt_err_dhat[is.na(mean_err), val_formatted := NA_character_]
  
  tab_wide <- dcast(dt_err_dhat, method ~ n_str, value.var = "val_formatted", drop = FALSE)
  
  tab_wide |>
    gt() |>
    tab_header(
      title = "Frobenius Norm Error (Evaluated at dhat)",
      subtitle = "Mean (SD) between estimated and true subspace using the estimated dimension"
    ) |>
    cols_label(method = "Method") |>
    sub_missing(missing_text = "-") |>
    opt_align_table_header(align = "left")
}

create_table_dhat <- function(dt_sub) {
  if(is.null(dt_sub)) return(NULL)
  dt_main <- dt_sub[method %in% main_methods]
  dt_main[, method := droplevels(method)]
  
  dim_methods <- c("MAVE", "Oracle-MAVE", "RA-MAVE", "DR-MAVE", "PO-MAVE", "RP-MAVE")
  dt_dim <- dt_main[error_type == "d0" & method %in% dim_methods]
  dt_dim[, method := droplevels(method)]
  
  correct_method_order <- levels(dt_dim$method)
  
  dt_dim_counts <- dt_dim[, .(count = .N), by = .(n_str, method, dhat)]
  dt_dim_counts[, prop := count / sum(count), by = .(n_str, method)]
  
  all_dhats <- unique(dt_dim_counts[!is.na(dhat)]$dhat)
  req_dhats <- sort(union(c(1, 2), all_dhats))
  
  grid <- CJ(n_str = levels(dt_dim$n_str), method = levels(dt_dim$method), dhat = req_dhats, sorted = FALSE)
  dt_dim_full <- merge(grid, dt_dim_counts[!is.na(dhat)], by = c("n_str", "method", "dhat"), all.x = TRUE)
  
  dt_dim_full[is.na(prop), prop := 0]
  dt_dim_full[, val_formatted := sprintf("%.2f", prop)]
  
  dt_dim_full[, method := factor(method, levels = correct_method_order)]
  dt_dim_full[, n_str := factor(n_str, levels = levels(dt_dim$n_str))]
  
  tab_wide <- dcast(dt_dim_full, n_str + method ~ dhat, value.var = "val_formatted", drop = FALSE)
  setorder(tab_wide, n_str, method)
  
  d_cols <- as.character(req_dhats)
  setnames(tab_wide, d_cols, paste0("d = ", d_cols))
  
  tab_wide |>
    gt(groupname_col = "n_str") |>
    row_group_order(groups = levels(dt_dim$n_str)) |> 
    tab_header(
      title = "Dimension Selection Distribution",
      subtitle = "Proportion of iterations selecting each dimension"
    ) |>
    cols_label(method = "Method") |>
    sub_missing(missing_text = "-") |>
    opt_align_table_header(align = "left")
}

create_table_runtime <- function(dt_sub) {
  if(is.null(dt_sub)) return(NULL)
  dt_main <- dt_sub[method %in% main_methods]
  dt_main[, method := droplevels(method)]
  
  dt_time <- dt_main[error_type == "d0", 
                     .(mean_time = mean(time, na.rm = TRUE) / 60,
                       sd_time   = sd(time, na.rm = TRUE) / 60), 
                     by = .(method, n_str)]
  
  dt_time[, val_formatted := sprintf("%.3f (%.3f)", mean_time, sd_time)]
  tab_wide <- dcast(dt_time, method ~ n_str, value.var = "val_formatted", drop = FALSE)
  
  tab_wide |>
    gt() |>
    tab_header(
      title = "Computation Time (Minutes)",
      subtitle = "Mean (SD) runtime per method across iterations"
    ) |>
    cols_label(method = "Method") |>
    sub_missing(missing_text = "-") |>
    opt_align_table_header(align = "left")
}

create_table_ee_comparison <- function(dt_sub, eval_type = "d0") {
  if(is.null(dt_sub)) return(NULL)
  comp_methods <- c("Oracle-MAVE", "Oracle-EE", "RA-MAVE", "RA-EE", "DR-MAVE", "DR-EE", 
                    "PO-MAVE", "PO-EE", "RP-MAVE", "RP-EE")
  
  dt_comp <- dt_sub[method %in% comp_methods & error_type == eval_type]
  if(nrow(dt_comp) == 0) return(NULL)
  
  dt_comp[, c("Estimator", "Algorithm") := tstrsplit(as.character(method), "-")]
  dt_comp[, Estimator := factor(Estimator, levels = c("Oracle", "RA", "DR", "PO", "RP"))]
  dt_comp[, Algorithm := factor(Algorithm, levels = c("MAVE", "EE"))]
  
  dt_error <- dt_comp[, .(mean_err = mean(frob_norm, na.rm = TRUE),
                          sd_err   = sd(frob_norm, na.rm = TRUE)), 
                      by = .(Estimator, Algorithm, n_str)]
  
  dt_error[, val_formatted := sprintf("%.3f (%.3f)", mean_err, sd_err)]
  dt_error[is.na(mean_err), val_formatted := NA_character_]
  
  tab_wide <- dcast(dt_error, Estimator ~ n_str + Algorithm, value.var = "val_formatted", sep = "_")
  
  gt_obj <- tab_wide |> gt()
  n_levels <- levels(dt_error$n_str)
  for (n_val in n_levels) {
    col_m <- paste0(n_val, "_MAVE")
    col_e <- paste0(n_val, "_EE")
    
    if (col_m %in% names(tab_wide) && col_e %in% names(tab_wide)) {
      gt_obj <- gt_obj |>
        tab_spanner(label = n_val, columns = c(col_m, col_e)) |>
        cols_label(.list = setNames(list("MAVE", "EE"), c(col_m, col_e)))
    }
  }
  
  gt_obj |>
    tab_header(
      title = paste0("MAVE vs. EE Comparison (Evaluated at ", eval_type, ")"),
      subtitle = "Frobenius Norm Error: Mean (SD)"
    ) |>
    cols_label(Estimator = "Method") |>
    sub_missing(missing_text = "-") |>
    opt_align_table_header(align = "left")
}

# =========================================================================
# 3. Table Generators (NEW: Diagonals and ERS Evaluation)
# =========================================================================

create_table_diag <- function(dt_diag, true_diag_vals = NULL) {
  if(is.null(dt_diag)) return(NULL)
  dt_main <- dt_diag[method %in% main_methods & method != "Full_X"]
  dt_main[, method := droplevels(method)]
  
  dt_agg <- dt_main[, .(mean_val = mean(P_diag_val, na.rm = TRUE),
                        sd_val   = sd(P_diag_val, na.rm = TRUE)), 
                    by = .(method, n_str, dimension)]
  
  dt_agg[, val_formatted := sprintf("%.3f (%.3f)", mean_val, sd_val)]
  
  tab_wide <- dcast(dt_agg, n_str + method ~ dimension, value.var = "val_formatted", drop = FALSE)
  setnames(tab_wide, as.character(unique(dt_agg$dimension)), paste0("Dim_", unique(dt_agg$dimension)))
  
  if (!is.null(true_diag_vals)) {
    n_levels <- levels(dt_main$n_str)
    true_rows <- lapply(n_levels, function(n_val) {
      row <- c(list(n_str = n_val, method = "True (Analytical)"), as.list(sprintf("%.3f (0.000)", true_diag_vals)))
      names(row) <- colnames(tab_wide)
      as.data.table(row)
    })
    tab_wide <- rbindlist(c(list(tab_wide), true_rows), use.names = TRUE)
    tab_wide[, method := factor(method, levels = c("True (Analytical)", levels(dt_main$method)))]
    setorder(tab_wide, n_str, method)
  }
  
  tab_wide |>
    gt(groupname_col = "n_str") |>
    tab_header(
      title = "Diagonal Elements of Projection Matrix P(\U03B2)",
      subtitle = "Mean (SD) across simulations"
    ) |>
    cols_label(method = "Method") |>
    sub_missing(missing_text = "-") |>
    opt_align_table_header(align = "left")
}

create_table_ers_pointwise <- function(dt_ers, metric = c("bias", "rmse", "coverage")) {
  if(is.null(dt_ers)) return(NULL)
  metric <- match.arg(metric)
  
  dt_main <- dt_ers[method %in% main_methods]
  dt_main[, method := droplevels(method)]
  
  # Calculate metrics across simulations FOR EACH specific evaluation point
  dt_pointwise <- dt_main[, .(
    mean_bias = mean(est - mu_true, na.rm = TRUE),       # Actual Bias (Not Absolute)
    sd_bias   = sd(est - mu_true, na.rm = TRUE),         # Empirical Standard Error
    rmse      = sqrt(mean((est - mu_true)^2, na.rm = TRUE)), # Exact RMSE at this point
    coverage  = mean(mu_true >= ci_lower & mu_true <= ci_upper, na.rm = TRUE) # Coverage Proportion
  ), by = .(eval_id, method, n_str)]
  
  # Format based on the selected metric
  if (metric == "bias") {
    dt_pointwise[, val_formatted := sprintf("%.3f (%.3f)", mean_bias, sd_bias)]
    title_text <- "Point-Wise ERS Estimation: Actual Bias"
    subtitle_text <- "Mean Error (SD of Error) across simulations at each coordinate"
  } else if (metric == "rmse") {
    dt_pointwise[, val_formatted := sprintf("%.3f", rmse)]
    title_text <- "Point-Wise ERS Estimation: RMSE"
    subtitle_text <- "Root Mean Squared Error across simulations at each coordinate"
  } else if (metric == "coverage") {
    dt_pointwise[, val_formatted := sprintf("%.3f", coverage)]
    title_text <- "Point-Wise ERS Estimation: 95% CI Coverage"
    subtitle_text <- "Proportion of simulations where CI covered the true mean at each coordinate"
  }
  
  dt_pointwise[grepl("NA", val_formatted), val_formatted := NA_character_]
  
  # Pivot Wide: Rows = Evaluation Point & Method, Columns = Sample Size
  tab_wide <- dcast(dt_pointwise, eval_id + method ~ n_str, value.var = "val_formatted", drop = FALSE)
  setorder(tab_wide, eval_id, method)
  
  tab_wide |>
    gt(groupname_col = "eval_id") |>
    tab_header(title = title_text, subtitle = subtitle_text) |>
    cols_label(method = "Method") |>
    sub_missing(missing_text = "-") |>
    opt_align_table_header(align = "left")
}

plot_ers_pointwise <- function(dt_ers, dt_true_z, metric = c("bias", "rmse", "coverage")) {
  if(is.null(dt_ers)) return(NULL)
  metric <- match.arg(metric)
  
  dt_main <- dt_ers[method %in% main_methods]
  dt_main[, method := droplevels(method)]
  
  # Step 1: Calculate point-wise metrics
  dt_pointwise <- dt_main[, .(
    bias     = mean(est - mu_true, na.rm = TRUE),
    rmse     = sqrt(mean((est - mu_true)^2, na.rm = TRUE)),
    coverage = mean(mu_true >= ci_lower & mu_true <= ci_upper, na.rm = TRUE)
  ), by = .(eval_id, method, n_str)]
  
  # Step 2: Merge the true Z coordinates based on the eval_id (1 to 100)
  dt_pointwise <- merge(dt_pointwise, dt_true_z, by = "eval_id", all.x = TRUE)
  
  # Save to disk
  out_dir <- here("results/jasa-initial-submission", NAME, EXPERIMENT)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  fwrite(dt_pointwise, here(out_dir, "ers_table.csv"))
  
  # Step 3: Build the ggplot
  p <- ggplot(dt_pointwise, aes(x = Z1_true, y = Z2_true)) +
    facet_grid(n_str ~ method) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.key.width = unit(2, "cm")
    ) +
    labs(x = expression(True~Z[1]), y = expression(True~Z[2]))
  
  # Map aesthetics based on the metric
  if (metric == "bias") {
    max_val <- max(abs(dt_pointwise$bias), na.rm = TRUE)
    p <- p + 
      geom_point(aes(color = bias), size = 2, alpha = 0.8) +
      scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, 
                            limits = c(-max_val, max_val), name = "Bias") +
      ggtitle("Point-Wise Bias across the True 2D Causal Subspace")
    
  } else if (metric == "rmse") {
    p <- p + 
      geom_point(aes(color = rmse), size = 2, alpha = 0.8) +
      scale_color_viridis_c(option = "magma", direction = -1, name = "RMSE") +
      ggtitle("Point-Wise RMSE across the True 2D Causal Subspace")
    
  } else if (metric == "coverage") {
    p <- p + 
      geom_point(aes(color = coverage), size = 2, alpha = 0.8) +
      scale_color_gradientn(
        colours = c("firebrick", "gray90", "dodgerblue"),
        values = c(0, 0.95, 1), # Anchors the gray90 exactly at 0.95
        limits = c(0, 1),
        breaks = c(0, 0.50, 0.80, 0.95),
        labels = c("0.0", "0.5", "0.8", "0.95"),
        name = "Coverage"
      ) +
      ggtitle("Point-Wise 95% CI Coverage across the True 2D Causal Subspace")
  }
  
  return(p)
}

create_table_ers_integrated <- function(dt_ers, metric = c("bias", "rmse", "coverage")) {
  if(is.null(dt_ers)) return(NULL)
  metric <- match.arg(metric)
  
  dt_main <- dt_ers[method %in% main_methods]
  dt_main[, method := droplevels(method)]
  
  # Step 1: Calculate the integrated error FOR EACH SIMULATION RUN across the 100 points
  dt_run <- dt_main[, .(
    run_bias = mean(abs(est - mu_true), na.rm = TRUE),               # Mean Absolute Integrated Bias
    run_rmse = sqrt(mean((est - mu_true)^2, na.rm = TRUE)),          # Root Integrated Squared Error (RISE)
    run_cov  = mean(mu_true >= ci_lower & mu_true <= ci_upper, na.rm = TRUE) # Curve Coverage Proportion
  ), by = .(ID, method, n_str)]
  
  # Step 2: Calculate Mean (SD) of these curve-level metrics across the simulations
  dt_agg <- dt_run[, .(
    mean_bias = mean(run_bias, na.rm = TRUE), sd_bias = sd(run_bias, na.rm = TRUE),
    mean_rmse = mean(run_rmse, na.rm = TRUE), sd_rmse = sd(run_rmse, na.rm = TRUE),
    mean_cov  = mean(run_cov, na.rm = TRUE),  sd_cov  = sd(run_cov, na.rm = TRUE)
  ), by = .(method, n_str)]
  
  # Format based on the selected metric
  if (metric == "bias") {
    dt_agg[, val_formatted := sprintf("%.3f (%.3f)", mean_bias, sd_bias)]
    title_text <- "Integrated ERS Estimation: Mean Absolute Bias"
    subtitle_text <- "Mean (SD) of the integrated absolute error across simulation runs"
  } else if (metric == "rmse") {
    dt_agg[, val_formatted := sprintf("%.3f (%.3f)", mean_rmse, sd_rmse)]
    title_text <- "Integrated ERS Estimation: Root Integrated Squared Error (RISE)"
    subtitle_text <- "Mean (SD) of the integrated RMSE across simulation runs"
  } else if (metric == "coverage") {
    dt_agg[, val_formatted := sprintf("%.3f (%.3f)", mean_cov, sd_cov)]
    title_text <- "Integrated ERS Estimation: Curve Coverage"
    subtitle_text <- "Mean (SD) proportion of the 100-point curve covered by the 95% CI per run"
  }
  
  dt_agg[grepl("NA", val_formatted), val_formatted := NA_character_]
  
  # Pivot Wide: Rows = Method, Columns = Sample Size
  tab_wide <- dcast(dt_agg, method ~ n_str, value.var = "val_formatted", drop = FALSE)
  
  tab_wide |>
    gt() |>
    tab_header(title = title_text, subtitle = subtitle_text) |>
    cols_label(method = "Method") |>
    sub_missing(missing_text = "-") |>
    opt_align_table_header(align = "left")
}

main <- function(){
  # Source the simulation DGP so we can perfectly reconstruct the true Z
  source(here("R/simulate_data.R"))
  
  sim_path <- here("outputs", "simulation", "jasa-initial-submission", NAME, EXPERIMENT)
  
  
  # Load the 3 distinct data streams
  dt_err  <- prep_sim_data(sim_path, "-error\\.csv$")
  dt_diag <- prep_sim_data(sim_path, "-diag\\.csv$")
  dt_ers  <- prep_sim_data(sim_path, "-ers\\.csv$")
  
  # -----------------------------------------------------------------------
  # 1. Reconstruct the True Z Evaluation Grid mathematically (seed 99999....should modify this later)
  # -----------------------------------------------------------------------
  message("Reconstructing true Z grid for EXPERIMENT: ", EXPERIMENT)
  set.seed(99999) # The exact hardcoded seed used in fit_simulation.R
  
  is_weak   <- (EXPERIMENT == "weak_dim")
  int_coef  <- if(EXPERIMENT == "additive") 0.0 else 5.0
  is_sparse <- if(NAME == "main_paper_final_results_nnet") TRUE else FALSE
  
  # Calling the DGP perfectly reconstructs the grid
  sim_grid <- simulate_causal_sdr_simple(n = 100, p = 10, q = 5, noise_sd = 0.5, 
                                         rho_X = 0.8, 
                                         interaction_coef = int_coef, 
                                         weak_dim_signal = is_weak, 
                                         sparse = is_sparse)
  
  dt_true_z <- data.table(
    eval_id = 1:100, 
    Z1_true = sim_grid$Z[, 1], 
    Z2_true = sim_grid$Z[, 2]
  )
  
  # -----------------------------------------------------------------------
  # 2. Generate Original Structural & Runtime Tables
  # -----------------------------------------------------------------------
  gt_table_frob_d0   <- create_table_frob_d0(dt_err)
  gt_table_frob_dhat <- create_table_frob_dhat(dt_err)
  gt_table_dhat      <- create_table_dhat(dt_err)
  gt_table_runtime   <- create_table_runtime(dt_err)
  gt_comp_d0         <- create_table_ee_comparison(dt_err, eval_type = "d0")
  gt_comp_dhat       <- create_table_ee_comparison(dt_err, eval_type = "dhat")
  
  # -----------------------------------------------------------------------
  # 3. Generate Diagonal Tables
  # -----------------------------------------------------------------------
  true_diag_vector <- NULL 
  gt_table_diag <- create_table_diag(dt_diag, true_diag_vals = true_diag_vector)
  
  # -----------------------------------------------------------------------
  # 4. Generate ERS Plots and Integrated Tables
  # -----------------------------------------------------------------------
  # Point-wise Plots (passing the dt_true_z mapping!)
  plot_ers_pw_bias <- plot_ers_pointwise(dt_ers, dt_true_z, metric = "bias")
  plot_ers_pw_rmse <- plot_ers_pointwise(dt_ers, dt_true_z, metric = "rmse")
  plot_ers_pw_cov  <- plot_ers_pointwise(dt_ers, dt_true_z, metric = "coverage")
  
  # Integrated Tables
  gt_table_ers_int_bias <- create_table_ers_integrated(dt_ers, metric = "bias")
  gt_table_ers_int_rmse <- create_table_ers_integrated(dt_ers, metric = "rmse")
  gt_table_ers_int_cov  <- create_table_ers_integrated(dt_ers, metric = "coverage")
  
  # -----------------------------------------------------------------------
  # Save Output
  # -----------------------------------------------------------------------
  out_dir <- here("results/jasa-initial-submission", NAME, EXPERIMENT)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Save Plots to PDF (width and height adjusted for wide facet grids)
  if(!is.null(plot_ers_pw_bias)) ggsave(file.path(out_dir, "plot_pw_bias.pdf"), plot_ers_pw_bias, width = 16, height = 10)
  if(!is.null(plot_ers_pw_rmse)) ggsave(file.path(out_dir, "plot_pw_rmse.pdf"), plot_ers_pw_rmse, width = 16, height = 10)
  if(!is.null(plot_ers_pw_cov))  ggsave(file.path(out_dir, "plot_pw_coverage.pdf"), plot_ers_pw_cov, width = 16, height = 10)
  
  table_list <- list(
    frob_d0      = gt_table_frob_d0,
    frob_dhat    = gt_table_frob_dhat,
    dhat         = gt_table_dhat,
    runtime      = gt_table_runtime,
    comp_d0      = gt_comp_d0,
    comp_dhat    = gt_comp_dhat,
    diag         = gt_table_diag,
    ers_int_bias = gt_table_ers_int_bias,
    ers_int_rmse = gt_table_ers_int_rmse,
    ers_int_cov  = gt_table_ers_int_cov
  )
  
  saveRDS(table_list, file.path(out_dir, "simulation_gt_tables.rds"))
  message("All tables and plots successfully compiled and saved to: ", out_dir)
}

# Main code ---------------------------------------------------------------
if(interactive()){
  #NAME       <- "main_paper_final_results_nonsparse_nnet_128x64x32"
  NAME       <- "main_paper_final_results_nnet"
  EXPERIMENT <- "interaction"
} else{
  NAME       <- as.character( commandArgs(trailingOnly = TRUE)[1] )
  EXPERIMENT <- as.character( commandArgs(trailingOnly = TRUE)[2] )
}

main()