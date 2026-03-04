library(data.table)
library(gt)
library(ggplot2)
library(here)

experiments <- c("baseline",
                 "additive-confounding",
                 "high-SNR",
                 "weak-dim-signal")

# =========================================================================
# 1. Data Processing Function
# =========================================================================

prep_sim_data <- function(sim_path) {
  # Read all files
  dt <- lapply(list.files(sim_path, full.names = TRUE), fread) |>
    rbindlist(fill = TRUE)
  
  # Define ALL methods we need for any table
  all_methods <- c(
    "PCA", "pCCA", "MAVE", "EE", "Oracle-MAVE", "Oracle-EE",
    "RA-MAVE", "RA-EE", "DR-MAVE", "DR-EE", 
    "PO-MAVE", "PO-EE", "RP-MAVE", "RP-EE"
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

# Define the core methods for the main tables globally so we don't repeat it
main_methods <- c("PCA", "pCCA", "MAVE", "Oracle-MAVE", 
                  "RA-MAVE", "DR-MAVE", "PO-MAVE", "RP-MAVE")

# =========================================================================
# 2. Table Generators (Main Tables)
# =========================================================================

create_table_frob_d0 <- function(dt_sub) {
  dt_main <- dt_sub[method %in% main_methods]
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
  dt_main <- dt_sub[method %in% main_methods]
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
  dt_main <- dt_sub[method %in% main_methods]
  dt_main[, method := droplevels(method)]
  
  dim_methods <- c("MAVE", "Oracle-MAVE", "RA-MAVE", "DR-MAVE", "PO-MAVE", "RP-MAVE")
  dt_dim <- dt_main[error_type == "d0" & method %in% dim_methods]
  dt_dim[, method := droplevels(method)]
  
  # Save the correct factor levels before they get converted to characters
  correct_method_order <- levels(dt_dim$method)
  
  dt_dim_counts <- dt_dim[, .(count = .N), by = .(n_str, method, dhat)]
  dt_dim_counts[, prop := count / sum(count), by = .(n_str, method)]
  
  all_dhats <- unique(dt_dim_counts[!is.na(dhat)]$dhat)
  req_dhats <- sort(union(c(1, 2), all_dhats))
  
  grid <- CJ(n_str = levels(dt_dim$n_str),
             method = levels(dt_dim$method),
             dhat = req_dhats, 
             sorted = FALSE)
  
  dt_dim_full <- merge(grid, dt_dim_counts[!is.na(dhat)], by = c("n_str", "method", "dhat"), all.x = TRUE)
  
  dt_dim_full[is.na(prop), prop := 0]
  dt_dim_full[, val_formatted := sprintf("%.2f", prop)]
  
  # CRITICAL FIX: Re-apply the strict factor levels before casting and sorting
  dt_dim_full[, method := factor(method, levels = correct_method_order)]
  dt_dim_full[, n_str := factor(n_str, levels = levels(dt_dim$n_str))]
  
  tab_wide <- dcast(dt_dim_full, n_str + method ~ dhat, value.var = "val_formatted", drop = FALSE)
  
  # Sort by sample size, then by the preserved method order
  setorder(tab_wide, n_str, method)
  
  d_cols <- as.character(req_dhats)
  new_d_cols <- paste0("d = ", d_cols)
  setnames(tab_wide, d_cols, new_d_cols)
  
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

# =========================================================================
# 3. Table Generator (MAVE vs EE Comparison)
# =========================================================================

create_table_ee_comparison <- function(dt_sub, eval_type = "d0") {
  
  # Target only the causal estimators
  comp_methods <- c("Oracle-MAVE", "Oracle-EE", "RA-MAVE", "RA-EE", "DR-MAVE", "DR-EE", 
                    "PO-MAVE", "PO-EE", "RP-MAVE", "RP-EE")
  
  dt_comp <- dt_sub[method %in% comp_methods & error_type == eval_type]
  
  # Split 'method' string into the Causal Estimator and the SDR Algorithm
  dt_comp[, c("Estimator", "Algorithm") := tstrsplit(as.character(method), "-")]
  dt_comp[, Estimator := factor(Estimator, levels = c("Oracle", "RA", "DR", "PO", "RP"))]
  dt_comp[, Algorithm := factor(Algorithm, levels = c("MAVE", "EE"))]
  
  dt_error <- dt_comp[, .(mean_err = mean(frob_norm, na.rm = TRUE),
                          sd_err   = sd(frob_norm, na.rm = TRUE)), 
                      by = .(Estimator, Algorithm, n_str)]
  
  dt_error[, val_formatted := sprintf("%.3f (%.3f)", mean_err, sd_err)]
  dt_error[is.na(mean_err), val_formatted := NA_character_]
  
  # Pivot to wide format: Rows = Estimator, Cols = n_str_Algorithm
  tab_wide <- dcast(dt_error, Estimator ~ n_str + Algorithm, value.var = "val_formatted", sep = "_")
  
  # Build GT table
  gt_obj <- tab_wide |> gt()
  
  # Dynamically add spanners for each sample size and rename the underlying columns to just "MAVE" / "EE"
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
# 4. Main Execution Script
# =========================================================================

# Define directory and load data
dir      <- paste0("jasa-initial-submission")
sim_path <- here("outputs", "simulation", dir)
dt_sub   <- prep_sim_data(sim_path)

# Generate and print main tables
gt_table_frob_d0 <- create_table_frob_d0(dt_sub)
print(gt_table_frob_d0)

gt_table_frob_dhat <- create_table_frob_dhat(dt_sub)
print(gt_table_frob_dhat)

gt_table_dhat <- create_table_dhat(dt_sub)
print(gt_table_dhat)

gt_table_runtime <- create_table_runtime(dt_sub)
print(gt_table_runtime)

# Generate and print EE Comparison tables
gt_comp_d0 <- create_table_ee_comparison(dt_sub, eval_type = "d0")
print(gt_comp_d0)

gt_comp_dhat <- create_table_ee_comparison(dt_sub, eval_type = "dhat")
print(gt_comp_dhat)

# Save all objects
table_list <- list(frob_d0   = gt_table_frob_d0,
                   frob_dhat = gt_table_frob_dhat,
                   dhat      = gt_table_dhat,
                   runtime   = gt_table_runtime,
                   comp_d0   = gt_comp_d0,
                   comp_dhat = gt_comp_dhat)

# Ensure results directory exists
out_dir <- here("results/jasa-initial-submission")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saveRDS(table_list, file.path(out_dir, "simulation_gt_tables.rds"))