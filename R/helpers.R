# ---- Helpers for analysis ----
# reshaping dataset from long to wide format, by adding attributes to each variable group
make_wide <- function(df, id_col, time_col, base_cols, outcome, ..., static = FALSE, waves = NULL) {
  
  require(rlang)
  require(dplyr)
  require(tidyr)

  outcome_name <- as_name(ensym(outcome))
  
  t_conf <- enquos(...)
  
  if (!is.null(waves)) {
    df <- df |> dplyr::filter({{time_col}} %in% waves)
  }

  t_max <- max(df |> dplyr::pull({{time_col}}))
  t_min <- min(df |> dplyr::pull({{time_col}}))

  df_out <- df |>
    dplyr::select({{id_col}}, {{time_col}}, {{base_cols}}, !!!t_conf, {{outcome}}) |>
    tidyr::pivot_wider(
      id_cols = c({{id_col}}, {{base_cols}}),
      names_from = {{time_col}},
      values_from = -c({{id_col}}, {{time_col}}, {{base_cols}})) |>
    dplyr::select(
#      {{id_col}},
      {{base_cols}},
      ends_with("0"), # wave 2
      ends_with("1"), # wave 3
      ends_with("2"), # wave 4
      ends_with("3"), # wave 5
      ends_with("4"), # wave 6
      ends_with("5"), # wave 7
      ends_with("6"), # wave 8
      ends_with("7"), # wave 9
      ends_with("8"), # wave 10
      ends_with("9"), # wave 11
      ends_with("10"), # wave 12
      ends_with("11"), # wave 13
      ends_with("12"), # wave 14
      ends_with("13") # wave 15
    )

  if (static) {
      intermediate_cols <- paste0(outcome_name, "_", (t_min + 1):(t_max - 1))
      df_out <- df_out |> dplyr::select(-any_of(intermediate_cols))
    }

  attr(df_out, "n_tvars") <- length(t_conf)
  attr(df_out, "n_base") <- length(df[1,] |> dplyr::select({{base_cols}}))
  attr(df_out, "baseline_vars") <- names(dplyr::select(df_out, {{base_cols}}))
  attr(df_out, "time_lagged")   <- names(dplyr::select(df_out, contains("lagged")))
  attr(df_out, "outcome_var")   <- paste0(outcome_name, "_", t_max)
  attr(df_out, "outcome_baseline") <- paste0(outcome_name, "_base")
  attr(df_out, "time_points") <- length(unique(df |> dplyr::pull({{time_col}})))
  attr(df_out, "intermediate_vars") <- names(dplyr::select(df_out, !contains("lagged"), -c({{base_cols}})))

  
  return(df_out)
}

# tag the exposure column(s) so downstream code (predictor matrix) knows which nodes are treatment.
set_exposure <- function(df, exposure) {
  require(rlang)
  exposure <- as_name(ensym(exposure))
  exposure_vars <- grep(paste0("^", exposure, "_\\d+$"), names(df), value = TRUE)
  if (length(exposure_vars) == 0L) stop("no exposure columns matched '", exposure, "_<wave>'")
  attr(df, "exposure_vars") <- exposure_vars
  return(df)
}

# function to replicate DAG arrows (counterfactual imputation for gFormulaMI)
make_counterfactual_matrix <- function(return_vals) {
  n_tvars <- attr(return_vals, "n_tvars")
  n_base <- attr(return_vals, "n_base")
  baseline_vars <- attr(return_vals, "baseline_vars")
  time_lagged <- attr(return_vals, "time_lagged")
  outcome_var <- attr(return_vals, "outcome_var")
  outcome_baseline <- attr(return_vals, "outcome_baseline")
  time_points <- attr(return_vals, "time_points")
  intermediate_vars <- attr(return_vals, "intermediate_vars")
  exposure <- attr(return_vals, "exposure_vars")

  n_vars <- length(return_vals) + 1
  
  p_mat <- matrix(0, ncol = n_vars, nrow = n_vars)
  
  nodes <- c(colnames(return_vals), "regime")
  dimnames(p_mat) <- list(nodes, nodes)

  # time-invariant confounders are lower-triangular among themselves
  # allow gFormulaMI to draw these from their observed values instead of the intercept
  base_ordered <- intersect(nodes, baseline_vars)
  for (i in seq_along(base_ordered)[-1]) {
    p_mat[base_ordered[i], base_ordered[seq_len(i - 1L)]] <- 1
  }

  # baseline outcome is predicted by all time-invariant confounders
  p_mat[outcome_baseline, setdiff(c(baseline_vars, outcome_baseline), outcome_baseline)] <- 1

  # exposure is predicted by all time-invariant confounders and baseline outcome
  # if it's time-varying, only the most recent exposure is predicted by baseline confounders and baseline outcome
  if (time_points <= 1){
    p_mat[exposure, c(baseline_vars, outcome_baseline)] <- 1
  } else {
    p_mat[exposure, c(baseline_vars)] <- 1
  }
  
  # time-lagged confounders are predicted by baseline outcome and all time-invariant confounders
  p_mat[time_lagged, c(baseline_vars)] <- 1
  
  # time-lagged confounders predict exposure (single-time) or lagged_t predicts exposure_t (time-varying)
  if (time_points <= 1) {
    p_mat[exposure, time_lagged] <- 1
  } else {
    for (t in seq_len(time_points - 1)) {
      exp_t <- exposure[endsWith(exposure, paste0("_", t))]
      lag_t <- time_lagged[endsWith(time_lagged, paste0("_lagged_", t))]
      p_mat[exp_t, lag_t] <- 1
    }
  }

  # time-varying covariates that are neither the exposure, the outcome, nor a lag. 
  covar_vars <- setdiff(intermediate_vars, c(exposure, outcome_var, outcome_baseline))

  if (time_points <= 1 && length(covar_vars) > 0) {
    # gFormulaMI runs mice with maxit = 1 and the default (column-order) visit
    # sequence, so a predictor is only useful if it precedes its target.
    for (v in covar_vars) {
      earlier <- nodes[seq_len(match(v, nodes) - 1L)]
      p_mat[v, intersect(c(baseline_vars, outcome_baseline, time_lagged), earlier)] <- 1
    }
    # they confound the exposure–outcome relationship
    p_mat[exposure, covar_vars]    <- 1
    p_mat[outcome_var, covar_vars] <- 1
  }

  # outcome is predicted by all time-invariant confounders, all time-lagged confounders, baseline outcome,
  # and exposure (single-time) or exposure_t (time-varying)
  p_mat[outcome_var, c(baseline_vars, time_lagged, outcome_baseline)] <- 1

  if (time_points <= 1) {
    p_mat[outcome_var, exposure] <- 1
  } else {
    for (t in seq_len(time_points - 1)) {
      exp_t <- exposure[endsWith(exposure, paste0("_", t))]
      p_mat[outcome_var, exp_t] <- 1
    }
  }

  # time-varying confounders are predicted by their immediate prior period
  # same for outcomes, predictors, lagged confounders
  for (t in seq_len(time_points - 1)) {
    prev_c <- intermediate_vars[endsWith(intermediate_vars, paste0("_", t))]
    curr_c <- intermediate_vars[endsWith(intermediate_vars, paste0("_", t + 1))]
    p_mat[curr_c, prev_c] <- 1
    
    prev_o <- paste0(outcome_var, "_", t)
    curr_o <- paste0(outcome_var, "_", t + 1)
    p_mat[curr_o, prev_o] <- 1

    prev_lag <- time_lagged[endsWith(time_lagged, paste0("_lagged_", t))]
    curr_lag <- time_lagged[endsWith(time_lagged, paste0("_lagged_", t + 1))]
    p_mat[curr_lag, prev_lag] <- 1
  }
  # regime node is predicted by all other variables (i.e., all confounders and the outcome)
  p_mat["regime", 1:(n_vars - 1)] <- 1

  return(p_mat)
}

# Make factor levels safe to paste into a model formula.
#
# mice's `formulas` mode (R/run_mice.R) builds design-matrix column names from
# factor level labels and re-parses them as formula terms, so a level containing
# a space, hyphen or "+" ("Non-white", "No benefits", "2+") either fails to parse
# or yields a mangled term.
#
# Rewrites a variable's levels ONLY when some level contains a character outside
# [A-Za-z0-9._].
#

sanitize_factor_levels <- function(df, skip = character()) {
  hostile <- "[^A-Za-z0-9._]"

  is_dt <- data.table::is.data.table(df)
  if (is_dt) df <- data.table::copy(df)   # never modify the caller's table by reference

  fct_cols <- setdiff(names(df)[purrr::map_lgl(df, is.factor)], skip)

  for (v in fct_cols) {
    lv <- levels(df[[v]])
    if (!any(grepl(hostile, lv))) next
    new_labels <- make.names(lv, unique = TRUE)
    new_v <- factor(new_labels[as.integer(df[[v]])],
                     levels  = sort(new_labels),
                     ordered = is.ordered(df[[v]]))
    if (is_dt) data.table::set(df, j = v, value = new_v) else df[[v]] <- new_v
  }

  df
}

# Imputation methods for the counterfactual (synthetic) block drawn by
# gFormulaMI::gFormulaImpute().
make_counterfactual_method <- function(return_vals) {

  probe <- rbind(as.data.frame(return_vals), NA)

  method <- mice::make.method(probe,
                              defaultMethod = c("pmm", "logreg", "polyreg", "polr"))

  # treatment columns are assigned from the regime rather than imputed, and are
  # complete in the observed block, so gFormulaMI expects them blank
  exposure <- attr(return_vals, "exposure_vars")
  if (is.null(exposure)) {
    stop("make_counterfactual_method: 'exposure_vars' attribute missing -- ",
         "pass wide_data through set_exposure() first.")
  }
  method[exposure] <- ""

  return(method)
}
