# Long to wide format; attributes record each variable group for downstream code.
make_wide <- function(df, id_col, time_col, base_cols, outcome, ..., waves = NULL) {
  require(rlang)
  require(dplyr)
  require(tidyr)

  outcome_name <- as_name(ensym(outcome))
  t_conf <- enquos(...)

  if (!is.null(waves)) {
    df <- df |> dplyr::filter({{time_col}} %in% waves)
  }

  t_max <- max(df |> dplyr::pull({{time_col}}))

  df_out <- df |>
    dplyr::select({{id_col}}, {{time_col}}, {{base_cols}}, !!!t_conf, {{outcome}}) |>
    tidyr::pivot_wider(
      id_cols = c({{id_col}}, {{base_cols}}),
      names_from = {{time_col}},
      values_from = -c({{id_col}}, {{time_col}}, {{base_cols}})) |>
    # drops id_col; groups columns by time suffix. ends_with() matches trailing
    # digits, so _10 sorts with _0.
    dplyr::select({{base_cols}}, ends_with(as.character(0:13)))

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

# Tag the exposure columns (<exposure>_<wave>) so downstream code knows the treatment nodes.
set_exposure <- function(df, exposure) {
  require(rlang)
  exposure <- as_name(ensym(exposure))
  exposure_vars <- grep(paste0("^", exposure, "_\\d+$"), names(df), value = TRUE)
  if (length(exposure_vars) == 0L) stop("no exposure columns matched '", exposure, "_<wave>'")
  attr(df, "exposure_vars") <- exposure_vars
  return(df)
}

# function to replicate DAG arrows (counterfactual imputation for gFormulaMI)
# DAG nodes: TI = baseline_vars, Y0 = outcome_baseline, Cl = time_lagged (E0 =
# econ_dist_bin_lagged among them), Job/Inc = the other time-varying covariates,
# E = exposure_vars, Y1 = outcome_var. The time-varying part repeats at every
# time point t, the integer after a column's last "_".
make_counterfactual_matrix <- function(return_vals) {
  baseline_vars     <- attr(return_vals, "baseline_vars")
  time_lagged       <- attr(return_vals, "time_lagged")
  outcome_var       <- attr(return_vals, "outcome_var")
  outcome_baseline  <- attr(return_vals, "outcome_baseline")
  intermediate_vars <- attr(return_vals, "intermediate_vars")
  exposure          <- attr(return_vals, "exposure_vars")

  nodes <- c(colnames(return_vals), "regime")
  p_mat <- matrix(0, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))

  wave  <- function(x) as.integer(sub(".*_", "", x))
  at    <- function(x, t) x[wave(x) == t]
  waves <- sort(unique(wave(c(time_lagged, intermediate_vars))))

  # outcome scale at earlier time points (sf12mcs_dv_0) is not a confounder: no
  # edges, so it is drawn intercept-only and predicts nothing
  outcome_stem <- sub("_\\d+$", "", outcome_var)
  outcome_prev <- setdiff(grep(paste0("^", outcome_stem, "_\\d+$"), intermediate_vars, value = TRUE),
                          outcome_var)
  # time-varying covariates that are neither the exposure, the outcome, nor a lag.
  covar_vars <- setdiff(intermediate_vars, c(exposure, outcome_var, outcome_prev))

  edge <- function(to, from) data.table::CJ(to = to, from = from)
  # gFormulaMI runs mice with maxit = 1 and the default (column-order) visit
  # sequence, so a predictor is only useful if it precedes its target.
  earlier_only <- function(e) e[match(from, nodes) < match(to, nodes)]

  within_t <- purrr::map(waves, \(t) data.table::rbindlist(list(
    # time-lagged confounders are predicted by baseline outcome and all
    # time-invariant confounders (TI, Y0 -> Cl)
    edge(at(time_lagged, t), baseline_vars),
    # TI, Y0, Cl, E0 -> Job, Inc
    earlier_only(edge(at(covar_vars, t), c(baseline_vars, outcome_baseline, time_lagged))),
    # TI, Y0, Cl, E0, Job, Inc -> E. Documents the DAG only: exposures are set by
    # the regime (method ""), so this row is never used and may point forward.
    edge(at(exposure, t), c(baseline_vars, outcome_baseline, at(time_lagged, t), at(covar_vars, t)))
  )))

  # time-varying confounders and exposure are predicted by their immediate prior
  # period, lagged confounders by the prior lags only (no E -> Cl in the DAG)
  across_t <- purrr::map2(waves[-1], utils::head(waves, -1), \(t, s) data.table::rbindlist(list(
    edge(c(at(covar_vars, t), at(exposure, t)), c(at(covar_vars, s), at(exposure, s))),
    edge(at(time_lagged, t), at(time_lagged, s))
  )))

  e <- data.table::rbindlist(c(
    # time-invariant confounders are lower-triangular among themselves
    # allow gFormulaMI to draw these from their observed values instead of the intercept
    list(earlier_only(edge(baseline_vars, baseline_vars)),
         # baseline outcome is predicted by all time-invariant confounders (TI -> Y0)
         edge(outcome_baseline, setdiff(baseline_vars, outcome_baseline))),
    within_t,
    across_t,
    # outcome is predicted by every confounder and exposure: TI, Y0, Cl, E0, Job, Inc, E -> Y1
    list(edge(outcome_var, c(baseline_vars, outcome_baseline, time_lagged, covar_vars, exposure)),
         # regime node is predicted by all other variables (i.e., all confounders and the outcome)
         edge("regime", setdiff(nodes, "regime")))
  ))
  p_mat[cbind(e$to, e$from)] <- 1

  return(p_mat)
}

# Make factor levels safe to paste into a model formula. mice's `formulas` mode
# (R/run_mice.R) builds design-matrix column names from level labels and
# re-parses them as formula terms, so a level with a space, hyphen or "+"
# ("Non-white", "No benefits", "2+") fails to parse or yields a mangled term.
# Rewrites a variable's levels ONLY when some level has a character outside
# [A-Za-z0-9._].
sanitize_factor_levels <- function(df, skip = character()) {
  hostile <- "[^A-Za-z0-9._]"

  is_dt <- data.table::is.data.table(df)
  if (is_dt) df <- data.table::copy(df)   # never modify the caller's table by reference

  fct_cols <- setdiff(names(df)[purrr::map_lgl(df, is.factor)], skip)

  purrr::walk(fct_cols, \(v) {
    lv <- levels(df[[v]])
    if (!any(grepl(hostile, lv))) return()
    new_labels <- make.names(lv, unique = TRUE)
    new_v <- factor(new_labels[as.integer(df[[v]])],
                    levels  = sort(new_labels),
                    ordered = is.ordered(df[[v]]))
    if (is_dt) data.table::set(df, j = v, value = new_v) else df[[v]] <<- new_v
  })

  df
}

# Imputation methods for the counterfactual (synthetic) block drawn by
# gFormulaMI::gFormulaImpute().
make_counterfactual_method <- function(return_vals) {
  # all-NA probe row: make.method() blanks complete columns, but the synthetic block is all NA
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
