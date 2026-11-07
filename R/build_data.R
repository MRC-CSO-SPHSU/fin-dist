build_data <- function(data, pre_wave = 1, base_wave = 2, start_wave = 3, final_wave = 4,
                       outcome = c("MCS", "PCS")) {
  outcome <- match.arg(outcome)
  # the lags at final_wave are start_wave values only for adjacent waves
  stopifnot(final_wave == start_wave + 1)

  # _base columns: covariates at base_wave, SF-12 at pre_wave; both move with the window
  base_vars <- c("sex_dv_fact", "race", "hiqual_dv_fact")
  pre_vars  <- c("sf12mcs_dv", "sf12pcs_dv")
  DT_base       <- data[wave == base_wave, c("pidp", base_vars), with = FALSE]
  DT_pre_window <- data[wave == pre_wave, c("pidp", pre_vars), with = FALSE]
  setnames(DT_base, base_vars, paste0(base_vars, "_base"))
  setnames(DT_pre_window, pre_vars, paste0(pre_vars, "_base"))

  data <- merge(DT_base, data, by = "pidp")
  data <- merge(DT_pre_window, data, by = "pidp")

  eligible_pidp <- intersect(
    data[wave == start_wave & response == 1, pidp],
    data[wave == final_wave & response == 1, pidp]
  )

  pop_DT <- data[pidp %in% eligible_pidp]

  data.table::setorder(pop_DT, pidp, wave)

  int_data_wide <- pop_DT[
      , `:=`(t0 = wave - as.integer(start_wave),
             pcs_lagged = shift(sf12pcs_dv, type = "lag"),
             mcs_lagged = shift(sf12mcs_dv, type = "lag"),
             dnc_lagged = shift(dnc_fact, type = "lag"),
             home_owner_lagged = shift(home_owner, type = "lag"),
             econ_benefits_lagged = shift(econ_benefits, type = "lag"),
             mastat_lagged = shift(mastat_dv, type = "lag"),
             econ_dist_bin_lagged = shift(econ_dist_bin, type = "lag")
            ),
        by = pidp
    ][wave %in% c(start_wave, final_wave)]

  outcome_col <- if (outcome == "MCS") "sf12mcs_dv" else "sf12pcs_dv"
  base_col    <- paste0(outcome_col, "_base")
  other_lag   <- if (outcome == "MCS") "pcs_lagged" else "mcs_lagged"

  # drop people missing the outcome scale at pre_wave
  final_data <- int_data_wide[!is.na(get(base_col))]

  # formula-safe factor levels before reshaping: downstream reformulate() breaks on some labels
  wide_data <- sanitize_factor_levels(final_data) |>
    make_wide(
      pidp,
      t0,
      base_cols = c(sex_dv_fact_base,
                    hiqual_dv_fact_base,
                    race_base,
                    !!rlang::sym(base_col)),
      outcome = !!rlang::sym(outcome_col),
      age_dv,
      gor_dv_fact,
      !!rlang::sym(other_lag),
      econ_dist_bin,
      econ_dist_bin_lagged,
      dnc_lagged,
      home_owner_lagged,
      econ_benefits_lagged,
      mastat_lagged,
      econ_emp_bin_fact,
      log_income,
      waves = 0:1
    ) |> data.table::as.data.table()

  # econ_dist_bin_lagged_1 is A0 (econ_dist_bin_0) itself
  wide_data[, econ_dist_bin_lagged_1 := NULL]
  data.table::setattr(wide_data, "time_lagged",
                      setdiff(attr(wide_data, "time_lagged"), "econ_dist_bin_lagged_1"))

  set_exposure(wide_data, exposure = econ_dist_bin)
}
