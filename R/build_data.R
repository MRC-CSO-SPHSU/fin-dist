build_data <- function(data, pre_wave = 2, target_wave = 3, base_wave = 1,
                       outcome = c("MCS", "PCS")) {
  library(dplyr)

  outcome <- match.arg(outcome)

  # Baseline data generation: Baseline will change depending of the analytic window we choose
  base_cols <- c("pidp", "wave", "sex_dv_fact", "race", "hiqual_dv_fact", "sf12mcs_dv", "sf12pcs_dv")
  DT_base <- data[wave == base_wave, ..base_cols]

  setnames(DT_base, setdiff(names(DT_base), c("pidp", "wave")),
           paste0(setdiff(names(DT_base), c("pidp", "wave")), "_base"))

  # Joining baseline DT with long DT
  data <- merge(DT_base[, !"wave"], data, by = c("pidp"))

  eligible_pidp <- intersect(
    data[wave == pre_wave & response == 1, pidp],
    data[wave == target_wave & response == 1, pidp]
  )

  pop_DT <- data[pidp %in% eligible_pidp]

  data.table::setorder(pop_DT, pidp, wave)

  int_data_wide <- pop_DT[
      , `:=`(t0 = wave - as.integer(target_wave),
             pcs_lagged = shift(sf12pcs_dv, type = "lag"),
             mcs_lagged = shift(sf12mcs_dv, type = "lag"),
             dnc_lagged = shift(dnc_fact, type = "lag"),
             home_owner_lagged = shift(home_owner, type = "lag"),
             econ_benefits_lagged = shift(econ_benefits, type = "lag"),
             mastat_lagged = shift(mastat_dv, type = "lag"),
             econ_dist_bin_lagged = shift(econ_dist_bin, type = "lag")
            ),
        by = pidp
    ][wave %in% target_wave]

  # Drop pidps with no usable MCS/PCS: baseline missing
  outcome_col <- if (outcome == "MCS") "sf12mcs_dv" else "sf12pcs_dv"
  base_col    <- paste0(outcome_col, "_base")

  bad_pidps <- int_data_wide[
    , .(base_na = any(is.na(get(base_col)))),
    by = pidp
  ][(base_na), pidp]

  final_data <- int_data_wide[!pidp %in% bad_pidps]
  
  # Formula-safe factor levels before reshaping. 
  # Downstream code uses reformulate, which may conflict with some var labels
  final_data <- final_data |>
    sanitize_factor_levels()

  if (outcome == "MCS") {
  wide_data <- final_data |>
    make_wide(
      pidp,
      t0,
      base_cols = c(sex_dv_fact_base,
                    hiqual_dv_fact_base,
                    race_base,
                    sf12mcs_dv_base),
      outcome = sf12mcs_dv,
      age_dv,
      gor_dv_fact,
      pcs_lagged,
      econ_dist_bin,
      econ_dist_bin_lagged,
      dnc_lagged,
      home_owner_lagged,
      econ_benefits_lagged,
      mastat_lagged,
      econ_emp_bin_fact,
      log_income,
      waves = 0L
    ) |> data.table::as.data.table()
  } else if (outcome == "PCS") {
    wide_data <- final_data |>
      make_wide(
        pidp,
        t0,
        base_cols = c(sex_dv_fact_base,
                      hiqual_dv_fact_base,
                      race_base,
                      sf12pcs_dv_base),
        outcome = sf12pcs_dv,
        age_dv,
        gor_dv_fact,
        mcs_lagged,
        econ_dist_bin,
        econ_dist_bin_lagged,
        dnc_lagged,
        home_owner_lagged,
        econ_benefits_lagged,
        mastat_lagged,
        econ_emp_bin_fact,
        log_income,
        waves = 0L
      ) |> data.table::as.data.table()
  }
  
  wide_data <- set_exposure(wide_data, 
                            exposure = econ_dist_bin)
  
  return(wide_data)
}
