# Shared confounder set (W) for the marginal TMLE and the DR-learner GATE arm.

confounders <- function(outcome = c("MCS", "PCS")) {
  outcome <- match.arg(outcome)
  c("sex_dv_fact_base",
    "hiqual_dv_fact_base",
    "race_base",
    if (outcome == "MCS") "sf12mcs_dv_base" else "sf12pcs_dv_base",
    "age_dv_0",
    "gor_dv_fact_0",
    "econ_dist_bin_lagged_0",
    if (outcome == "MCS") "pcs_lagged_0" else "mcs_lagged_0",
    "dnc_lagged_0",
    "home_owner_lagged_0",
    "econ_benefits_lagged_0",
    "mastat_lagged_0",
    "econ_emp_bin_fact_0",
    "log_income_0")
}
