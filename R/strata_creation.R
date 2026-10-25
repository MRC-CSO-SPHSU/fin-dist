# Intersectional strata (sex x race x education) for the DR-learner GATE arm.
make_strata <- function(df) {
  needed <- c("sex_dv_fact_base", "race_base", "hiqual_dv_fact_base")
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols)) {
    stop("make_strata: missing stratum column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  dplyr::mutate(df,
    # intersectional strata id: unique cell
    strata_id = 100 * as.numeric(sex_dv_fact_base) +
                 10 * as.numeric(race_base) +
                      as.numeric(hiqual_dv_fact_base),
    strata_id = as.factor(strata_id),
    # human-readable stratum name, e.g. "Female:White:High"
    strata_label = interaction(sex_dv_fact_base, race_base, hiqual_dv_fact_base,
                               sep = ":", drop = TRUE)
  )
}
