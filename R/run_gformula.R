# g-formula by multiple imputation (gFormulaMI): impute synthetic data under each
# regime, fit `outcome ~ estimand` in each dataset and pool with syntheticPool().
run_gformula <- function(wide_mids, wide_data_mi, estimand, outcome_scale, regimes, M = 50) {
  outcome <- attr(wide_data_mi, "outcome_var")
  stopifnot(startsWith(outcome, c(MCS = "sf12mcs", PCS = "sf12pcs")[[outcome_scale]]))
  predictor_matrix <- make_counterfactual_matrix(wide_data_mi)
  method <- make_counterfactual_method(wide_data_mi)

  imps <- gFormulaMI::gFormulaImpute(
    data             = wide_mids,
    M                = M,
    trtVars          = attr(wide_data_mi, "exposure_vars"),
    trtRegimes       = unname(regimes),
    predictorMatrix  = predictor_matrix,
    method           = method,
    silent           = TRUE
  )
  fits    <- with(imps, lm(as.formula(paste(outcome, "~", estimand))))
  outvals <- gFormulaMI::syntheticPool(fits)

  # gFormulaImpute codes regime k (list order) as factor level k: row
  # "factor(regime)k" is regime k's mean, or with an intercept its contrast with
  # regime 1, whose mean is then "(Intercept)"
  rn <- rownames(outvals)
  k  <- match(rn, paste0("factor(regime)", seq_along(regimes)))
  k[rn == "(Intercept)"] <- 1L
  stopifnot(!anyNA(k))
  is_contrast  <- "(Intercept)" %in% rn & rn != "(Intercept)"
  intervention <- ifelse(is_contrast,
                         paste(names(regimes)[k], "-", names(regimes)[1]),
                         names(regimes)[k])

  outvals |>
    tibble::as_tibble() |>
    dplyr::transmute(
      mi_effect = Estimate,
      mi_se     = sqrt(Total),
      mi_ll     = `95% CI L`,
      mi_ul     = `95% CI U`,
      intervention = intervention
    )
}
