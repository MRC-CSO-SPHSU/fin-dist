# g-formula MI predictor matrix and imputation of counterfactual outcomes

run_gformula <- function(wide_mids, wide_data_mi, estimand, outcome_scale, M = 50) {
  estimand <- as.character(estimand)

  regimes <- list(0,1)

  predictor_matrix <- make_counterfactual_matrix(wide_data_mi)

  outcome_final <- if (outcome_scale == "MCS") "sf12mcs_dv_0" else "sf12pcs_dv_0"

  method <- make_counterfactual_method(wide_data_mi)

  imps <- gFormulaMI::gFormulaImpute(
    data             = wide_mids,
    M                = M,
    trtVars          = c("econ_dist_bin_0"),
    trtRegimes       = regimes,
    predictorMatrix  = predictor_matrix,
    method           = method,
    silent           = TRUE
  )

  fits <- imps %$%
    lm(as.formula(paste(outcome_final, "~", estimand)))

  outvals <- gFormulaMI::syntheticPool(fits)

  regimes_1 <- tibble::tibble(
    intervention = regimes |>
      purrr::reduce(c)
  )

  out <- outvals |>
    tibble::as_tibble() |>
    tibble::rownames_to_column("Intervention") |>
    dplyr::transmute(
      mi_effect = Estimate,
      mi_se     = sqrt(Total),
      mi_ll     = `95% CI L`,
      mi_ul     = `95% CI U`
    ) |>
    dplyr::bind_cols(regimes_1)

  return(out)
}
