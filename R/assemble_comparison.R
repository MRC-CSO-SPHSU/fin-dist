assemble_comparison <- function(ltmle_results, mi_results, mi_ate_results) {
  ltmle_of_type <- function(want) {
    ltmle_results |>
      dplyr::filter(type == want) |>
      dplyr::rename_with(~ paste0("tmle_", .x), c("effect", "se", "ll", "ul")) |>
      dplyr::select(!type)
  }

  comparison_marginal <- mi_results |>
    dplyr::relocate("intervention") |>
    dplyr::left_join(dplyr::rename(ltmle_of_type("mean"), intervention = estimand),
                     by = "intervention")

  # contrast rows are labelled "k - ref"; this drops the reference-mean (intercept) row
  comparison_ate <- mi_ate_results |>
    dplyr::filter(grepl(" - ", intervention, fixed = TRUE)) |>
    dplyr::rename(estimand = intervention) |>
    dplyr::relocate("estimand") |>
    dplyr::left_join(ltmle_of_type("ate"), by = "estimand")

  list(marginal = comparison_marginal, ate = comparison_ate)
}
