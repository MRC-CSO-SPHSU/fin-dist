assemble_comparison <- function(tmle_results, mi_results, mi_ate_results) {
  tmle_marginal_means <- tmle_results |>
    dplyr::rename(intervention = estimand) |>
    dplyr::filter(type == "mean") |>
    dplyr::rename_with(~ paste0("tmle_", .x), c("effect", "se", "ll", "ul")) |>
    dplyr::select(!type)

  comparison_marginal <- mi_results |>
    dplyr::mutate(intervention = dplyr::case_when(intervention == 0 ~ "No distress",
                                                  intervention == 1 ~ "Distress")) |>
    dplyr::relocate("intervention") |>
    dplyr::left_join(tmle_marginal_means, by = "intervention")

  tmle_ate_results <- tmle_results |>
    dplyr::filter(estimand == "ATE") |>
    dplyr::mutate(estimand = dplyr::recode_values(estimand, "ATE" ~ "Distress")) |>
    dplyr::rename_with(~ paste0("tmle_", .x), c("effect", "se", "ll", "ul")) |>
    dplyr::select(!type)

  comparison_ate <- mi_ate_results |>
    dplyr::filter(intervention == 1) |>
    dplyr::mutate(intervention = dplyr::recode_values(intervention, 1 ~ "Distress")) |>
    dplyr::rename(estimand = intervention) |>
    dplyr::relocate("estimand") |>
    dplyr::left_join(tmle_ate_results, by = "estimand")

  return(list(
    marginal = comparison_marginal,
    ate = comparison_ate
  ))
}