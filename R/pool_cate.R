# Pool DR-learner GATE results across imputations with Rubin's rules. Mirrors pool_tmle()
# Stratum membership is computed on completed data and can in principle vary across imputations

pool_cate <- function(cate_one) {
  m <- length(cate_one)

  gate_all <- purrr::map_dfr(cate_one, "gate")
  blp_all  <- purrr::map_dfr(cate_one, "blp")
  ate_all  <- purrr::map_dfr(cate_one, "ate")
  wald_all <- purrr::map_dfr(cate_one, "wald")

  pool_pair <- function(Q, U) {
    pooled <- mice::pool.scalar(Q = Q, U = U, n = Inf, k = 1)
    se <- sqrt(pooled$t)
    tibble::tibble(
      estimate = pooled$qbar,
      se       = se,
      ll       = pooled$qbar - 1.96 * se,
      ul       = pooled$qbar + 1.96 * se
    )
  }

  gate <- gate_all |>
    dplyr::group_by(strata_id, strata_label) |>
    dplyr::group_modify(function(g, key) {
      dplyr::bind_cols(
        pool_pair(g$estimate, g$se^2),
        tibble::tibble(
          n_j_min          = min(g$n_j),
          n_j_max          = max(g$n_j),
          pct_exposed      = mean(g$pct_exposed),
          min_g            = min(g$min_g),
          median_g         = stats::median(g$median_g),
          share_g_at_bound = mean(g$share_g_at_bound)
        )
      )
    }) |>
    dplyr::ungroup()

  blp <- blp_all |>
    dplyr::group_by(term) |>
    dplyr::group_modify(function(g, key) pool_pair(g$estimate, g$se^2)) |>
    dplyr::ungroup()

  ate <- pool_pair(ate_all$estimate, ate_all$se^2)

  list(gate = gate, blp = blp, ate = ate, wald = wald_all, m = m)
}
