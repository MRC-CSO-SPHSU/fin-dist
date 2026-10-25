# Pool single-point TMLE estimates across imputations with Rubin's rules.
pool_tmle <- function(tmle_one) {
  # Expects the list of per-imputation fits produced by a map()/iteration="list"
  # target. A bare `tmle` object means the upstream target was not branched.
  if (inherits(tmle_one, "tmle")) {
    stop("pool_tmle() got a single `tmle` fit, not a list of per-imputation fits -- ",
         "is the upstream target missing pattern = map(tmle_imp_idx), iteration = \"list\"?",
         call. = FALSE)
  }

  estimands <- tibble::tribble(
    ~slot, ~estimand,     ~type,
    "EY0", "No distress", "mean",
    "EY1", "Distress",    "mean",
    "ATE", "ATE",         "ate"
  )

  purrr::pmap_dfr(estimands, function(slot, estimand, type) {
    Q <- purrr::map_dbl(tmle_one, \(fit) fit$estimates[[slot]]$psi)
    U <- purrr::map_dbl(tmle_one, \(fit) fit$estimates[[slot]]$var.psi)
    pooled <- mice::pool.scalar(Q = Q, U = U, n = Inf, k = 1)
    se <- sqrt(pooled$t)
    tibble::tibble(
      estimand = estimand,
      type     = type,
      effect   = pooled$qbar,
      se       = se,
      ll       = pooled$qbar - 1.96 * se,
      ul       = pooled$qbar + 1.96 * se
    )
  })
}
