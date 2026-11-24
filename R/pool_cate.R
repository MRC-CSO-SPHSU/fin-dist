# Pool DR-learner results across imputations with Rubin's rules. Mirrors pool_ltmle(),
# which does the pooling; diag is a per-stratum overlap check, not a pooled estimate.
# A stratum can be absent from some imputations (imputed sex/race/hiqual can empty a cell),
# so each stratum is pooled over the imputations that have it.
pool_cate <- function(cate_one) {
  if (!is.null(cate_one$gate)) {
    stop("pool_cate() got a single estimate_cate() result, not a list of them -- ",
         "is the upstream target missing pattern = map(ltmle_imp_idx), iteration = \"list\"?",
         call. = FALSE)
  }
  # NA when every value is NA, without min()/max() warnings
  na_safe <- \(f) \(x) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)

  diag <- data.table::rbindlist(purrr::map(cate_one, "diag"))[, .(
    n_j_min          = min(n_j),
    n_j_max          = max(n_j),
    n_follow_min     = min(n_follow),
    pct_follow       = mean(pct_follow),
    min_g            = na_safe(min)(min_g),
    median_g         = na_safe(stats::median)(median_g),
    share_g_at_bound = na_safe(mean)(share_g_at_bound)
  ), by = .(strata_id, strata_label, regime)]

  list(
    gate    = pool_ltmle(purrr::map(cate_one, "gate"),
                         by = c("strata_id", "strata_label", "estimand", "type")),
    overall = pool_ltmle(purrr::map(cate_one, "overall")),
    diag    = diag,
    m       = length(cate_one)
  )
}
