# Rubin's rules over the list of per-imputation tables from a
# map()/iteration = "list" target.
pool_ltmle <- function(ltmle_one) {
  if (is.data.frame(ltmle_one)) {
    stop("pool_ltmle() got a single per-imputation table, not a list of them -- ",
         "is the upstream target missing pattern = map(ltmle_imp_idx), iteration = \"list\"?",
         call. = FALSE)
  }
  data.table::rbindlist(ltmle_one)[, {
    p  <- mice::pool.scalar(Q = estimate, U = variance, n = Inf, k = 1)
    se <- sqrt(p$t)
    .(effect = p$qbar, se = se, ll = p$qbar - 1.96 * se, ul = p$qbar + 1.96 * se)
  }, by = .(estimand, type)]
}
