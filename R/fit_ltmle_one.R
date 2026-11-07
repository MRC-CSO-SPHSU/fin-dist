# Two-time-point LTMLE (ltmle::ltmleMSM) for one imputed dataset: mean SF-12
# under each static distress regime c(A0, A1), and each regime's contrast with
# the first (reference) regime.
#
# One ltmleMSM() call, not one ltmle() call per regime, so the g-models are
# fitted once and shared. The MSM is saturated (one indicator per regime, equal
# weights), so its estimates equal those of separate ltmle() fits. Returns only
# a compact table: the fit object is ~100 MB at n = 20k.
fit_ltmle_one <- function(wide_mids, imp_idx, sl_libs, outcome, regimes,
                          sl_cv_control, gbound) {
  # mice::complete() with a vector `action` silently stacks the imputations, so
  # an unbranched target would fit one TMLE on m * n rows instead of erroring.
  stopifnot(
    "imp_idx must be a single imputation index -- is the target missing pattern = map(ltmle_imp_idx)?" =
      length(imp_idx) == 1L
  )
  # ltmle checks only the names it is given, so an unnamed list gets through and
  # SuperLearner then reads it by position.
  stopifnot(
    "sl_cv_control must be a named list (V, stratifyCV, shuffle)" =
      is.list(sl_cv_control) && !is.null(names(sl_cv_control))
  )
  # SuperLearner finds learners by name in globalenv(); fresh batchtools workers
  # do not have the custom wrappers there.
  register_sl_wrappers()

  # W, all measured pre-exposure. Shared with estimate_cate() via
  # R/confounders.R: do not inline a copy here.
  confounders <- confounders(outcome)
  W <- grep("_1$", confounders, value = TRUE, invert = TRUE)
  L <- grep("_1$", confounders, value = TRUE)
  A <- c("econ_dist_bin_0", "econ_dist_bin_1")
  Y <- if (outcome == "MCS") "sf12mcs_dv_1" else "sf12pcs_dv_1"

  dat <- data.table::as.data.table(mice::complete(wide_mids, action = imp_idx))
  dat[, (A) := purrr::map(.SD, \(a) as.integer(as.character(a))), .SDcols = A]
  # ltmle reads the node order from the column order: W, A0, L1, A1, Y.
  d <- as.data.frame(dat[, c(W, A[1], L, A[2], Y), with = FALSE])
  n <- nrow(d)
  K <- length(regimes)

  fit <- suppressMessages(ltmle::ltmleMSM(
    data             = d,
    Anodes           = A,
    Lnodes           = L,
    Ynodes           = Y,
    regimes          = simplify2array(purrr::map(regimes, \(a) matrix(a, n, length(a), byrow = TRUE))),
    working.msm      = paste("Y ~ -1 +", paste(names(regimes), collapse = " + ")),
    summary.measures = array(diag(K), c(K, K, 1),
                             dimnames = list(names(regimes), names(regimes), NULL)),
    msm.weights      = NULL,
    gbounds          = gbound,
    SL.library       = sl_libs,
    SL.cvControl     = sl_cv_control,
    variance.method  = "ic",
    estimate.time    = FALSE
  ))
  ltmle_msm_estimates(fit, regimes, imp_idx)
}

# Regime means and contrasts vs the first regime, on the SF-12 scale, from a
# saturated ltmleMSM fit. fit$beta is logit(mean) on the [0,1] scale ltmle
# rescales Y to; the delta method carries fit$IC back to the SF-12 scale.
ltmle_msm_estimates <- function(fit, regimes, imp_idx) {
  K  <- length(regimes)
  # FALSE when Y was already in [0,1] and ltmle did not rescale it.
  yr <- if (isFALSE(fit$transformOutcome)) c(0, 1) else attr(fit$transformOutcome, "Yrange")
  m   <- stats::plogis(unname(fit$beta))
  psi <- m * diff(yr) + yr[1]
  IC  <- sweep(fit$IC, 2, m * (1 - m) * diff(yr), `*`)
  D   <- cbind(IC, IC[, -1, drop = FALSE] - IC[, 1])
  data.table::data.table(
    imp      = as.integer(imp_idx),
    estimand = c(names(regimes), paste(names(regimes)[-1], "-", names(regimes)[1])),
    type     = rep(c("mean", "ate"), c(K, K - 1L)),
    estimate = c(psi, psi[-1] - psi[1]),
    # var(IC) / n, as summary.ltmle() computes it.
    variance = unname(diag(stats::cov(D))) / nrow(D)
  )
}
