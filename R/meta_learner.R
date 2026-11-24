# DR-learner for time-varying treaments over sex x race x education: score
# primitives, then estimate_cate(), which fits the nuisances and consumes them.

# Lower bound on the cumulative regime probability. Overridable per run through `gbound`.
CATE_EPS <- 0.01

bound_cum_g <- function(cum_g, eps = CATE_EPS) {
  pmin(pmax(cum_g, eps), 1)
}

# Sequential doubly-robust pseudo-outcome for one regime (a0, a1):
#   phi = Q1 + r0/pi0 * (Q2 - Q1) + r1/pi1 * (Y - Q2)
# r0 = 1{A0 = a0}, r1 = 1{A0 = a0, A1 = a1}; pi0 and pi1 are the (bounded) cumulative
# probabilities of following the regime through time 0 and through time 1.

seq_dr_pseudo_outcome <- function(r0, r1, Y, Q1, Q2, pi0, pi1) {
  args <- list(r0 = r0, r1 = r1, Y = Y, Q1 = Q1, Q2 = Q2, pi0 = pi0, pi1 = pi1)
  if (length(unique(lengths(args))) != 1L) {
    stop("seq_dr_pseudo_outcome: arguments differ in length: ",
         paste(names(args), lengths(args), sep = " = ", collapse = ", "))
  }
  purrr::iwalk(args, \(x, nm) if (anyNA(x)) {
    stop("seq_dr_pseudo_outcome: ", nm, " contains NA/NaN.")
  })
  if (any(pi0 <= 0 | pi0 > 1 | pi1 <= 0 | pi1 > 1)) {
    stop("seq_dr_pseudo_outcome: pi0 and pi1 must lie in (0, 1]; ",
         "bound them with bound_cum_g() first.")
  }
  if (any(r1 & !r0)) {
    stop("seq_dr_pseudo_outcome: r1 is set for rows where r0 is not.")
  }
  Q1 + r0 / pi0 * (Q2 - Q1) + r1 / pi1 * (Y - Q2)
}

# Stratum means of the columns of D (one estimand per column) with the
# influence-curve variance var(IC)/n, IC_i = 1{S_i = j}(D_ik - theta_jk) / p_j.
# Closed form: sum_j (D - theta_j)^2 * n / ((n - 1) * n_j^2).
stratum_estimates <- function(D, strata_id) {
  n <- nrow(D)
  stopifnot(length(strata_id) == n)
  long <- data.table::data.table(
    strata_id = rep(as.character(strata_id), times = ncol(D)),
    estimand  = factor(rep(colnames(D), each = n), levels = colnames(D)),
    d         = as.vector(D)
  )
  long[, .(estimate = mean(d),
           variance = sum((d - mean(d))^2) * n / ((n - 1) * .N^2)),
       keyby = .(strata_id, estimand)][, estimand := as.character(estimand)][]
}

# Per-stratum positivity diagnostics for the cumulative g of one regime, over the
# rows that followed it. An empty cell gets NA, not a warning from min()/median().
cum_g_diagnostics <- function(strata_id, strata_label, follow, pi1, pi1_raw, gbound) {
  data.table::data.table(strata_id = as.character(strata_id), strata_label, follow, pi1, pi1_raw)[
    , .(strata_label = strata_label[1L],
        n_j          = .N,
        n_follow     = sum(follow),
        pct_follow   = 100 * sum(follow) / .N,
        min_g        = if (any(follow)) min(pi1[follow]) else NA_real_,
        median_g     = if (any(follow)) stats::median(pi1[follow]) else NA_real_,
        share_g_at_bound = if (any(follow)) mean(pi1_raw[follow] <= gbound) else NA_real_),
    keyby = strata_id]
}

# Out-of-fold ensemble predictions from a SuperLearner fit. $Z holds the
# cross-validated library predictions and $coef the ensemble weights,
# How the two combine is method-specific: method.NNLS (gaussian fits) is the plain Z %*% coef,
# whereas method.NNloglik (binomial propensity) combines on the logit
# scale, plogis(trimLogit(Z) %*% coef).

sl_oof_predictions <- function(fit, n_expected, label) {
  if (is.null(fit$Z) || nrow(fit$Z) != n_expected) {
    stop("estimate_cate [", label, "]: nrow($Z) = ",
         if (is.null(fit$Z)) "NULL" else nrow(fit$Z),
         " but the training set had ", n_expected, " rows.")
  }
  if (length(fit$coef) != ncol(fit$Z)) {
    stop("estimate_cate [", label, "]: length($coef) = ", length(fit$coef),
         " != ncol($Z) = ", ncol(fit$Z), ".")
  }
  if (all(fit$coef == 0)) {
    stop("estimate_cate [", label, "]: every metalearner coefficient is zero; ",
         "computePred() would return all-zero predictions.")
  }
  oof <- drop(fit$method$computePred(fit$Z, fit$coef, control = fit$control))
  if (anyNA(oof)) {
    stop("estimate_cate [", label, "]: out-of-fold ensemble prediction contains NA.")
  }
  oof
}

# fit_ltmle_one()'s node split, repeated because editing that file re-runs every
# LTMLE branch. A test keeps the two equal.
cate_nodes <- function(outcome = c("MCS", "PCS")) {
  outcome <- match.arg(outcome)
  x <- confounders(outcome)
  list(W = grep("_1$", x, value = TRUE, invert = TRUE),
       L = grep("_1$", x, value = TRUE),
       A = c("econ_dist_bin_0", "econ_dist_bin_1"),
       Y = if (outcome == "MCS") "sf12mcs_dv_1" else "sf12pcs_dv_1")
}

# Honest predictions without outer folds: out-of-fold on the training rows, the
# full fit's prediction on the rest. Returns predictions only, so the fit (large
# with the real library) is freed on return.
sl_honest_predict <- function(y, X, train, family, sl_libs, sl_cv_control, label) {
  fit <- SuperLearner::SuperLearner(
    Y = y[train], X = X[train, , drop = FALSE], family = family,
    SL.library = sl_libs, cvControl = sl_cv_control,
    method = if (family$family == "binomial") "method.NNloglik" else "method.NNLS")
  pred <- numeric(length(y))
  pred[train] <- sl_oof_predictions(fit, sum(train), label)
  if (!all(train)) {
    pred[!train] <- as.numeric(stats::predict(fit, newdata = X[!train, , drop = FALSE],
                                              onlySL = TRUE)$pred)
  }
  pred
}

# Sequential DR-learner for one imputed dataset: per stratum and overall, the
# mean SF-12 under each static regime c(A0, A1) and its difference with the first
# (reference) regime.
estimate_cate <- function(wide_mids, imp_idx, sl_libs, outcome, regimes,
                          sl_cv_control, gbound = CATE_EPS) {
  stopifnot(
    "imp_idx must be a single imputation index: is the target missing pattern = map(ltmle_imp_idx)?" =
      length(imp_idx) == 1L,
    "gbound must be a single finite number in (0, 0.5)" =
      is.numeric(gbound) && length(gbound) == 1L && is.finite(gbound) &&
      gbound > 0 && gbound < 0.5,
    "regimes must be a list of c(A0, A1) 0/1 pairs with unique, non-empty names" =
      is.list(regimes) && rlang::is_named(regimes) && !anyDuplicated(names(regimes)) &&
      all(purrr::map_lgl(regimes, \(r) is.numeric(r) && length(r) == 2L && all(r %in% 0:1))),
    # SuperLearner would read an unnamed list by position
    "sl_cv_control must be a named list (V, stratifyCV, shuffle)" =
      is.list(sl_cv_control) && rlang::is_named(sl_cv_control)
  )
  # SuperLearner finds learners in the global environment. Batchtool workers finds them here with this function.
  register_sl_wrappers()

  nd  <- cate_nodes(outcome)
  dat <- as.data.frame(mice::complete(wide_mids, action = imp_idx))
  stopifnot("the completed data has NA in a node column" = !anyNA(dat[unlist(nd)]))
  n  <- nrow(dat)
  A0 <- as.integer(as.character(dat[[nd$A[1]]]))
  A1 <- as.integer(as.character(dat[[nd$A[2]]]))
  Y  <- dat[[nd$Y]]

  design <- \(cols) {
    X <- as.data.frame(stats::model.matrix(~ ., dat[cols])[, -1, drop = FALSE])
    names(X) <- make.names(names(X), unique = TRUE)
    X
  }
  XW  <- design(nd$W)
  XH  <- design(c(nd$W, nd$L))
  Xg1 <- data.frame(XH, econ_dist_bin_0 = A0)

  s <- make_strata(dat)
  strata_id    <- as.character(s$strata_id)
  strata_label <- as.character(s$strata_label)

  everyone <- rep(TRUE, n)
  g0hat <- sl_honest_predict(A0, XW,  everyone, stats::binomial(), sl_libs, sl_cv_control, "g0")
  g1hat <- sl_honest_predict(A1, Xg1, everyone, stats::binomial(), sl_libs, sl_cv_control, "g1")

  per_regime <- purrr::imap(regimes, \(a, k) {
    r0 <- A0 == a[1]
    r1 <- r0 & A1 == a[2]
    if (!any(r1)) {
      stop("estimate_cate: nobody follows regime '", k, "' in imputation ", imp_idx,
           ", so there is nothing to fit Q2 on.")
    }
    Q2 <- sl_honest_predict(Y,  XH, r1, stats::gaussian(), sl_libs, sl_cv_control, paste("Q2", k))
    Q1 <- sl_honest_predict(Q2, XW, r0, stats::gaussian(), sl_libs, sl_cv_control, paste("Q1", k))
    pi0_raw <- if (a[1] == 1) g0hat else 1 - g0hat
    pi1_raw <- pi0_raw * (if (a[2] == 1) g1hat else 1 - g1hat)
    pi1     <- bound_cum_g(pi1_raw, gbound)
    list(phi  = seq_dr_pseudo_outcome(r0, r1, Y, Q1, Q2, bound_cum_g(pi0_raw, gbound), pi1),
         diag = cum_g_diagnostics(strata_id, strata_label, r1, pi1, pi1_raw, gbound))
  })

  phi <- do.call(cbind, purrr::map(per_regime, "phi"))
  D   <- cbind(phi, phi[, -1, drop = FALSE] - phi[, 1])
  colnames(D) <- c(names(regimes),
                   paste(names(regimes)[-1], "-", names(regimes)[1], recycle0 = TRUE))

  estimates <- \(by) stratum_estimates(D, by)[, `:=`(
    imp  = as.integer(imp_idx),
    type = data.table::fifelse(estimand %in% names(regimes), "mean", "ate"))]
  labels   <- unique(data.table::data.table(strata_id, strata_label))
  gate     <- estimates(strata_id)[labels, on = "strata_id", strata_label := i.strata_label]
  overall  <- estimates(rep("overall", n))[, strata_id := NULL]
  diag_tbl <- data.table::rbindlist(purrr::map(per_regime, "diag"), idcol = "regime")[
    , imp := as.integer(imp_idx)]

  list(
    gate    = data.table::setcolorder(gate, c("imp", "strata_id", "strata_label", "estimand",
                                              "type", "estimate", "variance")),
    overall = data.table::setcolorder(overall, c("imp", "estimand", "type", "estimate", "variance")),
    diag    = data.table::setcolorder(diag_tbl, c("imp", "strata_id", "strata_label", "regime"))
  )
}
