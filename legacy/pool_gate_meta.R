# Pool the MAIHDA decomposition across imputations with Rubin's rules.
# SHRINK-THEN-POOL: the whole second stage runs inside each imputation and
# Rubin's rules are applied afterwards, because tau2 is a complete-data
# estimand. Pooling the GATEs first instead shrinks tau2 (averaged estimates
# have less spread, Rubin SEs are larger) and cuts PCV by roughly 4x --
# measured, not hypothetical (design 2026-07-31, 3).
#
# Mirrors pool_cate()/pool_tmle(): pool.scalar(n = Inf, k = 1), 1.96 cutoffs.

# tau2 is variance-scaled and bounded below by zero, so it is pooled on the log
# scale. eps keeps a boundary tau2 = 0 out of log(0); it is small against the
# observed tau2 range (0.70-5.37) so it does not bias the back-transform.
POOL_TAU2_EPS <- 1e-6

.pool_pair <- function(Q, U) {
  pooled <- mice::pool.scalar(Q = Q, U = U, n = Inf, k = 1)
  se <- sqrt(pooled$t)
  tibble::tibble(estimate = pooled$qbar, se = se,
                 ll = pooled$qbar - 1.96 * se, ul = pooled$qbar + 1.96 * se)
}

# Delta method on log(tau2 + eps): se_log = se_tau2 / (tau2 + eps).
#
# CAUTION: back-transforming a mean of logs is a GEOMETRIC mean, so a single
# imputation with tau2 = 0 contributes log(1e-6) = -13.8 and drags the pooled
# tau2 to near zero however large the others are. Measured on a 3-imputation
# fixture with one boundary hit: per-imputation tau2_main were 0.347, 0.101, 0;
# pooled tau2_main came out 0.0033 and PCV 0.998 (entirely additive). tau2 was
# interior in all 50 real imputations, so this does not bite today, but it
# would the moment the imputation model or m changes. Hence the warning in
# pool_gate_meta().
#
# tau2 is a variance (non-negative), so clamp both estimate and ll at zero.
# The ul can be Inf when every imputation hits the boundary: se_log = eps/eps
# is enormous and the interval is genuinely unbounded. This is the honest
# answer; the n_tau2_zero warning tells the user why.
.pool_tau2 <- function(tau2, se_tau2, eps = POOL_TAU2_EPS) {
  Q <- log(tau2 + eps)
  U <- (se_tau2 / (tau2 + eps))^2
  pooled <- mice::pool.scalar(Q = Q, U = U, n = Inf, k = 1)
  se <- sqrt(pooled$t)
  tibble::tibble(estimate = pmax(exp(pooled$qbar) - eps, 0),
                 ll = pmax(exp(pooled$qbar - 1.96 * se) - eps, 0),
                 ul = exp(pooled$qbar + 1.96 * se) - eps)
}

pool_gate_meta <- function(gate_meta_one) {
  m <- length(gate_meta_one)

  scalars_all <- purrr::map_dfr(gate_meta_one, "scalars")
  coefs_all   <- purrr::map_dfr(gate_meta_one, "coefs")
  cells_all   <- purrr::map_dfr(gate_meta_one, "cells")
  loco_all    <- purrr::map_dfr(gate_meta_one, "loco")

  # A boundary tau2 = 0 in ANY imputation collapses the geometric mean (see
  # .pool_tau2). Loud, because the failure direction is "no heterogeneity found"
  # -- the most dangerous way for a decomposition to be silently wrong. This
  # applies to BOTH layers that call .pool_tau2(): the 2*m full-model fits
  # (scalars_all) AND the 2*12*m leave-one-cell-out fits (loco_all). A real-data
  # check found n_tau2_zero == 0 but n_tau2_zero_loco == 7 of 600 -- the LOCO
  # layer hit the boundary while the full-model layer looked clean, which is
  # exactly why both must be counted and reported (see Fix 1b: the LOCO block
  # itself no longer uses log-scale pooling, but the count below still matters
  # as a documented, reviewable fact about the per-imputation fits).
  n_tau2_zero      <- sum(scalars_all$tau2_null == 0) + sum(scalars_all$tau2_main == 0)
  n_tau2_zero_loco <- sum(loco_all$tau2_null_d == 0) + sum(loco_all$tau2_main_d == 0)
  if (n_tau2_zero > 0 || n_tau2_zero_loco > 0) {
    hit_layers <- c(
      if (n_tau2_zero > 0)      "full-model" else NULL,
      if (n_tau2_zero_loco > 0) "leave-one-cell-out" else NULL
    )
    warning("pool_gate_meta: tau2 hit the zero boundary in the ",
            paste(hit_layers, collapse = " and "), " fit(s): ",
            n_tau2_zero, " of ", 2 * m, " full-model fits, ",
            n_tau2_zero_loco, " of ", 2 * 12 * m, " leave-one-cell-out fits. ",
            "The full-model tau2 is pooled on the log scale (a geometric mean), ",
            "so a boundary hit there pulls pooled tau2 toward zero and PCV ",
            "toward 1 -- inspect gate_meta_results$tests before reporting PCV. ",
            "The LOCO tau2 is pooled arithmetically (see the loco block) and so ",
            "is robust to a boundary hit there, but the count is reported for ",
            "transparency.")
  }

  t2_null <- .pool_tau2(scalars_all$tau2_null, scalars_all$se_tau2_null)
  t2_main <- .pool_tau2(scalars_all$tau2_main, scalars_all$se_tau2_main)
  mu_res  <- .pool_pair(scalars_all$mu, scalars_all$se_mu^2)

  # vt depends only on v_j and X, never on tau2, so averaging it is benign.
  vt_null_bar <- mean(scalars_all$vt_null)
  vt_main_bar <- mean(scalars_all$vt_main)

  # Extract scalar values before tibble constructor to avoid shadowing in
  # evaluation context (first assignment to 'mu' would rebind the variable).
  t2n_est <- t2_null$estimate[[1]]
  t2n_ll  <- t2_null$ll[[1]]
  t2n_ul  <- t2_null$ul[[1]]
  t2m_est <- t2_main$estimate[[1]]
  t2m_ll  <- t2_main$ll[[1]]
  t2m_ul  <- t2_main$ul[[1]]
  mu_est  <- mu_res$estimate[[1]]
  mu_se   <- mu_res$se[[1]]
  mu_ll   <- mu_res$ll[[1]]
  mu_ul   <- mu_res$ul[[1]]

  scalars <- tibble::tibble(
    tau2_null    = t2n_est,
    tau2_null_ll = t2n_ll, tau2_null_ul = t2n_ul,
    tau2_main    = t2m_est,
    tau2_main_ll = t2m_ll, tau2_main_ul = t2m_ul,
    # PCV from the POOLED tau2 values, never the mean of per-imputation PCVs:
    # individual imputations can give a negative PCV at J = 12 (REML lets the
    # additive model fit worse than the null). Reported as a point estimate
    # with no interval -- its sampling uncertainty is not credibly quantifiable
    # from 12 cells.
    pcv     = .pcv(t2n_est, t2m_est),
    I2_null = 100 * t2n_est / (vt_null_bar + t2n_est),
    I2_main = 100 * t2m_est / (vt_main_bar + t2m_est),
    mu      = mu_est, se_mu = mu_se, mu_ll = mu_ll, mu_ul = mu_ul
  )

  coefs <- coefs_all |>
    dplyr::group_by(term) |>
    dplyr::group_modify(function(d, key) .pool_pair(d$estimate, d$se^2)) |>
    dplyr::ungroup()

  pool_col <- function(d, est, se, prefix) {
    out <- .pool_pair(d[[est]], d[[se]]^2)
    stats::setNames(out, paste0(prefix, c("", "_se", "_ll", "_ul")))
  }

  cells <- cells_all |>
    dplyr::group_by(strata_id, strata_label) |>
    dplyr::group_modify(function(d, key) {
      dplyr::bind_cols(
        pool_col(d, "estimate",      "se",           "gate"),
        pool_col(d, "additive_pred", "se_pred",      "additive"),
        pool_col(d, "blup",          "se_blup",      "blup"),
        pool_col(d, "del_resid",     "se_del_resid", "resid"),
        tibble::tibble(
          n_j_min          = min(d$n_j),
          n_j_max          = max(d$n_j),
          pct_exposed      = mean(d$pct_exposed),
          min_g            = min(d$min_g),
          median_g         = stats::median(d$median_g),
          share_g_at_bound = mean(d$share_g_at_bound)
        )
      )
    }) |>
    dplyr::ungroup()

  loco <- loco_all |>
    dplyr::group_by(dropped_id, dropped_label) |>
    dplyr::group_modify(function(d, key) {
      # Arithmetic, NOT the log-scale geometric pooling used for the headline
      # tau2. tau2_*_d hits exactly 0 in real data (measured: 7 of 600 fits,
      # 6 of them on one cell), and log(0 + eps) = -13.8 would collapse the
      # geometric mean and inflate pcv_d -- 0.958 against a defensible 0.70.
      # This table is a descriptive sensitivity check with no interval, so the
      # arithmetic mean is both adequate and boundary-safe.
      #
      # NOT calling .pcv() here (unlike the scalars block above): tn can
      # legitimately be exactly 0 for a dropped cell, and .pcv() has no
      # zero-denominator guard -- it would return NaN where this block
      # deliberately returns NA_real_ instead.
      tn <- mean(d$tau2_null_d)
      tm <- mean(d$tau2_main_d)
      tibble::tibble(
        tau2_null_d = tn,
        tau2_main_d = tm,
        pcv_d       = if (tn == 0) NA_real_ else (tn - tm) / tn,
        # per-imputation median, as a robustness column. Guard each element's
        # division the same way as pcv_d: a per-imputation tau2_null_d == 0
        # is undefined (0/0), not zero or infinite, so it becomes NA and is
        # dropped from the median rather than corrupting it with NaN/Inf.
        pcv_d_median = stats::median(
          ifelse(d$tau2_null_d == 0, NA_real_,
                 (d$tau2_null_d - d$tau2_main_d) / d$tau2_null_d),
          na.rm = TRUE
        ),
        n_zero_d = sum(d$tau2_null_d == 0) + sum(d$tau2_main_d == 0)
      )
    }) |>
    dplyr::ungroup()

  # Q and F statistics do not Rubin-pool. Return the per-imputation rows and let
  # the report summarise median and range -- the convention pool_cate() already
  # uses for its Wald test.
  tests <- scalars_all |>
    dplyr::select(imp, QE_null, QE_main, QM, QMp, pcv)

  list(scalars = scalars, coefs = coefs, cells = cells,
       loco = loco, tests = tests, m = m, n_tau2_zero = n_tau2_zero,
       n_tau2_zero_loco = n_tau2_zero_loco)
}
