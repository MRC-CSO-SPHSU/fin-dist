# DR-learner CATEs over sex x race x education -- score primitives.
# estimate_cate() below (added with the full rewrite) consumes these.

# Default bound for clever variable (propensity) positivity. Overridable per run
# through the `gbound` argument of estimate_cate() (set from _targets.R).
CATE_EPS <- 0.01

# Function to apply the bound over IPTW scores
bound_propensity <- function(g_hat, eps = CATE_EPS) {
  pmin(pmax(g_hat, eps), 1 - eps)
}

# Function to compute the doubly-robust pseudo-outcome for the DR-learner CATE arm.
dr_pseudo_outcome <- function(A, Y, mu0_hat, mu1_hat, g_hat) {
  stopifnot(length(A) == length(Y),
            length(mu0_hat) == length(Y),
            length(mu1_hat) == length(Y),
            length(g_hat) == length(Y))
  if (anyNA(g_hat)) {
    stop("dr_pseudo_outcome: g_hat contains NA/NaN -- a degenerate ",
         "SuperLearner propensity fit. Cannot bound or use it.")
  }
  if (any(g_hat <= 0 | g_hat >= 1)) {
    stop("dr_pseudo_outcome: g_hat must be strictly inside (0, 1); ",
         "bound it with bound_propensity() first.")
  }
  (mu0_hat - mu1_hat) +
    (1 - A) / (1 - g_hat) * (Y - mu0_hat) -
     A      /      g_hat  * (Y - mu1_hat)
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

# DR-learner CATEs over sex x race x education, one imputed dataset. 
# Pooled across imputations by pool_cate().

estimate_cate <- function(wide_mids, imp_idx, sl_libs, outcome,
                          gbound = CATE_EPS) {

  imp_idx <- as.integer(imp_idx)

  # Symmetric bound: g_hat is clipped to [gbound, 1 - gbound], so anything at or
  # above 0.5 would collapse every propensity onto a single point.
  stopifnot(
    "gbound must be a single finite number in (0, 0.5)" =
      is.numeric(gbound) && length(gbound) == 1L &&
      is.finite(gbound) && gbound > 0 && gbound < 0.5
  )

  # SuperLearner resolves learner names by string against globalenv() on fresh
  # batchtools workers 
  register_sl_wrappers()

  ## -- 3.1 complete data call -----------------------------------------------
  dat <- mice::complete(wide_mids, action = imp_idx)

  ## -- 3.2 confounders helper call (see helpers.R) ---------------
  W <- confounders(outcome)

  ## -- 3.3 numeric design matrix (SL.glmnet/SL.xgboost cannot take factors) --
  Wmat <- as.data.frame(model.matrix(~ ., dat[W])[, -1])
  names(Wmat) <- make.names(names(Wmat), unique = TRUE)

  stopifnot(nrow(Wmat) == nrow(dat))

  ## -- 3.4 treatment and outcome ---------------------------------------------
  A <- as.integer(as.character(dat$econ_dist_bin_0))  # refactoring as it needs to be an integer for SuperLearner
  Y <- dat[[if (outcome == "MCS") "sf12mcs_dv_0" else "sf12pcs_dv_0"]]
  n <- nrow(dat)

  ## -- 3.5 main parameters for outcome regression and propensity score modeling --------

  mu0 <- SuperLearner::SuperLearner(
    Y = Y[A == 0], X = Wmat[A == 0, , drop = FALSE], family = gaussian(),
    SL.library = sl_libs, cvControl = list(V = 5))
  mu1 <- SuperLearner::SuperLearner(
    Y = Y[A == 1], X = Wmat[A == 1, , drop = FALSE], family = gaussian(),
    SL.library = sl_libs, cvControl = list(V = 5))
  gm  <- SuperLearner::SuperLearner(
    Y = A, X = Wmat, family = binomial(),
    SL.library = sl_libs, method = "method.NNloglik",
    cvControl = list(V = 5))

  ## -- 3.6 honest predictions without outer folds ----------------------------

  mu0_hat <- numeric(n)
  mu1_hat <- numeric(n)
  mu0_hat[A == 0] <- sl_oof_predictions(mu0, sum(A == 0), "mu0")
  mu0_hat[A == 1] <- predict(mu0, newdata = Wmat[A == 1, , drop = FALSE])$pred
  mu1_hat[A == 1] <- sl_oof_predictions(mu1, sum(A == 1), "mu1")
  mu1_hat[A == 0] <- predict(mu1, newdata = Wmat[A == 0, , drop = FALSE])$pred

  stopifnot(!anyNA(mu0_hat), !anyNA(mu1_hat), !anyNA(Y))
  g_hat <- sl_oof_predictions(gm, n, "g")

  ## -- 3.7 positivity bound and pseudo-outcome --
  g_hat <- bound_propensity(g_hat, eps = gbound)
  psi   <- dr_pseudo_outcome(A, Y, mu0_hat, mu1_hat, g_hat)

  ## -- 3.8 strata creation -----------------------------------
  dat_s   <- make_strata(dat)
  proj_df <- data.frame(psi          = psi,
                        strata_id    = dat_s$strata_id,
                        strata_label = dat_s$strata_label)

  proj  <- estimatr::lm_robust(psi ~ 0 + strata_id, data = proj_df,
                               se_type = "HC2")
  cells <- names(coef(proj))
  lh    <- car::linearHypothesis(proj, paste0(cells[-1], " = ", cells[1]),
                                 test = "F")
  wald  <- tibble::tibble(imp = imp_idx,
                          F   = lh$F[2],
                          df1 = lh$Df[2],
                          df2 = lh$Res.Df[2],
                          p   = lh[2, "Pr(>F)"])

  ## -- 3.9 BLP global heterogeneity test (conventional Y(1)-Y(0) scale) ------
  m_hat <- g_hat * mu1_hat + (1 - g_hat) * mu0_hat  
  tau   <- -proj$fitted.values                      # undo the 3.7 flip, once
  blp_fit <- lm(I(Y - m_hat) ~ 0 + I(A - g_hat) +
                  I((tau - mean(tau)) * (A - g_hat)))
  blp_ct  <- lmtest::coeftest(blp_fit, vcov = sandwich::vcovHC, type = "HC3")
  blp <- tibble::tibble(
    imp      = imp_idx,
    term     = c("beta1_ate", "beta2_hte"),
    estimate = unname(blp_ct[, "Estimate"]),
    se       = unname(blp_ct[, "Std. Error"]),
    p        = unname(blp_ct[, "Pr(>|t|)"])
  )

  ## -- 3.10 per-cell diagnostics ---------------------------------------------
  diag_tbl <- proj_df |>
    dplyr::mutate(A = A, g = g_hat) |>
    dplyr::group_by(strata_id, strata_label) |>
    dplyr::summarise(
      n_j              = dplyr::n(),
      n_exposed        = sum(A),
      pct_exposed      = 100 * mean(A),
      min_g            = min(g),
      median_g         = stats::median(g),
      share_g_at_bound = mean(g <= gbound | g >= 1 - gbound),
      .groups          = "drop"
    ) |>
    dplyr::mutate(strata_id = as.character(strata_id))

  est_tbl <- tibble::tibble(
    strata_id = sub("^strata_id", "", cells),
    estimate  = unname(coef(proj)),
    se        = unname(proj$std.error)
  )

  gate <- diag_tbl |>
    dplyr::inner_join(est_tbl, by = "strata_id") |>
    dplyr::mutate(imp = imp_idx, .before = 1)

  if (nrow(gate) != length(cells)) {
    stop("estimate_cate: projection produced ", length(cells),
         " cells but the diagnostics join kept ", nrow(gate),
         " -- cell ids disagree between coef names and make_strata().")
  }

  ## -- 3.11 output ------------------------------------------------------------
  # Scale note: ate$estimate is mean(psi), the Y(0)-Y(1) scale;
  # blp$estimate[blp$term == "beta1_ate"] is the same population ATE but on
  # the conventional Y(1)-Y(0) scale. The two therefore always print with opposite signs. 
  # This is by design, not a bug.
  list(
    gate = gate,
    blp  = blp,
    ate  = tibble::tibble(imp      = imp_idx,
                          estimate = mean(psi),
                          se       = stats::sd(psi) / sqrt(n)),
    wald = wald
  )
}
