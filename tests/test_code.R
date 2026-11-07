# Tests for the pipeline functions in R/.
here::i_am("tests/test_code.R")

## reproducing targets environment
pacman::p_load(testthat,
  data.table,
  dplyr,
  tidyr,
  tibble,
  purrr,
  rlang,
  here,
  mice,
  lattice,
  SuperLearner,
  estimatr,
  car,
  lmtest,
  sandwich,
  ltmle,
  gFormulaMI,
  magrittr)

## shared constants. Defined up here because test_file(desc = ) runs only the
## code above the matching test.
# the 27 MCS wide_data columns, in build_data() order: _base, then _0
# (start_wave) and _1 (final_wave)
wide_cols_mcs <- c(
  "sex_dv_fact_base", "hiqual_dv_fact_base", "race_base", "sf12mcs_dv_base",
  "age_dv_0", "gor_dv_fact_0", "pcs_lagged_0", "econ_dist_bin_0",
  "econ_dist_bin_lagged_0", "dnc_lagged_0", "home_owner_lagged_0",
  "econ_benefits_lagged_0", "mastat_lagged_0", "econ_emp_bin_fact_0",
  "log_income_0", "sf12mcs_dv_0",
  "age_dv_1", "gor_dv_fact_1", "pcs_lagged_1", "econ_dist_bin_1",
  "dnc_lagged_1", "home_owner_lagged_1", "econ_benefits_lagged_1",
  "mastat_lagged_1", "econ_emp_bin_fact_1", "log_income_1", "sf12mcs_dv_1"
)
# ltmle_regimes and sl_cv_control from _targets.R: c(A0, A1), first = reference
toy_regimes <- list(never = c(0L, 0L), before = c(1L, 0L),
                    now = c(0L, 1L), permanent = c(1L, 1L))
toy_cv <- list(V = 5L, stratifyCV = FALSE, shuffle = TRUE)
# fit_ltmle_one()/pool_ltmle() rows: 4 regime means, then 3 contrasts
ltmle_estimands <- c("never", "before", "now", "permanent",
                     "before - never", "now - never", "permanent - never")

### build_data: the two-wave wide columns
test_that("build_data produces expected output", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))
  source(here::here("R", "confounders.R"))

  pop_data <- import_data(force = FALSE) |> clean_data() |> preproc_data()

  wide_data <- build_data(pop_data, pre_wave = 1, base_wave = 2, start_wave = 3,
                          final_wave = 4, outcome = "MCS")
  expect_true(is.data.frame(wide_data))

  # exact names and order: ltmle reads the node order from the columns
  expect_identical(names(wide_data), wide_cols_mcs)
  expect_equal(attr(wide_data, "outcome_var"), "sf12mcs_dv_1")
  expect_equal(attr(wide_data, "exposure_vars"), c("econ_dist_bin_0", "econ_dist_bin_1"))
  # econ_dist_bin_lagged_1 is A0 itself: dropped, from the attribute too
  expect_false("econ_dist_bin_lagged_1" %in% names(wide_data))
  expect_false("econ_dist_bin_lagged_1" %in% attr(wide_data, "time_lagged"))
  expect_true(all(c(confounders("MCS"), time_varying_confounders("MCS")) %in%
                    names(wide_data)))
})

### make_counterfactual_method: the methods gFormulaImpute() draws the
### synthetic block with. Synthetic frame: no raw .dta needed.
test_that("make_counterfactual_method covers every column and keeps factors categorical", {
  source(here::here("R", "helpers.R"))

  wide <- data.frame(
    sex_dv_base         = factor(c("Female", "Male")),                    # complete
    hiqual_dv_fact_base = factor(c("High", "Medium"),
                                 levels = c("High", "Medium", "Low")),    # >2 levels
    sf12mcs_dv_base     = c(40, 60),                                      # complete numeric
    gor_dv_fact_0       = factor(c("R1", "R2"), levels = paste0("R", 1:12)),
    econ_dist_bin_0     = factor(c("0", "1")),                            # exposure
    pcs_lagged_0        = c(45, NA),
    sf12mcs_dv_0        = c(NA, 55)
  )
  wide <- set_exposure(wide, exposure = econ_dist_bin)

  method <- make_counterfactual_method(wide)

  # gFormulaImpute() appends regime = "" and hands this straight to mice, whose
  # check.method() indexes by name and blanks by position: one entry per column,
  # in column order.
  expect_equal(names(method), names(wide))

  # treatment is set from the regime, never imputed
  expect_equal(unname(method["econ_dist_bin_0"]), "")

  # every other column must be drawable: the synthetic block is all NA, so a ""
  # here would leave the column missing in the counterfactual population. This
  # is what the all-NA probe row inside the function buys.
  expect_false(any(method[setdiff(names(wide), "econ_dist_bin_0")] == ""))

  # pmm, not gFormulaMI's default unbounded "norm", for the continuous nodes
  expect_equal(unname(method[c("sf12mcs_dv_0", "pcs_lagged_0")]), c("pmm", "pmm"))

  # >2-level factors must NOT be pmm: pmm regresses the integer level codes and,
  # under gFormulaImpute's maxit = 1, collapses the synthetic marginal onto a
  # few donors (levels have gone missing entirely).
  expect_equal(unname(method[c("hiqual_dv_fact_base", "gor_dv_fact_0")]),
               c("polyreg", "polyreg"))
  expect_equal(unname(method["sex_dv_base"]), "logreg")

  # a wide frame that never went through set_exposure() fails loudly rather than
  # silently letting the treatment column be imputed from a model
  bare <- wide
  attr(bare, "exposure_vars") <- NULL
  expect_error(make_counterfactual_method(bare), "exposure_vars")
})

### run_mice and make_counterfactual_matrix on the real data
test_that("run_mice produces expected output", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))

  pop_data  <- import_data(force = FALSE) |> clean_data() |> preproc_data()
  wide_data <- build_data(pop_data, pre_wave = 1, base_wave = 2, start_wave = 3,
                          final_wave = 4, outcome = "MCS")

  wide_data <- set_exposure(wide_data, exposure = "econ_dist_bin")
  predictor_matrix <- make_counterfactual_matrix(wide_data)

  # mini-mice call to check structure only
  mids <- run_mice(wide_data, m = 2, maxit = 2, seed = 42, outcome = "MCS")

  expect_s3_class(mids, "mids")
  expect_equal(mids$m, 2L)

  # mice diagnostics: only trust them with at least 15 iterations
  source(here::here("tests", "diagnose_mice.R"))
  diagnose_mice(mids, wide_data,
                plot_dir = here::here("tests", "figs", "wide", "mice_diag"))

  # no exposure/outcome interaction dropped. loggedEvents need not be empty:
  # mice logs age_dv_1 (~ age_dv_0 + 1) as "collinear", which is harmless.
  expect_true(assert_congenial_terms(mids, a_col = c("econ_dist_bin_0", "econ_dist_bin_1"),
                                     y_col = "sf12mcs_dv_1"))

  expect_true("regime" %in% rownames(predictor_matrix))
  off_diag <- predictor_matrix["regime", colnames(predictor_matrix) != "regime"]
  expect_true(all(off_diag == 1))
  expect_equal(unname(predictor_matrix["regime", "regime"]), 0)
})

test_that("fit_ltmle_one runs on the real 2-wave data", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "fit_ltmle_one.R"))

  pop_data  <- import_data(force = FALSE) |> clean_data() |> preproc_data()
  wide_data <- build_data(pop_data, pre_wave = 1, base_wave = 2, start_wave = 3,
                          final_wave = 4, outcome = "MCS")

  mids <- run_mice(wide_data, m = 1, maxit = 1, seed = 42, outcome = "MCS")

  res <- fit_ltmle_one(mids,
                       imp_idx       = 1,
                       sl_libs       = c("SL.mean", "SL.glm"),
                       outcome       = "MCS",
                       regimes       = toy_regimes,
                       sl_cv_control = toy_cv,
                       gbound        = 0.001)

  # a compact table, not the ltmleMSM fit: 4 regime means + 3 contrasts
  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), 7L)
  expect_named(res, c("imp", "estimand", "type", "estimate", "variance"))
  expect_equal(res$estimand, ltmle_estimands)
  expect_equal(res$type, rep(c("mean", "ate"), c(4, 3)))
  expect_true(all(is.finite(res$estimate)))
  expect_true(all(res$variance > 0))
})

test_that("will a new variable exist after imputation", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))

  pop_data  <- import_data(force = FALSE) |> clean_data() |> preproc_data()
  wide_data <- build_data(pop_data, pre_wave = 1, base_wave = 2, start_wave = 3,
                          final_wave = 4, outcome = "MCS")

  wide_data <- set_exposure(wide_data, exposure = "econ_dist_bin")
  predictor_matrix <- make_counterfactual_matrix(wide_data)

  # mini-mice call to check structure only
  mids <- run_mice(wide_data, m = 2, maxit = 2, seed = 42, outcome = "MCS")

  expect_s3_class(mids, "mids")
  expect_equal(mids$m, 2L)

  as_long_mids <- mice::complete(mids, action = "long", include = TRUE)

  as_long_mids <- as_long_mids |> mutate(new_var = paste0("var_", as.character(sex_dv_fact_base)))

  imp_new <- mice::as.mids(as_long_mids)
  expect_true("new_var" %in% names(mice::complete(imp_new, 1)))

  # no exposure/outcome interaction dropped (a "collinear" age_dv_1 event is harmless)
  expect_true(assert_congenial_terms(mids, a_col = c("econ_dist_bin_0", "econ_dist_bin_1"),
                                     y_col = "sf12mcs_dv_1"))

  expect_true("regime" %in% rownames(predictor_matrix))
  off_diag <- predictor_matrix["regime", colnames(predictor_matrix) != "regime"]
  expect_true(all(off_diag == 1))
  expect_equal(unname(predictor_matrix["regime", "regime"]), 0)
})

test_that("how many strata will be for MAIHDA", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))

  pop_data  <- import_data(force = FALSE) |> clean_data() |> preproc_data()
  wide_data <- build_data(pop_data, pre_wave = 1, base_wave = 2, start_wave = 3,
                          final_wave = 4, outcome = "MCS")

  # mini-mice call to check structure only
  mids <- run_mice(wide_data, m = 5, maxit = 2, seed = 42, outcome = "MCS")

  expect_s3_class(mids, "mids")

  # no exposure/outcome interaction dropped (a "collinear" age_dv_1 event is harmless)
  expect_true(assert_congenial_terms(mids, a_col = c("econ_dist_bin_0", "econ_dist_bin_1"),
                                     y_col = "sf12mcs_dv_1"))

  # strata id and label per imputed dataset
  fits <- purrr::map(seq_len(mids$m), function(i) {
    mice::complete(mids, action = i) |>
      dplyr::mutate(
        # intersectional strata id: unique cell
        strata_id = 100 * as.numeric(sex_dv_fact_base) +
                     10 * as.numeric(race_base) +
                          as.numeric(hiqual_dv_fact_base),
        strata_id = as.factor(strata_id),
        # human-readable stratum name, e.g. "female:white:high"
        strata_label = interaction(sex_dv_fact_base, race_base, hiqual_dv_fact_base,
                                   sep = ":", drop = TRUE)
      )
  })
  # counts per stratum, per imputed dataset (not asserted)
  strata_tables <- purrr::map(fits, function(dat_i) {
    dplyr::count(dat_i, strata_id, strata_label)
  })
})

## ---- DR-learner GATE arm (design 2026-07-30) --------------------------------

### shared confounder set: byte-identical between fit_ltmle_one and estimate_cate
test_that("confounders returns the exact ltmle baseline (W) set", {
  source(here::here("R", "confounders.R"))

  expect_identical(
    confounders("MCS"),
    c("sex_dv_fact_base",
      "hiqual_dv_fact_base",
      "race_base",
      "sf12mcs_dv_base",
      "age_dv_0",
      "gor_dv_fact_0",
      "econ_dist_bin_lagged_0",
      "pcs_lagged_0",
      "dnc_lagged_0",
      "home_owner_lagged_0",
      "econ_benefits_lagged_0",
      "mastat_lagged_0",
      "econ_emp_bin_fact_0",
      "log_income_0")
  )
  # PCS swaps exactly the two outcome-dependent entries
  expect_identical(setdiff(confounders("PCS"), confounders("MCS")),
                   c("sf12pcs_dv_base", "mcs_lagged_0"))
  expect_identical(setdiff(confounders("MCS"), confounders("PCS")),
                   c("sf12mcs_dv_base", "pcs_lagged_0"))
  expect_error(confounders("XXX"))
})

### L1: the confounders measured between A0 and A1
test_that("time_varying_confounders returns the exact ltmle L1 set", {
  source(here::here("R", "confounders.R"))

  expect_identical(
    time_varying_confounders("MCS"),
    c("age_dv_1",
      "gor_dv_fact_1",
      "pcs_lagged_1",
      "dnc_lagged_1",
      "home_owner_lagged_1",
      "econ_benefits_lagged_1",
      "mastat_lagged_1",
      "econ_emp_bin_fact_1",
      "log_income_1")
  )
  # PCS swaps only the other-scale lag; econ_dist_bin_lagged_1 (== A0) never enters
  expect_identical(
    time_varying_confounders("PCS"),
    c("age_dv_1",
      "gor_dv_fact_1",
      "mcs_lagged_1",
      "dnc_lagged_1",
      "home_owner_lagged_1",
      "econ_benefits_lagged_1",
      "mastat_lagged_1",
      "econ_emp_bin_fact_1",
      "log_income_1")
  )
  expect_error(time_varying_confounders("XXX"))
})

### make_strata: pure df -> df, same id arithmetic as the old strata_creation()
test_that("make_strata builds the 12-cell id and label", {
  source(here::here("R", "strata_creation.R"))

  set.seed(1)
  toy <- data.frame(
    sex_dv_base         = factor(sample(c("Female", "Male"), 240, TRUE)),
    race_base           = factor(sample(c("Non.white", "White"), 240, TRUE)),
    hiqual_dv_fact_base = factor(sample(c("High", "Low", "Medium"), 240, TRUE)),
    other_col           = rnorm(240)
  )

  out <- make_strata(toy)

  expect_equal(nrow(out), 240L)
  expect_true(all(c("strata_id", "strata_label") %in% names(out)))
  expect_s3_class(out$strata_id, "factor")
  expect_s3_class(out$strata_label, "factor")
  expect_equal(nlevels(out$strata_id), 12L)
  # id arithmetic: 100*sex + 10*race + hiqual on the factor level indices
  expect_equal(as.character(out$strata_id),
               as.character(100 * as.numeric(toy$sex_dv_base) +
                             10 * as.numeric(toy$race_base) +
                                  as.numeric(toy$hiqual_dv_fact_base)))
  # row order and other columns untouched
  expect_equal(out$other_col, toy$other_col)
  # label format sex:race:educ
  expect_true(all(grepl("^[^:]+:[^:]+:[^:]+$", as.character(out$strata_label))))
})

test_that("make_strata errors when a stratum column is absent", {
  source(here::here("R", "strata_creation.R"))
  expect_error(make_strata(data.frame(sex_dv_base = factor("Female"))), "missing")
})

### DR score primitives (design 3.6-3.7; 7: AIPW identity, sign, Z alignment)
test_that("dr_pseudo_outcome: mean(psi) recovers the true ATE with known nuisances", {
  source(here::here("R", "meta_learner.R"))

  set.seed(101)
  n   <- 20000
  x   <- rnorm(n)
  g   <- plogis(-0.5 + 0.8 * x)
  A   <- rbinom(n, 1, g)
  mu0 <- 50 + 3 * x
  mu1 <- mu0 - 2                       # true conventional ATE = -2
  Y   <- ifelse(A == 1, mu1, mu0) + rnorm(n, 0, 2)

  psi <- dr_pseudo_outcome(A, Y, mu0, mu1, bound_propensity(g))

  # psi is on the Y(0) - Y(1) scale, so the truth here is +2
  expect_equal(mean(psi), 2, tolerance = 0.1)
})

test_that("dr_pseudo_outcome: sign is Y(0) - Y(1), opposite to mu1 - mu0", {
  source(here::here("R", "meta_learner.R"))

  set.seed(102)
  n   <- 500
  g   <- rep(0.4, n)
  A   <- rbinom(n, 1, g)
  mu0 <- rnorm(n, 50, 5)
  mu1 <- mu0 - 3
  Y   <- ifelse(A == 1, mu1, mu0)      # zero residual: IPW terms vanish exactly

  psi <- dr_pseudo_outcome(A, Y, mu0, mu1, g)

  expect_equal(mean(psi), -(mean(mu1) - mean(mu0)), tolerance = 1e-12)
  expect_true(sign(mean(psi)) == -sign(mean(mu1) - mean(mu0)))
})

test_that("dr_pseudo_outcome refuses unbounded propensities", {
  source(here::here("R", "meta_learner.R"))
  expect_error(
    dr_pseudo_outcome(c(0L, 1L), c(1, 2), c(1, 1), c(1, 1), c(0, 0.5)),
    "bound"
  )
})

### carried Task-4 finding: NA g_hat must fail diagnostically, not with R's
### generic "missing value where TRUE/FALSE needed" from the bare comparison
test_that("dr_pseudo_outcome refuses NA propensities with a diagnostic message", {
  source(here::here("R", "meta_learner.R"))
  expect_error(
    dr_pseudo_outcome(c(0L, 1L), c(1, 2), c(1, 1), c(1, 1), c(NA, 0.5)),
    "NA"
  )
})

test_that("bound_propensity clamps to [0.025, 0.975]", {
  source(here::here("R", "meta_learner.R"))
  expect_equal(bound_propensity(c(0, 0.5, 1)), c(0.025, 0.5, 0.975))
})

test_that("SuperLearner $Z is out-of-fold, input-row-ordered and nameless", {
  source(here::here("R", "meta_learner.R"))

  set.seed(103)
  n  <- 300
  X  <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Yc <- 2 * X$x1 - X$x2 + rnorm(n)
  fit <- SuperLearner::SuperLearner(Y = Yc, X = X, family = gaussian(),
                                    SL.library = c("SL.mean", "SL.glm"),
                                    cvControl = list(V = 5))

  # $Z carries no dimnames: index by libraryNames position, never by name
  expect_null(dimnames(fit$Z))
  idx <- which(fit$libraryNames == "SL.glm_All")
  expect_length(idx, 1L)

  # manual V-fold refit reproduces Z exactly, in input row order (design 3.6's
  # empirical check, kept as a regression test against SuperLearner changes)
  manual <- rep(NA_real_, n)
  manual[unlist(fit$validRows)] <- unlist(purrr::map(fit$validRows, \(rows) {
    gfit <- glm(y ~ ., data = data.frame(y = Yc[-rows], X[-rows, ]),
                family = gaussian())
    predict(gfit, newdata = X[rows, ], type = "response")
  }))
  expect_equal(max(abs(manual - fit$Z[, idx])), 0, tolerance = 1e-8)

  # the hard-assertion wrapper: right length through, wrong length stops
  expect_length(sl_oof_predictions(fit, n, "test"), n)
  expect_error(sl_oof_predictions(fit, n - 1, "test"), "rows")

  # NNLS combines linearly, so here the wrapper *is* Z %*% coef
  expect_equal(sl_oof_predictions(fit, n, "test"),
               drop(fit$Z %*% fit$coef), tolerance = 1e-12)
})

test_that("sl_oof_predictions combines on the fit's own method scale", {
  source(here::here("R", "meta_learner.R"))

  set.seed(104)
  n <- 300
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  A <- rbinom(n, 1, plogis(1.2 * X$x1 - 0.8 * X$x2))
  gm <- SuperLearner::SuperLearner(
    Y = A, X = X, family = binomial(),
    SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction"),
    method = "method.NNloglik", cvControl = list(V = 5))

  # method.NNloglik combines on the logit scale, NOT Z %*% coef. Force a split
  # weight vector: the fitted ensemble often collapses onto one learner, where
  # the two formulas coincide and the regression would not bite.
  gm$coef <- stats::setNames(c(0.4, 0.3, 0.3), names(gm$coef))

  # trim must come from $control (default 0.001), not trimLogit()'s own 1e-5
  expect_equal(
    sl_oof_predictions(gm, n, "g"),
    drop(stats::plogis(
      SuperLearner::trimLogit(gm$Z, trim = gm$control$trimLogit) %*% gm$coef)),
    tolerance = 1e-12
  )
  # and that this is materially different from the linear mix
  expect_gt(max(abs(sl_oof_predictions(gm, n, "g") - drop(gm$Z %*% gm$coef))), 0.01)

  # an all-zero weight vector stops instead of silently returning zeros
  gm$coef <- stats::setNames(c(0, 0, 0), names(gm$coef))
  expect_error(sl_oof_predictions(gm, n, "g"), "coefficient")
})

### toy generators for the DR-learner tests (estimate_cate, pool_cate)
# Synthetic wide frame shaped like the real one: all 14 columns of
# confounders("MCS"), the exposure and the outcome. Balanced 12-cell design,
# every A x S cell populated (deterministic 40% exposed per cell).
# tau_by_cell: length-12, the treatment effect (conventional Y(1) - Y(0)
# scale) per stratum, indexed in strata_id order (ids sorted ascending).
make_toy_wide_df <- function(n_per_cell = 100, tau_by_cell = rep(-2, 12), seed = 1) {
  set.seed(seed)
  grid <- expand.grid(
    sex_dv_base         = factor(c("Female", "Male")),
    race_base           = factor(c("Non.white", "White")),
    hiqual_dv_fact_base = factor(c("High", "Low", "Medium")),
    KEEP.OUT.ATTRS = FALSE
  )
  df <- grid[rep(seq_len(12), each = n_per_cell), ]
  n  <- nrow(df)

  df$sf12mcs_dv_base        <- rnorm(n, 50, 8)
  df$age_dv_0               <- rnorm(n, 45, 10)
  df$gor_dv_fact_0          <- factor(sample(c("London", "North.East", "Wales"), n, TRUE))
  df$econ_dist_bin_lagged_0 <- factor(rbinom(n, 1, 0.2))
  df$pcs_lagged_0           <- rnorm(n, 50, 8)
  df$dnc_lagged_0           <- factor(sample(c("Zero", "One", "Two.plus"), n, TRUE))
  df$home_owner_lagged_0    <- factor(rbinom(n, 1, 0.6))
  df$econ_benefits_lagged_0 <- factor(rbinom(n, 1, 0.3))
  df$mastat_lagged_0        <- factor(sample(c("Partnered", "Not.partnered"), n, TRUE))
  df$econ_emp_bin_fact_0    <- factor(sample(c("Employed", "Not.employed"), n, TRUE))
  df$log_income_0           <- rnorm(n, 10, 1)

  # exposure: exactly 40% treated inside every cell (all 24 A x S cells filled)
  n_tr <- round(0.4 * n_per_cell)
  a <- purrr::list_c(purrr::map(seq_len(12), \(j) {
    sample(rep(c(1L, 0L), c(n_tr, n_per_cell - n_tr)))
  }))
  df$econ_dist_bin_0 <- factor(a)

  # map each row's stratum to its rank in strata_id (= sorted 100s+10r+h) order
  sid      <- 100 * as.numeric(df$sex_dv_base) +
               10 * as.numeric(df$race_base) +
                    as.numeric(df$hiqual_dv_fact_base)
  sid_rank <- match(sid, sort(unique(sid)))

  df$sf12mcs_dv_0 <- 0.5 * df$sf12mcs_dv_base +
    2 * (df$sex_dv_base == "Female") +
    tau_by_cell[sid_rank] * a +
    rnorm(n, 0, 1)

  rownames(df) <- NULL
  df
}

# Punch a few holes and impute so estimate_cate sees a genuine mids object.
make_toy_wide_mids <- function(df, m = 2, seed = 7) {
  set.seed(seed)
  df_na <- df
  df_na$sf12mcs_dv_0[sample(nrow(df), round(0.05 * nrow(df)))] <- NA
  df_na$log_income_0[sample(nrow(df), round(0.05 * nrow(df)))] <- NA
  mice::mice(df_na, m = m, maxit = 1, seed = seed, printFlag = FALSE)
}

### toy generators for the two-time-point LTMLE tests (fit_ltmle_one,
### run_mice formulas, make_counterfactual_matrix, toy pipeline)
# Synthetic frame shaped like build_data()'s MCS output: the 27 wide_cols_mcs
# columns with their types, plus the make_wide()/set_exposure() attributes.
# 12 balanced S cells. A0 ~ W; L1 ~ A0 + W (only income and employment respond
# to A0); A1 ~ A0 + L1; Y ~ W + L1 + A0 + A1. P(A0) ~ 0.28 and P(A1) ~ 0.23,
# so all 48 A0 x A1 x S cells are populated (min 33 at 400 per cell).
# attr(, "truth"): sample-average regime means E[Y(a0, a1)], with the
# counterfactual L1 drawn from the same noise.
make_toy_2wave_df <- function(n_per_cell = 400, seed = 1) {
  set.seed(seed)
  d <- data.table::as.data.table(expand.grid(
    sex_dv_fact_base    = factor(c("Female", "Male")),
    hiqual_dv_fact_base = factor(c("High", "Medium", "Low"), levels = c("High", "Medium", "Low")),
    race_base           = factor(c("Non.white", "White")),
    KEEP.OUT.ATTRS = FALSE
  ))[rep(seq_len(12), each = n_per_cell)]
  n    <- nrow(d)
  lab  <- function(x, lv) factor(lv[x + 1L], levels = lv)   # 0/1 -> labelled factor
  flip <- function(x, p) { i <- runif(n) < p; x[i] <- 1L - x[i]; x }
  male <- d$sex_dv_fact_base == "Male"
  low  <- d$hiqual_dv_fact_base == "Low"

  base   <- 50 + 2 * male + rnorm(n, 0, 8)
  age0   <- as.integer(round(runif(n, 20, 70)))
  gor0   <- sample(c("London", "North.East", "Wales"), n, TRUE)
  pcs0   <- 50 - 0.1 * (age0 - 45) + rnorm(n, 0, 8)
  dist_l <- rbinom(n, 1, 0.25)
  dnc0   <- sample(c("One", "Two.or.more", "Zero"), n, TRUE)
  own0   <- rbinom(n, 1, 0.6)
  ben0   <- rbinom(n, 1, 0.3)
  mar0   <- rbinom(n, 1, 0.6)
  emp0   <- rbinom(n, 1, 0.7)
  inc0   <- 10 - 0.3 * ben0 + 0.3 * emp0 + rnorm(n, 0, 0.6)

  a0 <- rbinom(n, 1, plogis(-1.4 + 1.2 * dist_l - 0.6 * (inc0 - 10) -
                              0.02 * (pcs0 - 50) + 0.3 * ben0 - 0.3 * own0 +
                              0.2 * (d$race_base == "Non.white") + 0.2 * low))

  # L1 as functions of a0 so the truth can reuse the same noise
  e_inc <- rnorm(n, 0, 0.3)
  u_emp <- runif(n)
  inc1  <- function(a) 3 + 0.7 * inc0 - 0.25 * a + e_inc
  emp1  <- function(a) as.integer(u_emp < plogis(-1 + 2.5 * emp0 - 0.6 * a))
  # the gap varies so mice does not log age_dv_1 as collinear (cor < .999)
  age1  <- age0 + sample(0:3, n, TRUE)
  gor1  <- ifelse(runif(n) < 0.03, sample(c("London", "North.East", "Wales"), n, TRUE), gor0)
  pcs1  <- 20 + 0.6 * pcs0 + rnorm(n, 0, 5)
  dnc1  <- ifelse(runif(n) < 0.1, sample(c("One", "Two.or.more", "Zero"), n, TRUE), dnc0)
  ben1  <- flip(ben0, 0.1)

  a1 <- rbinom(n, 1, plogis(-1.5 + 1.1 * a0 - 0.6 * (inc1(a0) - 10) -
                              0.02 * (pcs1 - 50) - 0.4 * emp1(a0) + 0.3 * ben1))

  # A0:A1 kept small: the glm Q-models and gFormulaMI's pmm are main-terms
  # only. A0 and A1 effects far apart, so "before" and "now" means differ by ~4
  # and a swapped c(A0, A1) order misses the truth.
  mu_y <- function(a0, a1) {
    50 + 0.3 * (base - 50) + 0.1 * (pcs1 - 50) + 1.5 * (inc1(a0) - 10) +
      1 * emp1(a0) + 1 * male - 0.5 * low - 1.5 * a0 - 6 * a1 - 1 * a0 * a1
  }

  d[, `:=`(
    sf12mcs_dv_base        = base,
    age_dv_0               = age0,
    gor_dv_fact_0          = factor(gor0),
    pcs_lagged_0           = pcs0,
    econ_dist_bin_0        = a0,
    econ_dist_bin_lagged_0 = dist_l,
    dnc_lagged_0           = factor(dnc0, levels = c("One", "Two.or.more", "Zero")),
    home_owner_lagged_0    = lab(own0, c("Renter", "Owner")),
    econ_benefits_lagged_0 = lab(1L - ben0, c("Benefits", "No.benefits")),
    mastat_lagged_0        = lab(mar0, c("Not.partnered", "Partnered")),
    econ_emp_bin_fact_0    = factor(emp0, levels = 0:1),
    log_income_0           = inc0,
    # start-wave outcome: an auxiliary only, it never feeds A1 or Y
    sf12mcs_dv_0           = 48 + 0.4 * (base - 50) - 2 * a0 + rnorm(n, 0, 6),
    age_dv_1               = age1,
    gor_dv_fact_1          = factor(gor1),
    pcs_lagged_1           = pcs1,
    econ_dist_bin_1        = a1,
    dnc_lagged_1           = factor(dnc1, levels = c("One", "Two.or.more", "Zero")),
    home_owner_lagged_1    = lab(flip(own0, 0.05), c("Renter", "Owner")),
    econ_benefits_lagged_1 = lab(1L - ben1, c("Benefits", "No.benefits")),
    mastat_lagged_1        = lab(flip(mar0, 0.05), c("Not.partnered", "Partnered")),
    econ_emp_bin_fact_1    = factor(emp1(a0), levels = 0:1),
    log_income_1           = inc1(a0),
    sf12mcs_dv_1           = mu_y(a0, a1) + rnorm(n, 0, 5)
  )]
  data.table::setcolorder(d, wide_cols_mcs)

  lagged <- grep("lagged", wide_cols_mcs, value = TRUE)
  purrr::iwalk(list(
    n_tvars           = 11L,
    n_base            = 4L,
    baseline_vars     = wide_cols_mcs[1:4],
    time_lagged       = lagged,
    outcome_var       = "sf12mcs_dv_1",
    outcome_baseline  = "sf12mcs_dv_base",
    time_points       = 2L,
    intermediate_vars = setdiff(wide_cols_mcs[-(1:4)], lagged),
    exposure_vars     = c("econ_dist_bin_0", "econ_dist_bin_1"),
    truth             = purrr::map_dbl(toy_regimes, \(r) mean(mu_y(r[1], r[2])))
  ), \(v, k) data.table::setattr(d, k, v))
  d
}

# Punch ~4% holes into both exposures, the outcome and a confounder, then
# impute with the real run_mice() so the mids carries its congenial formulas.
make_toy_2wave_mids <- function(df, m = 2, seed = 7) {
  set.seed(seed)
  df_na <- data.table::copy(df)
  purrr::walk(c("econ_dist_bin_0", "econ_dist_bin_1", "sf12mcs_dv_1", "log_income_0"),
              \(v) data.table::set(df_na, i = sample(nrow(df_na), round(0.04 * nrow(df_na))),
                                   j = v, value = NA))
  run_mice(df_na, m = m, maxit = 1, seed = seed, outcome = "MCS")
}

### design 7 "factor handling": W contains factors.
### Distinct per-cell tau (design review #1a): a constant tau_by_cell (e.g. the
### generator's default rep(-2, 12)) cannot tell a correct sign from a flipped
### one, since mean(psi) == n_j-weighted mean of the cell GATEs is an algebraic
### identity for *any* psi, negated included. Mutation-tested by the
### design-2026-07-30 reviewer: a negated-psi mutant and a mu0/mu1 arm-swap
### mutant both passed every original assertion. Distinct, signed per-cell
### effects plus an explicit sign/recovery check close that gap.
test_that("estimate_cate runs on a small mids with factor confounders", {
  source(here::here("R", "confounders.R"))
  source(here::here("R", "strata_creation.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "meta_learner.R"))

  tau_by_cell <- c(-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6)
  df   <- make_toy_wide_df(n_per_cell = 100, tau_by_cell = tau_by_cell, seed = 11)
  mids <- make_toy_wide_mids(df, m = 2, seed = 11)

  res <- estimate_cate(mids, imp_idx = 1, sl_libs = c("SL.mean", "SL.glm"),
                       outcome = "MCS")

  expect_named(res, c("gate", "blp", "ate", "wald"))
  expect_equal(nrow(res$gate), 12L)
  expect_setequal(
    names(res$gate),
    c("imp", "strata_id", "strata_label", "n_j", "n_exposed", "pct_exposed",
      "min_g", "median_g", "share_g_at_bound", "estimate", "se")
  )
  expect_false(anyNA(res$gate$estimate))
  expect_true(all(res$gate$se > 0))
  expect_equal(sum(res$gate$n_j), nrow(df))
  expect_true(all(res$gate$min_g >= CATE_EPS & res$gate$min_g <= 1 - CATE_EPS))

  expect_equal(nrow(res$blp), 2L)
  expect_equal(res$blp$term, c("beta1_ate", "beta2_hte"))

  expect_equal(res$wald$df1, 11)

  # exact identity: mean(psi) is the n_j-weighted mean of the cell GATEs
  expect_equal(res$ate$estimate,
               sum(res$gate$estimate * res$gate$n_j) / sum(res$gate$n_j),
               tolerance = 1e-10)

  # sign and per-cell recovery: gate$estimate is Y(0)-Y(1), i.e. -tau_by_cell.
  # A sign-flip mutant misses by ~2x|tau| (>= 2 everywhere here). Measured max
  # error under the true code: 0.88 (no-imputation-noise floor ~0.36, through
  # make_toy_wide_mids() ~0.70), so tolerance = 1.0 is margin, not vacuous.
  expect_equal(res$gate$estimate, -tau_by_cell, tolerance = 1.0)
  # beta2_hte pins near +1 when the tau proxy tracks the true heterogeneity on
  # the conventional scale; it flips to about -1 under a sign error.
  expect_gt(res$blp$estimate[res$blp$term == "beta2_hte"], 0.5)
})

### design review #1b: make_toy_wide_df() randomizes A at a fixed 40% per
### cell, so ghat is exact and AIPW is invariant to which fitted model plays
### mu0 vs mu1: double robustness under an exactly correct g, not coincidence
### (confirmed algebraically and empirically: an arm-swap mutant passed the
### test above). Catching an arm swap needs (i) A confounded by a continuous W
### the tiny SL.glm/SL.mean library cannot capture (ghat NOT exact) and (ii) a
### genuine W x A outcome interaction (mu0/mu1 differ as functions of W, not
### by a location shift); otherwise the swap has nothing to bite on. A
### SEPARATE generator: make_toy_wide_df() stays unchanged for Tasks 6/12.
# A ~ quadratic logit in age (b1*age_c + b2*age_c^2): linear-only SL.glm
# cannot recover it, so ghat stays biased at any n. Y's treatment effect is
# linear in age (tau0 + tau_age*(age-45)), a W x A interaction SL.glm CAN fit
# for mu0/mu1, so the arm swap is the only error source. Calibrated
# empirically: swapping mu0_hat/mu1_hat inflates mean(psi) by ~+0.9 to +1.1
# over the truth; correctly assigned arms recover it within ~0.15.
make_toy_confounded_df <- function(n_per_cell = 60, tau0 = -3, tau_age = -0.4,
                                    b1 = 2.0, b2 = -3.0, seed = 31) {
  set.seed(seed)
  grid <- expand.grid(
    sex_dv_base         = factor(c("Female", "Male")),
    race_base           = factor(c("Non.white", "White")),
    hiqual_dv_fact_base = factor(c("High", "Low", "Medium")),
    KEEP.OUT.ATTRS = FALSE
  )
  df <- grid[rep(seq_len(12), each = n_per_cell), ]
  n  <- nrow(df)

  df$sf12mcs_dv_base        <- rnorm(n, 50, 8)
  df$age_dv_0               <- rnorm(n, 45, 10)
  df$gor_dv_fact_0          <- factor(sample(c("London", "North.East", "Wales"), n, TRUE))
  df$econ_dist_bin_lagged_0 <- factor(rbinom(n, 1, 0.2))
  df$pcs_lagged_0           <- rnorm(n, 50, 8)
  df$dnc_lagged_0           <- factor(sample(c("Zero", "One", "Two.plus"), n, TRUE))
  df$home_owner_lagged_0    <- factor(rbinom(n, 1, 0.6))
  df$econ_benefits_lagged_0 <- factor(rbinom(n, 1, 0.3))
  df$mastat_lagged_0        <- factor(sample(c("Partnered", "Not.partnered"), n, TRUE))
  df$econ_emp_bin_fact_0    <- factor(sample(c("Employed", "Not.employed"), n, TRUE))
  df$log_income_0           <- rnorm(n, 10, 1)

  age_c  <- (df$age_dv_0 - 45) / 10
  g_true <- plogis(-0.4 + b1 * age_c + b2 * age_c^2)
  a      <- rbinom(n, 1, g_true)
  df$econ_dist_bin_0 <- factor(a)

  df$sf12mcs_dv_0 <- 0.5 * df$sf12mcs_dv_base +
    2 * (df$sex_dv_base == "Female") +
    0.1 * (df$age_dv_0 - 45) +
    a * (tau0 + tau_age * (df$age_dv_0 - 45)) +
    rnorm(n, 0, 1)

  rownames(df) <- NULL
  df
}

test_that("estimate_cate recovers the ATE under a confounded, W x A heterogeneous DGP", {
  source(here::here("R", "confounders.R"))
  source(here::here("R", "strata_creation.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "meta_learner.R"))

  tau0    <- -3
  tau_age <- -0.4
  df   <- make_toy_confounded_df(n_per_cell = 60, tau0 = tau0, tau_age = tau_age, seed = 31)
  mids <- make_toy_wide_mids(df, m = 2, seed = 31)

  res <- estimate_cate(mids, imp_idx = 1, sl_libs = c("SL.mean", "SL.glm"),
                       outcome = "MCS")

  # true ATE on the Y(0)-Y(1) scale over this sample's realized ages
  # (tau0/tau_age fully determine it). No per-cell truth: age is iid of
  # strata, so all 12 cells share it up to sampling noise.
  truth <- -mean(tau0 + tau_age * (df$age_dv_0 - 45))
  # expect_equal(tolerance =) on a *scalar* is RELATIVE (all.equal.numeric),
  # not absolute: confirmed empirically it would silently pass the mutant this
  # test exists to catch. So assert the absolute gap: correct arms land within
  # ~0.05-0.15 at this n/seed; a mu0/mu1 arm swap adds ~+1.0 (measured). 0.5
  # sits well inside that gap.
  expect_lt(abs(res$ate$estimate - truth), 0.5)
})

### design 7 "recovery": inject known cell effects; pooled GATE_j recovers them
test_that("pool_cate recovers injected cell-specific effects", {
  source(here::here("R", "confounders.R"))
  source(here::here("R", "strata_creation.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "meta_learner.R"))
  source(here::here("R", "pool_cate.R"))

  tau  <- seq(-5, 3, length.out = 12)   # conventional Y(1)-Y(0), one per stratum
  df   <- make_toy_wide_df(n_per_cell = 150, tau_by_cell = tau, seed = 21)
  mids <- make_toy_wide_mids(df, m = 2, seed = 21)

  fits   <- purrr::map(1:2, \(i) {
    estimate_cate(mids, imp_idx = i, sl_libs = c("SL.mean", "SL.glm"),
                  outcome = "MCS")
  })
  pooled <- pool_cate(fits)

  expect_named(pooled, c("gate", "blp", "ate", "wald", "m"))
  expect_equal(pooled$m, 2L)
  expect_equal(nrow(pooled$gate), 12L)
  expect_equal(nrow(pooled$blp), 2L)
  expect_equal(nrow(pooled$wald), 2L)   # one row per imputation, unpooled

  # GATE_j is on the Y(0) - Y(1) scale: recover -tau, in strata_id order
  expect_lt(max(abs(pooled$gate$estimate - (-tau))), 1)
  expect_gt(cor(pooled$gate$estimate, -tau), 0.9)

  expect_true(all(pooled$gate$ll < pooled$gate$estimate))
  expect_true(all(pooled$gate$ul > pooled$gate$estimate))
  expect_true(all(pooled$gate$n_j_min <= pooled$gate$n_j_max))

  # marginal mean(psi) pooled too, and near the n-weighted mean of -tau
  expect_equal(nrow(pooled$ate), 1L)
  expect_lt(abs(pooled$ate$estimate - mean(-tau)), 0.5)
})

### sanitize_factor_levels: fix formula-hostile levels, leave numeric codes alone
test_that("sanitize_factor_levels fixes only levels with formula-hostile characters", {
  source(here::here("R", "helpers.R"))

  toy <- data.frame(
    race_base            = factor(c("White", "Non-white", "White")),
    dnc_lagged           = factor(c("Zero", "One", "2+")),
    econ_benefits_lagged = factor(c("No benefits", "Benefits", "Benefits")),
    hiqual_dv_fact_base  = factor(c("High", "Medium", "Low")),
    econ_dist_bin_0      = factor(c("0", "1", "0")),
    age_dv_0             = c(30, 40, 50),
    stringsAsFactors     = FALSE
  )

  out <- sanitize_factor_levels(toy)

  # hostile levels rewritten
  expect_equal(levels(out$race_base),            c("Non.white", "White"))
  expect_equal(levels(out$dnc_lagged),           c("One", "X2.", "Zero"))
  expect_equal(levels(out$econ_benefits_lagged), c("Benefits", "No.benefits"))

  # already-clean factor untouched
  expect_equal(levels(out$hiqual_dv_fact_base), levels(toy$hiqual_dv_fact_base))

  # THE LANDMINE: numeric-coded exposure must survive as.integer(as.character(.))
  expect_equal(levels(out$econ_dist_bin_0), c("0", "1"))
  expect_equal(as.integer(as.character(out$econ_dist_bin_0)), c(0L, 1L, 0L))

  # structure preserved
  expect_equal(nrow(out), 3L)
  expect_equal(names(out), names(toy))
  expect_equal(out$age_dv_0, toy$age_dv_0)
  # values, not just levels, still line up
  expect_equal(as.character(out$race_base), c("White", "Non.white", "White"))
})

test_that("sanitize_factor_levels honours skip and handles data.table without warning", {
  source(here::here("R", "helpers.R"))

  dt <- data.table::data.table(
    race_base = factor(c("White", "Non-white")),
    keep_me   = factor(c("a b", "c d"))
  )
  out <- expect_no_warning(sanitize_factor_levels(dt, skip = "keep_me"))

  expect_equal(levels(out$race_base), c("Non.white", "White"))
  expect_equal(levels(out$keep_me), c("a b", "c d"))   # skipped
  expect_true(data.table::is.data.table(out))
  # input not modified by reference
  expect_equal(levels(dt$race_base), c("Non-white", "White"))
})

### review finding: factor(new_labels[...]) without `levels =` takes its levels
### from observed rows only, so a defined-but-unused level silently vanishes
### instead of being relabelled and kept, breaking the "levels replaced"
### contract. Guards the `levels = sort(unique(new_labels))` fix.
test_that("sanitize_factor_levels keeps a defined-but-unused level, relabelled", {
  source(here::here("R", "helpers.R"))

  # "2+" is a defined level with zero rows in this particular frame
  toy <- data.frame(
    dnc_lagged = factor(c("Zero", "One", "Zero"), levels = c("Zero", "One", "2+"))
  )
  expect_equal(nlevels(toy$dnc_lagged), 3L)

  out <- sanitize_factor_levels(toy)

  # level count preserved: the unused level is relabelled, not dropped
  expect_equal(nlevels(out$dnc_lagged), 3L)
  expect_equal(levels(out$dnc_lagged), c("One", "X2.", "Zero"))
  # observed values still line up with the original
  expect_equal(as.character(out$dnc_lagged), c("Zero", "One", "Zero"))
})

### same root cause, sharper failure mode: on a zero-row frame every hostile
### column would otherwise come out with zero levels (factor(character(0)))
test_that("sanitize_factor_levels keeps full level sets on a zero-row frame", {
  source(here::here("R", "helpers.R"))

  toy <- data.frame(
    race_base  = factor(character(0), levels = c("White", "Non-white")),
    dnc_lagged = factor(character(0), levels = c("Zero", "One", "2+"))
  )
  expect_equal(nrow(toy), 0L)

  out <- sanitize_factor_levels(toy)

  expect_equal(nrow(out), 0L)
  expect_equal(levels(out$race_base), c("Non.white", "White"))
  expect_equal(levels(out$dnc_lagged), c("One", "X2.", "Zero"))
})

### ---- assert_congenial_terms ---------------------------------------------
### mids$loggedEvents$out is a *display* string, not a data structure:
### mice:::remove.lindep() collapses every predictor it drops for one
### (it, im, dep) into a single comma-joined `out`. A programmatic check must
### re-parse it per term and cope with the two exits where remove.lindep names
### no term at all.

# minimal one-row loggedEvents in mice's own column layout
logged_event <- function(out, dep = "sf12mcs_dv_1", meth = "collinear") {
  list(loggedEvents = data.frame(it = 1L, im = 1L, dep = dep, meth = meth,
                                 out = out, stringsAsFactors = FALSE))
}

test_that("assert_congenial_terms requires both node names", {
  source(here::here("R", "run_mice.R"))

  # No defaults: a_col/y_col are always resolved by the caller, so a default
  # exists only to be silently wrong on a window whose suffixes differ.
  expect_error(
    assert_congenial_terms(logged_event("econ_dist_bin_0:race_baseWhite"),
                           a_col = c("econ_dist_bin_0", "econ_dist_bin_1")),
    "y_col"
  )
  expect_error(
    assert_congenial_terms(logged_event("econ_dist_bin_0:race_baseWhite"),
                           y_col = "sf12mcs_dv_1"),
    "a_col"
  )
})

test_that("assert_congenial_terms stops when the Y model drops an A:S interaction", {
  source(here::here("R", "run_mice.R"))

  A <- c("econ_dist_bin_0", "econ_dist_bin_1"); Y <- "sf12mcs_dv_1"

  err <- expect_error(
    assert_congenial_terms(logged_event("econ_dist_bin_0:race_baseWhite"), A, Y),
    "congenial"
  )
  # the dropped term must appear in the message -- deciding rung 1 vs rung 2
  # needs to know *which* term went
  expect_true(grepl("econ_dist_bin_0:race_baseWhite",
                    conditionMessage(err), fixed = TRUE))
  expect_true(grepl("1 exposure-interaction term(s)",
                    conditionMessage(err), fixed = TRUE))
  # A1 is guarded as well as A0
  expect_error(
    assert_congenial_terms(logged_event("econ_dist_bin_1:sex_dv_fact_baseMale"), A, Y),
    "congenial"
  )
})

test_that("assert_congenial_terms stops when the Y model drops an A0:A1:S interaction", {
  source(here::here("R", "run_mice.R"))

  # the saturated A0 x A1 x S cells: the term names both exposures, yet it is
  # one dropped term
  err <- expect_error(
    assert_congenial_terms(
      logged_event("econ_dist_bin_0:econ_dist_bin_1:race_baseWhite"),
      a_col = c("econ_dist_bin_0", "econ_dist_bin_1"), y_col = "sf12mcs_dv_1"),
    "congenial"
  )
  expect_true(grepl("1 exposure-interaction term(s)",
                    conditionMessage(err), fixed = TRUE))
})

test_that("assert_congenial_terms stops when the A model drops a Y:S interaction", {
  source(here::here("R", "run_mice.R"))

  # run_mice puts y * A_other * S in each A model, so a break there is logged
  # under the OUTCOME name with dep = the exposure. Matching only on a_col
  # misses it. y enters centred, so the label carries an I() wrapper.
  expect_error(
    assert_congenial_terms(logged_event("I(sf12mcs_dv_1 - 47.44):sex_dv_fact_baseMale",
                                        dep = "econ_dist_bin_0"),
                           a_col = c("econ_dist_bin_0", "econ_dist_bin_1"),
                           y_col = "sf12mcs_dv_1"),
    "congenial"
  )
})

test_that("assert_congenial_terms counts dropped terms, not logged rows", {
  source(here::here("R", "run_mice.R"))

  err <- expect_error(
    assert_congenial_terms(
      logged_event(paste("econ_dist_bin_0:sex_dv_fact_baseMale",
                         "econ_dist_bin_0:race_baseWhite",
                         "econ_dist_bin_0:hiqual_dv_fact_baseLow",
                         sep = ", ")),
      a_col = c("econ_dist_bin_0", "econ_dist_bin_1"), y_col = "sf12mcs_dv_1"),
    "congenial"
  )
  expect_true(grepl("3 exposure-interaction term(s)",
                    conditionMessage(err), fixed = TRUE))
})

test_that("assert_congenial_terms catches remove.lindep's unnamed wipeout", {
  source(here::here("R", "run_mice.R"))

  # remove.lindep() gives up without naming any term when nothing survives.
  # Every interaction went with it, but the message carries no ":" and no node
  # name, so a term-name scan alone reports the run as clean.
  expect_error(
    assert_congenial_terms(
      logged_event("All predictors are constant or have too high correlation."),
      a_col = c("econ_dist_bin_0", "econ_dist_bin_1"), y_col = "sf12mcs_dv_1"),
    "wipeout"
  )
})

test_that("assert_congenial_terms ignores a main-effect drop sharing a row with an unrelated interaction", {
  source(here::here("R", "run_mice.R"))

  # The A main effect was dropped from some other variable's model, and an
  # unrelated S:S term went with it. Matching the whole row for "a_col AND :"
  # calls this an A-interaction drop and halts a run that is fine.
  expect_true(
    assert_congenial_terms(
      logged_event("econ_dist_bin_0, sex_dv_fact_baseMale:race_baseWhite",
                   dep = "log_income_0"),
      a_col = c("econ_dist_bin_0", "econ_dist_bin_1"), y_col = "sf12mcs_dv_1")
  )
})

test_that("assert_congenial_terms matches on the node names it is given, not on wave 0", {
  source(here::here("R", "run_mice.R"))

  # Other suffixes, e.g. _3/_4. The guard must follow the names it is handed.
  wave4 <- logged_event("econ_dist_bin_4:race_baseWhite", dep = "sf12mcs_dv_4")
  expect_error(
    assert_congenial_terms(wave4, a_col = c("econ_dist_bin_3", "econ_dist_bin_4"),
                           y_col = "sf12mcs_dv_4"),
    "congenial"
  )
  # PCS resolves to a different outcome stem, same wave
  expect_error(
    assert_congenial_terms(logged_event("sf12pcs_dv_4:sex_dv_fact_baseMale",
                                        dep = "econ_dist_bin_4"),
                           a_col = c("econ_dist_bin_3", "econ_dist_bin_4"),
                           y_col = "sf12pcs_dv_4"),
    "congenial"
  )
})

test_that("assert_congenial_terms passes on clean or unrelated logged events", {
  source(here::here("R", "run_mice.R"))

  A <- c("econ_dist_bin_0", "econ_dist_bin_1"); Y <- "sf12mcs_dv_1"

  expect_true(assert_congenial_terms(list(loggedEvents = NULL), A, Y))
  expect_true(assert_congenial_terms(
    logged_event("gor_dv_fact_0", dep = "log_income_0", meth = "pmm"), A, Y))
  # the real-data setup event: age_dv_1 (~ age_dv_0 + 1) dropped as collinear
  expect_true(assert_congenial_terms(logged_event("age_dv_1", dep = ""), A, Y))
  # a main-effect drop of the exposure itself is not an interaction drop
  expect_true(assert_congenial_terms(
    logged_event("econ_dist_bin_lagged_0", meth = "pmm"), A, Y))
  # ...and neither is a lagged-exposure interaction: econ_dist_bin_lagged_0 is a
  # different node from the A nodes, even though it shares the stem
  expect_true(assert_congenial_terms(
    logged_event("econ_dist_bin_lagged_0:gor_dv_fact_0North.East", meth = "pmm"), A, Y))
})

### run_mice: formulas carry the saturated A0 x A1 x S interaction
test_that("run_mice builds Y and A formulas with the saturated strata interaction", {
  source(here::here("R", "run_mice.R"))

  # m = 1, maxit = 1: model structure, not convergence
  mids <- make_toy_2wave_mids(make_toy_2wave_df(n_per_cell = 100), m = 1)

  expect_s3_class(mids, "mids")

  rhs_terms <- function(v) attr(terms(mids$formulas[[v]]), "term.labels")

  # Term-label variable ORDER in terms() follows each variable's main-effect
  # position in the whole formula, not the order written in the interaction,
  # and the A models' y enters centred, as "I(sf12mcs_dv_1 - 47.44)". So match
  # on the *set* of variables in a term rather than an exact string.
  has_term <- function(term_labels, vars) {
    any(purrr::map_lgl(term_labels, \(tl) setequal(all.vars(str2lang(tl)), vars)))
  }
  S <- c("sex_dv_fact_base", "race_base", "hiqual_dv_fact_base")

  # Y model: saturated A0 x A1 x S
  expect_true(has_term(rhs_terms("sf12mcs_dv_1"), c("econ_dist_bin_0", "econ_dist_bin_1", S)))
  # A models: y x A_other x S, and never the exposure being imputed
  expect_true(has_term(rhs_terms("econ_dist_bin_0"), c("sf12mcs_dv_1", "econ_dist_bin_1", S)))
  expect_true(has_term(rhs_terms("econ_dist_bin_1"), c("sf12mcs_dv_1", "econ_dist_bin_0", S)))
  expect_false("econ_dist_bin_0" %in% all.vars(mids$formulas[["econ_dist_bin_0"]][[3]]))
  expect_false("econ_dist_bin_1" %in% all.vars(mids$formulas[["econ_dist_bin_1"]][[3]]))
})

### run_mice: methods_by_group is keyed without the _base suffix, so it must still reach S
test_that("run_mice imputes sex, hiqual and race with their methods_by_group methods", {
  source(here::here("R", "run_mice.R"))

  S  <- c("sex_dv_fact_base", "hiqual_dv_fact_base", "race_base")
  df <- make_toy_2wave_df(n_per_cell = 100)
  # mice blanks the method of a complete column, so give each S column some NAs
  set.seed(3)
  purrr::walk(S, \(v) data.table::set(df, i = sample(nrow(df), 20L), j = v, value = NA))
  mids <- make_toy_2wave_mids(df, m = 1)

  expect_identical(unname(mids$method[S]), c("logreg", "polr", "logreg"))
  expect_false(anyNA(mice::complete(mids, 1)[S]))
})

## ---- two-time-point LTMLE (toy data) -----------------------------------------

### one saturated ltmleMSM call must reproduce one ltmle() call per regime
test_that("fit_ltmle_one equals separate ltmle fits per regime", {
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "fit_ltmle_one.R"))

  mids <- make_toy_2wave_mids(make_toy_2wave_df(), m = 1)
  # "glm" bypasses SuperLearner, so both sides are deterministic
  r <- fit_ltmle_one(mids, 1, "glm", "MCS", toy_regimes, toy_cv, 0.001)

  expect_s3_class(r, "data.table")
  expect_named(r, c("imp", "estimand", "type", "estimate", "variance"))
  expect_equal(r$estimand, ltmle_estimands)
  expect_equal(r$type, rep(c("mean", "ate"), c(4, 3)))
  expect_equal(r$imp, rep(1L, 7))

  # the ltmle data the contract specifies: W, A0, L1, A1, Y in node order, A integer
  A <- c("econ_dist_bin_0", "econ_dist_bin_1")
  L <- time_varying_confounders("MCS")
  d <- as.data.frame(mice::complete(mids, 1))[, c(confounders("MCS"), A[1], L, A[2], "sf12mcs_dv_1")]
  d[A] <- purrr::map(d[A], \(a) as.integer(as.character(a)))
  ltmle_sep <- function(abar) {
    suppressMessages(ltmle::ltmle(d, Anodes = A, Lnodes = L, Ynodes = "sf12mcs_dv_1",
                                  abar = abar, SL.library = "glm", gbounds = c(0.001, 1),
                                  variance.method = "ic", estimate.time = FALSE))
  }

  purrr::iwalk(toy_regimes, \(abar, k) {
    s <- summary(ltmle_sep(abar))$treatment
    expect_equal(r$estimate[r$estimand == k], unname(s$estimate), tolerance = 1e-6)
    expect_equal(sqrt(r$variance[r$estimand == k]), s$std.dev, tolerance = 1e-6)
  })
  # contrasts: ltmle(abar = list(treatment, control)), control = the reference
  purrr::iwalk(toy_regimes[-1], \(abar, k) {
    s   <- summary(ltmle_sep(list(abar, toy_regimes$never)))$effect.measures$ATE
    lbl <- paste(k, "-", "never")
    expect_equal(r$estimate[r$estimand == lbl], s$estimate, tolerance = 1e-6)
    expect_equal(sqrt(r$variance[r$estimand == lbl]), s$std.dev, tolerance = 1e-6)
  })
})

test_that("fit_ltmle_one refuses stacked imputations and an unnamed cv control", {
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "fit_ltmle_one.R"))

  mids <- make_toy_2wave_mids(make_toy_2wave_df(n_per_cell = 100), m = 2)

  # complete(action = 1:2) would silently stack the imputations
  expect_error(fit_ltmle_one(mids, 1:2, "glm", "MCS", toy_regimes, toy_cv, 0.001),
               "imp_idx")
  # ltmle checks only the names it is given, so SuperLearner would read an
  # unnamed list by position
  expect_error(fit_ltmle_one(mids, 1, "glm", "MCS", toy_regimes,
                             list(5L, FALSE, TRUE, NULL), 0.001),
               "sl_cv_control")
})

test_that("pool_ltmle applies Rubin's rules and refuses a bare table", {
  source(here::here("R", "pool_ltmle.R"))

  one_imp <- function(i, estimate, variance) {
    data.table::data.table(imp = i, estimand = ltmle_estimands,
                           type = rep(c("mean", "ate"), c(4, 3)),
                           estimate = estimate, variance = variance)
  }
  tabs <- list(
    one_imp(1L, c(51.2, 47.8, 47.3, 42.9, -3.4, -3.9, -8.3),
            c(0.011, 0.039, 0.050, 0.063, 0.047, 0.057, 0.071)),
    one_imp(2L, c(51.5, 47.4, 47.7, 42.6, -4.1, -3.8, -8.9),
            c(0.012, 0.036, 0.054, 0.060, 0.044, 0.061, 0.075))
  )
  p <- pool_ltmle(tabs)

  expect_named(p, c("estimand", "type", "effect", "se", "ll", "ul"))
  expect_equal(p$estimand, ltmle_estimands)
  expect_equal(p$type, rep(c("mean", "ate"), c(4, 3)))
  purrr::walk(ltmle_estimands, \(k) {
    rr <- mice::pool.scalar(Q = purrr::map_dbl(tabs, \(t) t$estimate[t$estimand == k]),
                            U = purrr::map_dbl(tabs, \(t) t$variance[t$estimand == k]),
                            n = Inf, k = 1)
    expect_equal(p$effect[p$estimand == k], rr$qbar)
    expect_equal(p$se[p$estimand == k], sqrt(rr$t))
  })
  expect_true(all(p$ll < p$effect & p$effect < p$ul))

  # a bare table means the upstream target was not branched over imputations
  expect_error(pool_ltmle(tabs[[1]]), "single")
})

test_that("make_counterfactual_matrix follows the DAG and the ltmle adjustment set", {
  source(here::here("R", "helpers.R"))
  source(here::here("R", "confounders.R"))

  df    <- make_toy_2wave_df(n_per_cell = 10)
  pm    <- make_counterfactual_matrix(df)
  nodes <- c(wide_cols_mcs, "regime")
  expect_equal(dimnames(pm), list(nodes, nodes))

  # the g-formula outcome model adjusts for exactly the ltmle W, L1, A0 and A1
  expect_setequal(colnames(pm)[pm["sf12mcs_dv_1", ] == 1],
                  c(confounders("MCS"), time_varying_confounders("MCS"),
                    "econ_dist_bin_0", "econ_dist_bin_1"))
  # the start-wave outcome is not a confounder: no parents, and it predicts
  # nothing (the regime row is all ones by construction)
  expect_true(all(pm["sf12mcs_dv_0", ] == 0))
  expect_true(all(pm[setdiff(nodes, "regime"), "sf12mcs_dv_0"] == 0))
  # pcs_lagged_1 is measured with A0, so A0 does not cause it; final-wave income
  # responds to A0 and to its own past
  expect_equal(pm["pcs_lagged_1", "econ_dist_bin_0"], 0)
  expect_equal(pm["log_income_1", "econ_dist_bin_0"], 1)
  expect_equal(pm["log_income_1", "log_income_0"], 1)
  # gFormulaImpute() draws in column order with maxit = 1, so every imputed
  # node's parents must come earlier (strictly lower-triangular). Exposure rows
  # are exempt: set by the regime, never imputed.
  later <- upper.tri(pm, diag = TRUE)
  later[nodes %in% attr(df, "exposure_vars"), ] <- FALSE
  expect_true(all(pm[later] == 0))
  expect_true(all(pm["regime", setdiff(nodes, "regime")] == 1))
  expect_equal(pm["regime", "regime"], 0)
})

test_that("toy pipeline: mice -> ltmle -> pool -> gformula -> comparison", {
  source(here::here("R", "helpers.R"))
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "fit_ltmle_one.R"))
  source(here::here("R", "pool_ltmle.R"))
  source(here::here("R", "run_gformula.R"))
  source(here::here("R", "assemble_comparison.R"))

  df <- make_toy_2wave_df(n_per_cell = 150)
  # gFormulaImpute() takes M = mids$m, and syntheticPool() stops on a
  # non-positive total variance, which m = 2 hits often
  mids <- make_toy_2wave_mids(df, m = 20)

  ltmle_results <- purrr::map(1:2, \(i) {
    fit_ltmle_one(mids, i, "glm", "MCS", toy_regimes, toy_cv, 0.001)
  }) |> pool_ltmle()

  # gFormulaImpute()'s mice drops age_dv_0 from the age_dv_1 model (cor > .99)
  # and warns on every imputation; muffle only that warning
  gformula <- function(estimand) withCallingHandlers(
    run_gformula(mids, df, estimand, "MCS", toy_regimes),
    warning = \(w) if (startsWith(conditionMessage(w), "Number of logged events"))
      invokeRestart("muffleWarning")
  )
  set.seed(20261022)
  mi_results     <- gformula("factor(regime) + 0")
  mi_ate_results <- gformula("factor(regime)")

  cmp <- assemble_comparison(ltmle_results, mi_results, mi_ate_results)

  expect_named(cmp, c("marginal", "ate"))
  expect_equal(cmp$marginal$intervention, names(toy_regimes))
  expect_false(anyNA(cmp$marginal$tmle_effect))
  expect_false(anyNA(cmp$marginal$mi_effect))
  expect_equal(cmp$ate$estimand, c("before - never", "now - never", "permanent - never"))
  expect_false(anyNA(cmp$ate$tmle_effect))
  expect_false(anyNA(cmp$ate$mi_effect))

  # both land near the DGP's regime means: the test is about wiring and labels
  truth <- attr(df, "truth")[cmp$marginal$intervention]
  expect_lt(max(abs(cmp$marginal$tmle_effect - truth)), 1.5)
  expect_lt(max(abs(cmp$marginal$mi_effect - truth)), 2.5)
})

### MAIHDA layer: the 2 x 2 x 3 design from the stratum labels
test_that(".gate_factors recovers the intersectional design from strata_label", {
  source(here::here("R", "fit_gate_meta.R"))

  g <- data.frame(
    strata_label = c("Female:Non.white:High", "Female:Non.white:Medium",
                     "Female:Non.white:Low",  "Female:White:High",
                     "Female:White:Medium",   "Female:White:Low",
                     "Male:Non.white:High",   "Male:Non.white:Medium",
                     "Male:Non.white:Low",    "Male:White:High",
                     "Male:White:Medium",     "Male:White:Low"),
    estimate = 1:12, se = rep(1, 12)
  )

  out <- .gate_factors(g)
  expect_true(all(c("sex", "race", "educ") %in% names(out)))
  expect_equal(nlevels(out$sex), 2L)
  expect_equal(nlevels(out$race), 2L)
  expect_equal(nlevels(out$educ), 3L)
  expect_equal(nrow(unique(out[c("sex", "race", "educ")])), 12L)
  expect_equal(ncol(model.matrix(~ sex + race + educ, out)), 5L)

  # A label that does not split into exactly three parts must ERROR, not
  # silently mislabel cells.
  bad <- g; bad$strata_label[3] <- "Female:Non.white"
  expect_error(.gate_factors(bad), "three")

  # So must a set that is not a full 12-cell design.
  short <- g[1:6, ]
  expect_error(.gate_factors(short), "12")
})

test_that("fit_gate_meta returns documented tibbles and separates additive from interaction", {
  source(here::here("R", "rma_reml.R"))
  source(here::here("R", "fit_gate_meta.R"))

  labs <- as.vector(outer(
    outer(c("Female", "Male"), c("Non.white", "White"), paste, sep = ":"),
    c("High", "Medium", "Low"), paste, sep = ":"))
  base <- data.frame(imp = 1L, strata_id = as.character(seq_along(labs)),
                     strata_label = labs, n_j = 500L, pct_exposed = 20,
                     min_g = 0.025, median_g = 0.1, share_g_at_bound = 0.05,
                     se = rep(0.4, 12))
  g <- .gate_factors(base)
  X <- model.matrix(~ sex + race + educ, g)

  # (a) ADDITIVE truth plus modest noise -> main effects explain nearly
  # everything -> PCV high. The noise matters: with EXACTLY additive estimates
  # the additive model's RSS is 0, so s2w = 0, every SE collapses to 0 and the
  # deleted residuals become 0/0. Verified values for this fixture:
  # tau2_null 1.0169, tau2_main 0.0431, PCV 0.9576.
  set.seed(3)
  add <- base
  add$estimate <- as.vector(X %*% c(2, -0.3, 2.0, -0.2, -0.6)) + rnorm(12, 0, 0.5)
  r_add <- fit_gate_meta(add)

  expect_named(r_add, c("scalars", "coefs", "cells", "loco"))

  # column sets in exact order, as fit_gate_meta builds them
  expect_equal(nrow(r_add$scalars), 1L)
  expect_named(r_add$scalars,
               c("imp", "tau2_null", "se_tau2_null", "tau2_main", "se_tau2_main",
                 "pcv", "I2_null", "I2_main", "vt_null", "vt_main", "mu", "se_mu",
                 "QE_null", "QE_main", "QM", "QMp", "s2w_main"),
               ignore.order = FALSE)
  expect_type(r_add$scalars$imp, "integer")

  # coefs: intercept + 4 main effects
  expect_equal(nrow(r_add$coefs), 5L)
  expect_named(r_add$coefs, c("imp", "term", "estimate", "se"),
               ignore.order = FALSE)
  expect_type(r_add$coefs$imp, "integer")

  # cells: one row per stratum
  expect_equal(nrow(r_add$cells), 12L)
  expect_named(r_add$cells,
               c("imp", "strata_id", "strata_label", "estimate", "se",
                 "additive_pred", "se_pred", "blup", "se_blup", "del_resid",
                 "se_del_resid", "z", "tau2_del", "n_j", "pct_exposed", "min_g",
                 "median_g", "share_g_at_bound"),
               ignore.order = FALSE)
  expect_type(r_add$cells$imp, "integer")

  # loco: one row per deleted cell
  expect_equal(nrow(r_add$loco), 12L)
  expect_named(r_add$loco,
               c("imp", "dropped_id", "dropped_label", "tau2_null_d",
                 "se_tau2_null_d", "tau2_main_d", "se_tau2_main_d", "pcv_d"),
               ignore.order = FALSE)
  expect_type(r_add$loco$imp, "integer")

  # f0/f1 attribution: mu and QE_null must come from the null (intercept-only)
  # model; refit it independently to check
  f0_chk <- rma_reml(add$estimate, add$se^2, model.matrix(~ 1, g), knha = TRUE)
  expect_equal(r_add$scalars$mu, f0_chk$beta[1])
  expect_equal(r_add$scalars$se_mu, f0_chk$se[1])
  expect_equal(r_add$scalars$QE_null, f0_chk$QE)

  # null model must leave more unexplained heterogeneity than additive model
  expect_gt(r_add$scalars$QE_null, r_add$scalars$QE_main)

  expect_gt(r_add$scalars$pcv, 0.9)

  # (b) inject a large INTERACTION into one cell -> PCV collapses.
  # Verified for this fixture: PCV 0.9576 -> 0.3397.
  int <- add; int$estimate[9] <- int$estimate[9] - 6
  r_int <- fit_gate_meta(int)

  expect_named(r_int$scalars,
               c("imp", "tau2_null", "se_tau2_null", "tau2_main", "se_tau2_main",
                 "pcv", "I2_null", "I2_main", "vt_null", "vt_main", "mu", "se_mu",
                 "QE_null", "QE_main", "QM", "QMp", "s2w_main"),
               ignore.order = FALSE)
  expect_equal(nrow(r_int$loco), 12L)
  expect_type(r_int$loco$imp, "integer")

  expect_lt(r_int$scalars$pcv, 0.5)
  expect_lt(r_int$scalars$pcv, r_add$scalars$pcv)

  # the injected cell is the one the leave-one-out curve fingers
  worst <- r_int$loco$dropped_label[which.max(r_int$loco$pcv_d)]
  expect_equal(worst, labs[9])
})

### pool_gate_meta: Rubin's rules, documented shape
test_that("pool_gate_meta pools across imputations and forms PCV from pooled tau2", {
  source(here::here("R", "rma_reml.R"))
  source(here::here("R", "fit_gate_meta.R"))
  source(here::here("R", "pool_gate_meta.R"))

  labs <- as.vector(outer(
    outer(c("Female", "Male"), c("Non.white", "White"), paste, sep = ":"),
    c("High", "Medium", "Low"), paste, sep = ":"))
  mk <- function(i, jitter) {
    b <- data.frame(imp = as.integer(i), strata_id = as.character(seq_along(labs)),
                    strata_label = labs, n_j = 500L + i, pct_exposed = 20,
                    min_g = 0.025, median_g = 0.1, share_g_at_bound = 0.05,
                    se = rep(0.4, 12))
    g <- .gate_factors(b)
    X <- model.matrix(~ sex + race + educ, g)
    b$estimate <- as.vector(X %*% c(2, -0.3, 2.0, -0.2, -0.6)) + jitter
    b
  }
  # Jitter sd 0.8 (> cell SE 0.4) keeps every full-model tau2 strictly interior
  # (n_tau2_zero == 0), isolating Rubin pooling from the full-model boundary
  # case (tested below). Under this seed one of the 72 leave-one-cell-out
  # refits still lands on tau2 = 0: the LOCO-layer pattern Fix 1a added
  # counting for, also seen in the real data (LOCO clips the boundary, the
  # full-model fits stay interior). The LOCO block pools arithmetically (see
  # pool_gate_meta.R), so the assertions below are unaffected; the warning is
  # suppressed here and covered by the next test.
  # Verified: pooled tau2_null 1.987, tau2_main 0.365, PCV 0.816, se_mu ratio
  # 1.135.
  set.seed(7)
  one <- purrr::map(1:3, \(i) fit_gate_meta(mk(i, rnorm(12, 0, 0.8))))

  p <- suppressWarnings(pool_gate_meta(one))

  expect_named(p, c("scalars", "coefs", "cells", "loco", "tests", "m",
                    "n_tau2_zero", "n_tau2_zero_loco"))
  expect_equal(p$n_tau2_zero, 0L)
  expect_equal(p$m, 3L)
  expect_equal(nrow(p$scalars), 1L)
  expect_equal(nrow(p$cells), 12L)
  expect_equal(nrow(p$loco), 12L)
  expect_equal(nrow(p$coefs), 5L)
  expect_equal(nrow(p$tests), 3L)          # one row per imputation, not pooled

  expect_true(is.finite(p$scalars$tau2_null))
  expect_true(is.finite(p$scalars$tau2_main))
  expect_gte(p$scalars$tau2_null, 0)
  expect_gte(p$scalars$tau2_main, 0)

  # Pinned values (Fix 3), promoted from the verified comment above: a
  # back-transform regression (e.g. dropping "- eps", or se_tau2/tau2 instead
  # of se_tau2/(tau2 + eps)) fails loudly instead of sliding through on
  # finiteness/shape checks alone.
  expect_equal(p$scalars$tau2_null, 1.987, tolerance = 1e-3)
  expect_equal(p$scalars$tau2_main, 0.365, tolerance = 1e-3)
  expect_equal(p$scalars$pcv,       0.816, tolerance = 1e-3)

  # PCV must be formed FROM the pooled tau2 values, not averaged across
  # imputations -- individual-imputation PCV can be negative at J = 12.
  expect_equal(p$scalars$pcv,
               (p$scalars$tau2_null - p$scalars$tau2_main) / p$scalars$tau2_null)

  # pooled SEs exceed the mean within-imputation SE (between-imputation variance)
  mean_within <- mean(purrr::map_dbl(one, ~ .x$scalars$se_mu))
  expect_gte(p$scalars$se_mu, mean_within * 0.999)
  expect_equal(p$scalars$se_mu / mean_within, 1.135, tolerance = 1e-3)

  expect_true(all(p$cells$n_j_min <= p$cells$n_j_max))
})

### log-scale pooling can collapse tau2: pool_gate_meta must warn
test_that("pool_gate_meta warns when any imputation hits the tau2 = 0 boundary", {
  source(here::here("R", "rma_reml.R"))
  source(here::here("R", "fit_gate_meta.R"))
  source(here::here("R", "pool_gate_meta.R"))

  labs <- as.vector(outer(
    outer(c("Female", "Male"), c("Non.white", "White"), paste, sep = ":"),
    c("High", "Medium", "Low"), paste, sep = ":"))
  mk <- function(i, sd_noise) {
    b <- data.frame(imp = as.integer(i), strata_id = as.character(seq_along(labs)),
                    strata_label = labs, n_j = 500L, pct_exposed = 20,
                    min_g = 0.025, median_g = 0.1, share_g_at_bound = 0.05,
                    se = rep(0.4, 12))
    g <- .gate_factors(b)
    X <- model.matrix(~ sex + race + educ, g)
    set.seed(100 + i)
    b$estimate <- as.vector(X %*% c(2, -0.3, 2.0, -0.2, -0.6)) + rnorm(12, 0, sd_noise)
    b
  }
  # sd well below the cell SE of 0.4 forces tau2_main to the zero boundary
  one <- purrr::map(1:3, \(i) fit_gate_meta(mk(i, 0.05)))
  expect_gt(sum(purrr::map_dbl(one, ~ .x$scalars$tau2_main) == 0), 0)

  expect_warning(pool_gate_meta(one), "boundary")

  # the non-negativity invariant still holds after the warning
  p <- suppressWarnings(pool_gate_meta(one))
  expect_gte(p$scalars$tau2_null, 0)
  expect_gte(p$scalars$tau2_main, 0)
  expect_gte(p$scalars$tau2_null_ll, 0)
  expect_gte(p$scalars$tau2_main_ll, 0)
  expect_gt(p$n_tau2_zero, 0)

  # Fix 4: Fix 1's defect (LOCO layer unguarded, only the full-model layer
  # counted) survived review because this test asserted only on p$scalars.
  # Assert on p$loco too.
  expect_true("n_zero_d" %in% names(p$loco))
  expect_gt(p$n_tau2_zero_loco, 0)
  expect_true(any(p$loco$n_zero_d > 0))
  # "either finite or NA, never NaN or Inf": is.na(NaN) is TRUE in R, so check
  # is.nan()/is.infinite() explicitly rather than relying on is.finite()/is.na().
  expect_true(all(!is.nan(p$loco$pcv_d) & !is.infinite(p$loco$pcv_d)))
})
