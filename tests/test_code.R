# testing files of the pipeline and specific functions
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
  sandwich)

## testing functions
### does build_data not crash and produce the expected single-wave wide columns
test_that("build_data produces expected output", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))

  pop_data <- import_data(force = FALSE) |> clean_data() |> preproc_data()

  wide_data <- build_data(pop_data, outcome = "MCS")
  expect_true(is.data.frame(wide_data))

  # single target wave: everything is _base or _0 (codebase-map 1)
  expected_cols <- c(
    "sex_dv_base", "hiqual_dv_fact_base", "race_base", "sf12mcs_dv_base",
    "sf12mcs_dv_0", "age_dv_0", "gor_dv_fact_0",
    "pcs_lagged_0", "econ_dist_bin_0", "econ_dist_bin_lagged_0",
    "dnc_lagged_0", "home_owner_lagged_0", "econ_benefits_lagged_0",
    "mastat_lagged_0", "econ_emp_bin_fact_0", "log_income_0"
  )
  expect_true(all(expected_cols %in% names(wide_data)))
})

### make_counterfactual_method: methods gFormulaImpute() draws the synthetic block with.
### Synthetic frame, so this runs without the raw .dta files.
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

### does run_mice not crash and generate mids object with expected variables
#### also if new make_counterfactual_matrix works fine
test_that("run_mice produces expected output", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))

  pop_data  <- import_data(force = FALSE) |> clean_data() |> preproc_data()
  wide_data <- build_data(pop_data)

  # setting exposure
  wide_data <- set_exposure(wide_data, exposure = "econ_dist_bin")

  # counterfactual predictor matrix
  predictor_matrix <- make_counterfactual_matrix(wide_data)
  
  # mini-mice call to check structure only
  mids <- run_mice(wide_data, m = 2, maxit = 2, seed = 42)

  # returns a proper mids object with the requested number of imputations
  expect_s3_class(mids, "mids")
  expect_equal(mids$m, 2L)

  # dianosing mice imputation, only check these results with at least 15 iterations
  source(here::here("tests", "diagnose_mice.R"))
  diagnose_mice(mids, wide_data,
                plot_dir = here::here("tests", "figs", "wide", "mice_diag"))

  # mice ran without logged events (warnings or errors during imputation)
  expect_null(mids$loggedEvents)

  # does counterfactual matrix have regime row and column
  expect_true("regime" %in% rownames(predictor_matrix))
  
  # are all regime row values 1 except for the regime column
  off_diag <- predictor_matrix["regime", colnames(predictor_matrix) != "regime"]
  expect_true(all(off_diag == 1))
  expect_equal(unname(predictor_matrix["regime", "regime"]), 0)
})

### does fit_tmle_one not crash
test_that("fit_tmle_one runs without error", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "fit_tmle_one.R"))

  pop_data  <- import_data(force = FALSE) |> clean_data() |> preproc_data()
  wide_data <- build_data(pop_data, outcome = "MCS")

  mids <- run_mice(wide_data, m = 2, maxit = 2, seed = 42)

  tmle_fit <- fit_tmle_one(mids,
                           imp_idx = 1,
                           sl_libs = c("SL.mean", "SL.glm"),
                           outcome = "MCS")

  # fit_tmle_one calls tmle::tmle(), which returns a "tmle" object
  expect_s3_class(tmle_fit, "tmle")
  expect_true(is.numeric(tmle_fit$estimates$ATE$psi))
})

## can i add an additional variable

test_that("will a new variable exist after imputation", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))

  pop_data  <- import_data(force = FALSE) |> clean_data() |> preproc_data()
  wide_data <- build_data(pop_data)

  # setting exposure
  wide_data <- set_exposure(wide_data, exposure = "econ_dist_bin")

  # counterfactual predictor matrix
  predictor_matrix <- make_counterfactual_matrix(wide_data)
  
  # mini-mice call to check structure only
  mids <- run_mice(wide_data, m = 2, maxit = 2, seed = 42)

  # returns a proper mids object with the requested number of imputations
  expect_s3_class(mids, "mids")
  expect_equal(mids$m, 2L)

  # adding new var
  as_long_mids <- mice::complete(mids, action = "long", include = TRUE)
  
  as_long_mids <- as_long_mids |> mutate(new_var = paste0("var_", as.character(sex_dv_base)))

  imp_new <- mice::as.mids(as_long_mids)
  
  # mice ran without logged events (warnings or errors during imputation)
  expect_null(mids$loggedEvents)

  # does counterfactual matrix have regime row and column
  expect_true("regime" %in% rownames(predictor_matrix))
  
  # are all regime row values 1 except for the regime column
  off_diag <- predictor_matrix["regime", colnames(predictor_matrix) != "regime"]
  expect_true(all(off_diag == 1))
  expect_equal(unname(predictor_matrix["regime", "regime"]), 0)
})


## verifying how many strata will be for MAIHDA
test_that("how many strata will be for MAIHDA", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))

  pop_data  <- import_data(force = FALSE) |> clean_data() |> preproc_data()
  wide_data <- build_data(pop_data, pre_wave = 2, target_wave = 3, outcome = "MCS")

  # mini-mice call to check structure only
  mids <- run_mice(wide_data, m = 5, maxit = 2, seed = 42)

  # counterfactual predictor matrix
  predictor_matrix <- make_counterfactual_matrix(wide_data)

  # returns a proper mids object with the requested number of imputations
  expect_s3_class(mids, "mids")

  # mice ran without logged events (warnings or errors during imputation)
  expect_null(mids$loggedEvents)

  # how many strata will be for MAIHDA in observed outcome data
  ## creating strata_id per each imputed dataset
  fits <- purrr::map(seq_len(mids$m), function(i) {
    dat_i <- mice::complete(mids, action = i) |>
      dplyr::mutate(
        # intersectional strata id: unique cell
        strata_id = 100 * as.numeric(sex_dv_base) +
                     10 * as.numeric(race_base) +
                          as.numeric(hiqual_dv_fact_base),
        strata_id = as.factor(strata_id),
        # human-readable stratum name, e.g. "female:white:high"
        strata_label = interaction(sex_dv_base, race_base, hiqual_dv_fact_base,
                                   sep = ":", drop = TRUE)
      ) }
    )
  # counts per stratum (id + label) for each imputed dataset
  strata_tables <- purrr::map(fits, function(dat_i) {
    dplyr::count(dat_i, strata_id, strata_label)
  })

})

## ---- DR-learner GATE arm (design 2026-07-30) --------------------------------

### shared confounder set: byte-identical between fit_tmle_one and estimate_cate
test_that("confounders returns the exact fit_tmle_one adjustment set", {
  source(here::here("R", "confounders.R"))

  expect_identical(
    confounders("MCS"),
    c("sex_dv_base",
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

  # manual V-fold refit reproduces Z exactly, in input row order (the design's
  # empirical verification of 3.6, kept as a regression test against
  # SuperLearner version changes)
  manual <- rep(NA_real_, n)
  for (v in seq_along(fit$validRows)) {
    rows <- fit$validRows[[v]]
    gfit <- glm(y ~ ., data = data.frame(y = Yc[-rows], X[-rows, ]),
                family = gaussian())
    manual[rows] <- predict(gfit, newdata = X[rows, ], type = "response")
  }
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

### toy generators for the DR-learner tests (used by the estimate_cate,
### pool_cate and congenial run_mice tests; defined once, here only)
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
  a <- unlist(lapply(seq_len(12), function(j) {
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

### design 7 "factor handling": estimate_cate runs on a small mids whose W
### contains factors
### Distinct per-cell tau (design review #1a): a same-constant-everywhere
### tau_by_cell (e.g. the generator's own default rep(-2, 12)) can't tell a
### correct sign from a flipped one, because mean(psi) == n_j-weighted mean
### of the cell GATEs is a pure algebraic identity that holds for *any* psi,
### negated included (mutation-tested by the design-2026-07-30 reviewer: a
### negated-psi mutant and a mu0/mu1-arm-swap mutant both left every original
### assertion here passing). Distinct, signed effects per cell plus an
### explicit sign/recovery check against them close that gap.
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
  # A sign-flip mutant misses this by ~2x|tau| (>= 2 everywhere here);
  # measured max error under the true code is 0.88 (no-imputation-noise floor
  # ~0.36, through make_toy_wide_mids() ~0.70), so tolerance = 1.0 is safe
  # margin without being vacuous.
  expect_equal(res$gate$estimate, -tau_by_cell, tolerance = 1.0)
  # beta2_hte pins near +1 when the tau proxy tracks the true heterogeneity on
  # the conventional scale; it flips to about -1 under a sign error.
  expect_gt(res$blp$estimate[res$blp$term == "beta2_hte"], 0.5)
})

### design review #1b: make_toy_wide_df() randomizes A at a fixed 40% inside
### every cell, so ghat is exactly recoverable and AIPW is then invariant to
### which fitted model plays mu0's role vs mu1's role (a real, not
### coincidental, consequence of double robustness when g is exactly
### correctly specified -- confirmed algebraically and empirically: an
### arm-swap mutant left the test above passing). Catching an arm swap needs
### (i) A confounded by a continuous W the tiny SL.glm/SL.mean library can't
### capture (so ghat is NOT exactly right) and (ii) a genuine W x A outcome
### interaction (so mu0/mu1 are different functions of W, not a location
### shift) -- otherwise the swap has nothing to bite on. This is a SEPARATE
### generator, not a change to make_toy_wide_df(), which Tasks 6/12 rely on
### unchanged.
# A depends on age through a quadratic logit (b1*age_c + b2*age_c^2): a
# linear-only learner (SL.glm) cannot recover this, so ghat is persistently
# biased regardless of sample size. Y's treatment effect is linear in age
# (tau0 + tau_age*(age-45)), a real W x A interaction SL.glm CAN represent for
# mu0/mu1, isolating the arm-assignment bug as the swap's only source of
# error. Calibrated empirically: swapping mu0_hat/mu1_hat inflates
# mean(psi) by ~+0.9 to +1.1 against the truth, while the correctly-assigned
# estimator recovers it within ~0.15.
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

  # true population ATE on the Y(0)-Y(1) scale, averaged over this sample's
  # realized ages (tau0/tau_age fully determine it; no per-cell truth here --
  # age is drawn iid of strata, so all 12 cells share the same population
  # truth up to sampling noise)
  truth <- -mean(tau0 + tau_age * (df$age_dv_0 - 45))
  # expect_equal(..., tolerance=) on a *scalar* is a RELATIVE tolerance
  # (all.equal.numeric semantics), not absolute -- confirmed empirically it
  # would silently pass a mutant this test exists to catch. Assert the
  # absolute gap directly instead: correctly-assigned arms recover the truth
  # within ~0.05-0.15 at this n/seed; a mu0/mu1 arm swap inflates it by
  # ~+1.0 (measured). 0.5 sits well inside that gap.
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

  fits   <- lapply(1:2, function(i) {
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

### review finding: a bare factor(new_labels[...]) with no `levels =` derives
### its level set from observed rows only, so a defined-but-unused level
### silently vanishes instead of being relabelled and kept -- breaking the
### "levels replaced" contract. Guards the `levels = sort(unique(new_labels))`
### fix.
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
### (it, im, dep) into a single comma-joined `out`. Anything checking it
### programmatically has to re-parse it per term, and has to cope with the two
### exits where remove.lindep names no term at all.

# minimal one-row loggedEvents in mice's own column layout
logged_event <- function(out, dep = "sf12mcs_dv_0", meth = "collinear") {
  list(loggedEvents = data.frame(it = 1L, im = 1L, dep = dep, meth = meth,
                                 out = out, stringsAsFactors = FALSE))
}

test_that("assert_congenial_terms requires both node names", {
  source(here::here("R", "run_mice.R"))

  # No defaults: a_col/y_col are always resolved by the caller, so a default
  # exists only to be silently wrong on a wave whose suffix is not _0.
  expect_error(
    assert_congenial_terms(logged_event("econ_dist_bin_01:race_baseNon.white"),
                           a_col = "econ_dist_bin_0"),
    "y_col"
  )
  expect_error(
    assert_congenial_terms(logged_event("econ_dist_bin_01:race_baseNon.white"),
                           y_col = "sf12mcs_dv_0"),
    "a_col"
  )
})

test_that("assert_congenial_terms stops when the Y model drops an A:S interaction", {
  source(here::here("R", "run_mice.R"))

  err <- expect_error(
    assert_congenial_terms(logged_event("econ_dist_bin_01:race_baseNon.white",
                                        dep = "sf12mcs_dv_0"),
                           a_col = "econ_dist_bin_0", y_col = "sf12mcs_dv_0"),
    "congenial"
  )
  # the dropped term must appear in the message -- deciding rung 1 vs rung 2
  # needs to know *which* term went
  expect_true(grepl("econ_dist_bin_01:race_baseNon.white",
                    conditionMessage(err), fixed = TRUE))
  expect_true(grepl("1 exposure-interaction term(s)",
                    conditionMessage(err), fixed = TRUE))
})

test_that("assert_congenial_terms stops when the A model drops a Y:S interaction", {
  source(here::here("R", "run_mice.R"))

  # run_mice puts y_col * S in the A model, so a break there is logged under the
  # OUTCOME name with dep = the exposure. Matching only on a_col misses it.
  expect_error(
    assert_congenial_terms(logged_event("sf12mcs_dv_0:sex_dv_baseMale",
                                        dep = "econ_dist_bin_0"),
                           a_col = "econ_dist_bin_0", y_col = "sf12mcs_dv_0"),
    "congenial"
  )
})

test_that("assert_congenial_terms counts dropped terms, not logged rows", {
  source(here::here("R", "run_mice.R"))

  err <- expect_error(
    assert_congenial_terms(
      logged_event(paste("econ_dist_bin_01:sex_dv_baseMale",
                         "econ_dist_bin_01:race_baseNon.white",
                         "econ_dist_bin_01:hiqual_dv_fact_baseLow",
                         sep = ", ")),
      a_col = "econ_dist_bin_0", y_col = "sf12mcs_dv_0"),
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
      a_col = "econ_dist_bin_0", y_col = "sf12mcs_dv_0"),
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
      logged_event("econ_dist_bin_0, sex_dv_baseMale:race_baseNon.white",
                   dep = "log_income_0"),
      a_col = "econ_dist_bin_0", y_col = "sf12mcs_dv_0")
  )
})

test_that("assert_congenial_terms matches on the node names it is given, not on wave 0", {
  source(here::here("R", "run_mice.R"))

  # target_wave = 7 -> suffix _4. The guard must follow the names it is handed.
  wave4 <- logged_event("econ_dist_bin_41:race_baseNon.white",
                        dep = "sf12mcs_dv_4")
  expect_error(
    assert_congenial_terms(wave4, a_col = "econ_dist_bin_4", y_col = "sf12mcs_dv_4"),
    "congenial"
  )
  # PCS resolves to a different outcome stem, same wave
  expect_error(
    assert_congenial_terms(logged_event("sf12pcs_dv_4:sex_dv_baseMale",
                                        dep = "econ_dist_bin_4"),
                           a_col = "econ_dist_bin_4", y_col = "sf12pcs_dv_4"),
    "congenial"
  )
})

test_that("assert_congenial_terms passes on clean or unrelated logged events", {
  source(here::here("R", "run_mice.R"))

  A <- "econ_dist_bin_0"; Y <- "sf12mcs_dv_0"

  expect_true(assert_congenial_terms(list(loggedEvents = NULL), A, Y))
  expect_true(assert_congenial_terms(
    logged_event("gor_dv_fact_0", dep = "log_income_0", meth = "pmm"), A, Y))
  # a main-effect drop of the exposure itself is not an interaction drop
  expect_true(assert_congenial_terms(
    logged_event("econ_dist_bin_lagged_0", meth = "pmm"), A, Y))
  # ...and neither is a lagged-exposure interaction: econ_dist_bin_lagged_0 is a
  # different node from the A node, even though it shares the stem
  expect_true(assert_congenial_terms(
    logged_event("econ_dist_bin_lagged_01:gor_dv_fact_0R2", meth = "pmm"), A, Y))
})

### run_mice: formulas carry the saturated A x S interaction
test_that("run_mice builds Y and A formulas with the saturated strata interaction", {
  source(here::here("R", "build_data.R"))
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "import_cleaning.R"))
  source(here::here("R", "helpers.R"))

  wide_data <- build_data(preproc_data(clean_data(import_data(force = TRUE))),
                          outcome = "MCS")

  # m = 1, maxit = 1: we are checking model structure, not convergence
  mids <- run_mice(wide_data, m = 1, maxit = 1, seed = 42, outcome = "MCS")

  expect_s3_class(mids, "mids")

  y_terms <- attr(terms(mids$formulas[["sf12mcs_dv_0"]]), "term.labels")
  a_terms <- attr(terms(mids$formulas[["econ_dist_bin_0"]]), "term.labels")

  # Term-label variable ORDER in terms() follows each variable's main-effect
  # position in the whole formula (set by make.formulas()'s column order), not
  # the order written in the interaction expression passed to reformulate() --
  # so "econ_dist_bin_0:sex_dv_base:race_base:hiqual_dv_fact_base" comes back
  # as "sex_dv_base:hiqual_dv_fact_base:race_base:econ_dist_bin_0" on real
  # data. Same interaction (model.matrix/mice do not care about label order),
  # so match on the *set* of variables in a term rather than an exact string.
  has_term <- function(term_labels, vars) {
    any(vapply(strsplit(term_labels, ":", fixed = TRUE),
               function(parts) setequal(parts, vars), logical(1)))
  }

  # Y model: the 4-way A x S term must be present
  expect_true(has_term(y_terms, c("econ_dist_bin_0", "sex_dv_base", "race_base", "hiqual_dv_fact_base")))
  # A model: saturated S, but no A term (A is the target there)
  expect_true(has_term(a_terms, c("sex_dv_base", "race_base", "hiqual_dv_fact_base")))
  expect_false(any(grepl("econ_dist_bin_0", a_terms, fixed = TRUE)))
})

### does the MAIHDA layer recover the 2 x 2 x 3 design from the stratum labels
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

### does fit_gate_meta return the documented shape, and does PCV behave
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

  # top-level structure: four named components
  expect_named(r_add, c("scalars", "coefs", "cells", "loco"))

  # scalars: full column set in exact order (as constructed in fit_gate_meta)
  expect_equal(nrow(r_add$scalars), 1L)
  expect_named(r_add$scalars,
               c("imp", "tau2_null", "se_tau2_null", "tau2_main", "se_tau2_main",
                 "pcv", "I2_null", "I2_main", "vt_null", "vt_main", "mu", "se_mu",
                 "QE_null", "QE_main", "QM", "QMp", "s2w_main"),
               ignore.order = FALSE)
  expect_type(r_add$scalars$imp, "integer")

  # coefs: full column set, 5 rows (intercept + 4 main effects)
  expect_equal(nrow(r_add$coefs), 5L)
  expect_named(r_add$coefs, c("imp", "term", "estimate", "se"),
               ignore.order = FALSE)
  expect_type(r_add$coefs$imp, "integer")

  # cells: full column set, 12 rows (one per stratum)
  expect_equal(nrow(r_add$cells), 12L)
  expect_named(r_add$cells,
               c("imp", "strata_id", "strata_label", "estimate", "se",
                 "additive_pred", "se_pred", "blup", "se_blup", "del_resid",
                 "se_del_resid", "z", "tau2_del", "n_j", "pct_exposed", "min_g",
                 "median_g", "share_g_at_bound"),
               ignore.order = FALSE)
  expect_type(r_add$cells$imp, "integer")

  # loco: full column set, 12 rows (one per deleted cell)
  expect_equal(nrow(r_add$loco), 12L)
  expect_named(r_add$loco,
               c("imp", "dropped_id", "dropped_label", "tau2_null_d",
                 "se_tau2_null_d", "tau2_main_d", "se_tau2_main_d", "pcv_d"),
               ignore.order = FALSE)
  expect_type(r_add$loco$imp, "integer")

  # f0/f1 attribution: mu and QE_null must come from the null (intercept-only) model
  # Fit independently to verify
  f0_chk <- rma_reml(add$estimate, add$se^2, model.matrix(~ 1, g), knha = TRUE)
  expect_equal(r_add$scalars$mu, f0_chk$beta[1])
  expect_equal(r_add$scalars$se_mu, f0_chk$se[1])
  expect_equal(r_add$scalars$QE_null, f0_chk$QE)

  # null model must leave more unexplained heterogeneity than additive model
  expect_gt(r_add$scalars$QE_null, r_add$scalars$QE_main)

  # PCV is high in the additive case
  expect_gt(r_add$scalars$pcv, 0.9)

  # (b) inject a large INTERACTION into one cell -> PCV collapses.
  # Verified for this fixture: PCV 0.9576 -> 0.3397.
  int <- add; int$estimate[9] <- int$estimate[9] - 6
  r_int <- fit_gate_meta(int)

  # same interface checks on the interaction case
  expect_named(r_int$scalars,
               c("imp", "tau2_null", "se_tau2_null", "tau2_main", "se_tau2_main",
                 "pcv", "I2_null", "I2_main", "vt_null", "vt_main", "mu", "se_mu",
                 "QE_null", "QE_main", "QM", "QMp", "s2w_main"),
               ignore.order = FALSE)
  expect_equal(nrow(r_int$loco), 12L)
  expect_type(r_int$loco$imp, "integer")

  # interaction behavior: PCV collapses
  expect_lt(r_int$scalars$pcv, 0.5)
  expect_lt(r_int$scalars$pcv, r_add$scalars$pcv)

  # the injected cell is the one the leave-one-out curve fingers
  worst <- r_int$loco$dropped_label[which.max(r_int$loco$pcv_d)]
  expect_equal(worst, labs[9])
})

### does pool_gate_meta apply Rubin's rules in the documented shape
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
  # Jitter sd 0.8 (> the cell SE of 0.4) keeps every full-model tau2 strictly
  # interior (n_tau2_zero == 0), so this fixture isolates Rubin pooling from
  # the full-model boundary case (tested separately below). Under this seed
  # one of the 72 leave-one-cell-out refits still lands on the tau2 = 0
  # boundary -- exactly the LOCO-layer pattern Fix 1a added counting for, and
  # the same pattern seen in the real data (LOCO clips the boundary while the
  # full-model fits stay interior). The LOCO block pools arithmetically (see
  # pool_gate_meta.R) so this does not affect the assertions below; the
  # resulting warning is suppressed here and covered by the next test.
  # Verified: pooled tau2_null 1.987, tau2_main 0.365, PCV 0.816, se_mu ratio
  # 1.135.
  set.seed(7)
  one <- lapply(1:3, function(i) fit_gate_meta(mk(i, rnorm(12, 0, 0.8))))

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

  # Pinned values (Fix 3): promoted from a verified-but-unasserted comment so
  # a regression in the back-transform (e.g. dropping "- eps", or writing
  # se_tau2/tau2 instead of se_tau2/(tau2 + eps)) fails loudly instead of
  # sliding through on finiteness/shape checks alone.
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

### does pool_gate_meta warn when log-scale pooling can collapse tau2
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
  one <- lapply(1:3, function(i) fit_gate_meta(mk(i, 0.05)))
  expect_gt(sum(purrr::map_dbl(one, ~ .x$scalars$tau2_main) == 0), 0)

  expect_warning(pool_gate_meta(one), "boundary")

  # After warning fires, also verify the non-negativity invariant holds
  p <- suppressWarnings(pool_gate_meta(one))
  expect_gte(p$scalars$tau2_null, 0)
  expect_gte(p$scalars$tau2_main, 0)
  expect_gte(p$scalars$tau2_null_ll, 0)
  expect_gte(p$scalars$tau2_main_ll, 0)
  expect_gt(p$n_tau2_zero, 0)

  # Fix 4: this is structurally why Fix 1's defect (the LOCO layer was
  # unguarded while only the full-model layer was counted) survived review --
  # the boundary test asserted only on p$scalars. Assert on p$loco too.
  expect_true("n_zero_d" %in% names(p$loco))
  expect_gt(p$n_tau2_zero_loco, 0)
  expect_true(any(p$loco$n_zero_d > 0))
  # "either finite or NA, never NaN or Inf": is.na(NaN) is TRUE in R, so check
  # is.nan()/is.infinite() explicitly rather than relying on is.finite()/is.na().
  expect_true(all(!is.nan(p$loco$pcv_d) & !is.infinite(p$loco$pcv_d)))
})
