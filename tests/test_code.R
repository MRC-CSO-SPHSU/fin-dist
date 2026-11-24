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
# ltmle_regimes, sl_cv_control and ltmle_bound from _targets.R: c(A0, A1), first = reference
toy_regimes <- list(never = c(0L, 0L), before = c(1L, 0L),
                    now = c(0L, 1L), permanent = c(1L, 1L))
toy_cv <- list(V = 5L, stratifyCV = FALSE, shuffle = TRUE)
# fit_ltmle_one() passes gbound straight to ltmle's gbounds: c(lower, upper) on the cumulative g
toy_gbound <- c(0.001, 0.999)
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
  expect_true(all(confounders("MCS") %in% names(wide_data)))
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
                       gbound        = toy_gbound)

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

### shared confounder set: W (base + start wave), then L1 (final wave)
test_that("confounders returns the exact ltmle adjustment set", {
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
      "log_income_0",
      "age_dv_1",
      "gor_dv_fact_1",
      "pcs_lagged_1",
      "dnc_lagged_1",
      "home_owner_lagged_1",
      "econ_benefits_lagged_1",
      "mastat_lagged_1",
      "econ_emp_bin_fact_1",
      "log_income_1")
  )
  # PCS swaps exactly the outcome-dependent entries
  expect_identical(setdiff(confounders("PCS"), confounders("MCS")),
                   c("sf12pcs_dv_base", "mcs_lagged_0", "mcs_lagged_1"))
  expect_identical(setdiff(confounders("MCS"), confounders("PCS")),
                   c("sf12mcs_dv_base", "pcs_lagged_0", "pcs_lagged_1"))
  # econ_dist_bin_lagged_1 is A0 itself: never a confounder
  expect_false("econ_dist_bin_lagged_1" %in% c(confounders("MCS"), confounders("PCS")))
  expect_error(confounders("XXX"))
})

### the W / L1 split fit_ltmle_one() and estimate_cate() take by the _1 suffix
test_that("confounders splits into W before A0 and L1 between A0 and A1", {
  source(here::here("R", "confounders.R"))

  purrr::walk(c("MCS", "PCS"), \(o) {
    x <- confounders(o)
    W <- grep("_1$", x, value = TRUE, invert = TRUE)
    L <- grep("_1$", x, value = TRUE)
    # W: base and start-wave columns only, in confounders() order
    expect_true(all(grepl("_(base|0)$", W)))
    expect_length(W, 14L)
    # L1: the start-wave confounders measured again at the final wave, same order,
    # minus the lagged exposure (its _1 copy is A0)
    expect_identical(L, sub("_0$", "_1", setdiff(grep("_0$", W, value = TRUE),
                                                 "econ_dist_bin_lagged_0")))
  })
})

### make_strata: pure df -> df, same id arithmetic as the old strata_creation()
test_that("make_strata builds the 12-cell id and label", {
  source(here::here("R", "strata_creation.R"))

  set.seed(1)
  toy <- data.frame(
    sex_dv_fact_base    = factor(sample(c("Female", "Male"), 240, TRUE)),
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
               as.character(100 * as.numeric(toy$sex_dv_fact_base) +
                             10 * as.numeric(toy$race_base) +
                                  as.numeric(toy$hiqual_dv_fact_base)))
  # row order and other columns untouched
  expect_equal(out$other_col, toy$other_col)
  # label format sex:race:educ
  expect_true(all(grepl("^[^:]+:[^:]+:[^:]+$", as.character(out$strata_label))))
})

test_that("make_strata errors when a stratum column is absent", {
  source(here::here("R", "strata_creation.R"))
  expect_error(make_strata(data.frame(sex_dv_fact_base = factor("Female"))), "missing")
})

### Two-time-point DR-learner primitives
# Closed-form two-time-point DGP: W -> A0 -> L1 -> A1 -> Y, L1 confounds A1.
# Q2 = E[Y | W, A0, L1, A1], Q1 = E[Q2 | W, A0] and E[Y^a] are known exactly.
make_toy_seq_truth <- function(n = 1e5, seed = 1) {
  set.seed(seed)
  W  <- rnorm(n)
  A0 <- rbinom(n, 1, plogis(-0.5 + 0.8 * W))
  L1 <- 0.5 * W - 0.7 * A0 + rnorm(n)
  A1 <- rbinom(n, 1, plogis(-1 + 0.6 * L1 + 0.8 * A0))
  Y  <- 50 + 2 * W + 1.5 * L1 - 2 * A0 - 4 * A1 + rnorm(n, 0, 2)
  list(Y = Y, A0 = A0, A1 = A1,
       Q2 = \(a0, a1) 50 + 2 * W + 1.5 * L1 - 2 * a0 - 4 * a1,
       Q1 = \(a0, a1) 50 + 2.75 * W - 3.05 * a0 - 4 * a1,
       g0 = \(a0) { p <- plogis(-0.5 + 0.8 * W); if (a0 == 1) p else 1 - p },
       g1 = \(a0, a1) { p <- plogis(-1 + 0.6 * L1 + 0.8 * a0); if (a1 == 1) p else 1 - p },
       truth = \(a0, a1) 50 - 3.05 * a0 - 4 * a1)
}

test_that("seq_dr_pseudo_outcome recovers E[Y^a] and is doubly robust", {
  source(here::here("R", "meta_learner.R"))
  d <- make_toy_seq_truth()
  n <- length(d$Y)
  res <- purrr::map(toy_regimes, \(r) {
    r0  <- d$A0 == r[1]
    r1  <- r0 & d$A1 == r[2]
    pi0 <- d$g0(r[1])
    pi1 <- pi0 * d$g1(r[1], r[2])
    phi <- \(Q1, Q2, p0, p1) seq_dr_pseudo_outcome(r0, r1, d$Y, Q1, Q2, p0, p1)
    est <- list(
      both   = phi(d$Q1(r[1], r[2]), d$Q2(r[1], r[2]), pi0, pi1),
      g_only = phi(rep(50, n), rep(50, n), pi0, pi1),
      Q_only = phi(d$Q1(r[1], r[2]), d$Q2(r[1], r[2]), rep(0.5, n), rep(0.25, n)),
      none   = phi(rep(50, n), rep(50, n), rep(0.5, n), rep(0.25, n))
    )
    purrr::map(est, \(p) c(err = mean(p) - d$truth(r[1], r[2]), se = sd(p) / sqrt(n)))
  })
  purrr::walk(res, \(b) purrr::walk(b[c("both", "g_only", "Q_only")],
                                    \(x) expect_lt(abs(x[["err"]]), 4 * x[["se"]])))
  # the checks bite: with both nuisances wrong some regime is far off
  expect_gt(max(purrr::map_dbl(res, \(b) abs(b$none[["err"]]))), 0.5)
})

test_that("seq_dr_pseudo_outcome refuses NA, unbounded pi and r1 without r0", {
  source(here::here("R", "meta_learner.R"))
  ok <- list(r0 = c(1L, 0L), r1 = c(1L, 0L), Y = c(1, 2), Q1 = c(1, 1), Q2 = c(1, 1),
             pi0 = c(0.5, 0.5), pi1 = c(0.5, 0.5))
  call_with <- \(...) do.call(seq_dr_pseudo_outcome, utils::modifyList(ok, list(...)))
  expect_length(call_with(), 2L)
  expect_error(call_with(pi0 = c(NA, 0.5)), "NA")
  expect_error(call_with(Q2 = c(1, NA)), "NA")
  expect_error(call_with(pi1 = c(0, 0.5)), "bound_cum_g")
  expect_error(call_with(pi0 = c(1.2, 0.5)), "bound_cum_g")
  expect_error(call_with(r0 = c(0L, 0L)), "r1")
  expect_error(call_with(Y = 1), "length")
})

test_that("bound_cum_g bounds below only, like ltmle gbounds = c(eps, 1)", {
  source(here::here("R", "meta_learner.R"))
  expect_equal(bound_cum_g(c(0, 0.0005, 0.5, 1), eps = 0.001), c(0.001, 0.001, 0.5, 1))
  expect_equal(bound_cum_g(0), CATE_EPS)
  # a near-zero cumulative g stays finite in the pseudo-outcome once bounded
  phi <- seq_dr_pseudo_outcome(1L, 1L, 40, 50, 45,
                               bound_cum_g(1e-8, 0.001), bound_cum_g(1e-12, 0.001))
  expect_equal(phi, 50 + (45 - 50) / 0.001 + (40 - 45) / 0.001)
})

test_that("stratum_estimates gives cell means with var(IC)/n, the HC0 sandwich", {
  source(here::here("R", "meta_learner.R"))
  set.seed(5)
  n <- 600
  s <- sample(c("111", "112", "221"), n, TRUE, prob = c(0.6, 0.3, 0.1))
  D <- cbind(never = rnorm(n, 50, 5), `now - never` = rnorm(n, -3, 8))
  out <- stratum_estimates(D, s)

  expect_named(out, c("strata_id", "estimand", "estimate", "variance"))
  expect_equal(out$strata_id, rep(c("111", "112", "221"), each = 2))
  expect_equal(out$estimand, rep(colnames(D), times = 3))
  purrr::walk(seq_len(nrow(out)), \(i) {
    j  <- s == out$strata_id[i]
    d  <- D[, out$estimand[i]]
    ic <- j * (d - mean(d[j])) / mean(j)
    expect_equal(out$estimate[i], mean(d[j]))
    expect_equal(out$variance[i], stats::var(ic) / n)
    hc0 <- diag(sandwich::vcovHC(lm(d ~ 0 + factor(s)), type = "HC0"))
    expect_equal(out$variance[i], unname(hc0[match(out$strata_id[i], sort(unique(s)))]) * n / (n - 1))
  })
  # one stratum holding everyone gives the whole-sample mean and var(D)/n
  all1 <- stratum_estimates(D, rep("overall", n))
  expect_equal(all1$estimate, unname(colMeans(D)))
  expect_equal(all1$variance, unname(purrr::map_dbl(asplit(D, 2), stats::var)) / n)
})

test_that("cum_g_diagnostics summarises followers and leaves an empty cell NA without warning", {
  source(here::here("R", "meta_learner.R"))
  sid     <- c("111", "111", "111", "112", "112")
  slab    <- c("F:W:H", "F:W:H", "F:W:H", "F:W:M", "F:W:M")
  follow  <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  pi1_raw <- c(0.0005, 0.2, 0.9, 0.3, 0.4)
  expect_no_warning(
    out <- cum_g_diagnostics(sid, slab, follow, bound_cum_g(pi1_raw, 0.001), pi1_raw, gbound = 0.001)
  )
  expect_named(out, c("strata_id", "strata_label", "n_j", "n_follow", "pct_follow",
                      "min_g", "median_g", "share_g_at_bound"))
  expect_equal(out$strata_id, c("111", "112"))
  expect_equal(out$strata_label, c("F:W:H", "F:W:M"))
  expect_equal(out$n_j, c(3, 2))
  expect_equal(out$n_follow, c(2, 0))
  expect_equal(out$pct_follow, c(100 * 2 / 3, 0))
  expect_equal(out$min_g, c(0.001, NA))
  expect_equal(out$median_g, c(stats::median(c(0.001, 0.2)), NA))
  expect_equal(out$share_g_at_bound, c(0.5, NA))
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

test_that("sl_honest_predict returns OOF on train rows and survives a learner that fails on the full fit", {
  source(here::here("R", "meta_learner.R"))

  set.seed(105)
  n  <- 300
  X  <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y  <- list(gaussian = 2 * X$x1 - X$x2 + rnorm(n),
             binomial = rbinom(n, 1, plogis(1.2 * X$x1 - 0.8 * X$x2)))
  fam   <- list(gaussian = stats::gaussian(), binomial = stats::binomial())
  train <- rep(c(TRUE, FALSE), c(200, 100))
  lib   <- c("SL.mean", "SL.glm")

  # train rows: out-of-fold; the rest: the full fit's ensemble prediction
  purrr::walk(names(y), \(f) {
    set.seed(11)
    got <- sl_honest_predict(y[[f]], X, train, fam[[f]], lib, toy_cv, f)
    set.seed(11)
    fit <- SuperLearner::SuperLearner(
      Y = y[[f]][train], X = X[train, , drop = FALSE], family = fam[[f]],
      SL.library = lib, cvControl = toy_cv,
      method = if (f == "binomial") "method.NNloglik" else "method.NNLS")
    expect_equal(got[train], sl_oof_predictions(fit, sum(train), f), tolerance = 1e-12)
    expect_equal(got[!train],
                 as.numeric(predict(fit, newdata = X[!train, , drop = FALSE], onlySL = TRUE)$pred),
                 tolerance = 1e-12)
  })

  # a learner that errors only on the full-size training set gets weight 0 and a
  # NULL fitLibrary entry; predict() without onlySL would then error. Assigned
  # into globalenv() as register_sl_wrappers() does.
  assign("SL.failfull", function(Y, X, newX, family, ...) {
    if (nrow(X) == sum(train)) stop("fails on the full fit")
    SuperLearner::SL.mean(Y, X, newX, family, ...)
  }, envir = globalenv())
  withr::defer(rm("SL.failfull", envir = globalenv()))
  got <- withCallingHandlers(
    sl_honest_predict(y$gaussian, X, train, stats::gaussian(), c("SL.glm", "SL.failfull"),
                      toy_cv, "failfull"),
    warning = \(w) if (grepl("Error in algorithm|already 0 for all failed", conditionMessage(w))) invokeRestart("muffleWarning"))
  expect_length(got, n)
  expect_true(all(is.finite(got)))
})

### toy generators for the two-time-point LTMLE tests (fit_ltmle_one,
### run_mice formulas, make_counterfactual_matrix, toy pipeline, estimate_cate)
# Synthetic frame shaped like build_data()'s MCS output: the 27 wide_cols_mcs
# columns with their types, plus the make_wide()/set_exposure() attributes.
# 12 balanced S cells. A0 ~ W; L1 ~ A0 + W (only income and employment respond
# to A0); A1 ~ A0 + L1; Y ~ W + L1 + A0 + A1. P(A0) ~ 0.28 and P(A1) ~ 0.23,
# so all 48 A0 x A1 x S cells are populated (min 33 at 400 per cell).
# attr(, "truth"): sample-average regime means E[Y(a0, a1)], with the
# counterfactual L1 drawn from the same noise.
# het scales stratum heterogeneity in the effects (-3 * a1 for men, +2 * a0 for
# low education). It draws nothing, so het = 0 is the toy above, unchanged.
# attr(, "truth_by_stratum"): the regime means and the contrasts with never per
# stratum, labelled as make_strata() labels them.
make_toy_2wave_df <- function(n_per_cell = 400, seed = 1, het = 0) {
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
      1 * emp1(a0) + 1 * male - 0.5 * low - 1.5 * a0 - 6 * a1 - 1 * a0 * a1 +
      het * (-3 * a1 * male + 2 * a0 * low)
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

  cell  <- interaction(d$sex_dv_fact_base, d$race_base, d$hiqual_dv_fact_base,
                       sep = ":", drop = TRUE)
  y_reg <- purrr::map(toy_regimes, \(r) mu_y(r[1], r[2]))
  y_dif <- purrr::map(y_reg[-1], \(y) y - y_reg[[1]])
  names(y_dif) <- paste(names(y_dif), "-", names(y_reg)[1])
  y_est <- c(y_reg, y_dif)

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
    truth             = purrr::map_dbl(toy_regimes, \(r) mean(mu_y(r[1], r[2]))),
    truth_by_stratum  = data.table::data.table(
      strata_label = rep(as.character(cell), length(y_est)),
      estimand     = rep(names(y_est), each = n),
      y            = unlist(y_est, use.names = FALSE)
    )[, .(truth = mean(y)), by = .(strata_label, estimand)]
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
  r <- fit_ltmle_one(mids, 1, "glm", "MCS", toy_regimes, toy_cv, toy_gbound)

  expect_s3_class(r, "data.table")
  expect_named(r, c("imp", "estimand", "type", "estimate", "variance"))
  expect_equal(r$estimand, ltmle_estimands)
  expect_equal(r$type, rep(c("mean", "ate"), c(4, 3)))
  expect_equal(r$imp, rep(1L, 7))

  # the ltmle data the contract specifies: W, A0, L1, A1, Y in node order, A integer
  A <- c("econ_dist_bin_0", "econ_dist_bin_1")
  x <- confounders("MCS")
  W <- grep("_1$", x, value = TRUE, invert = TRUE)
  L <- grep("_1$", x, value = TRUE)
  d <- as.data.frame(mice::complete(mids, 1))[, c(W, A[1], L, A[2], "sf12mcs_dv_1")]
  d[A] <- purrr::map(d[A], \(a) as.integer(as.character(a)))
  ltmle_sep <- function(abar) {
    suppressMessages(ltmle::ltmle(d, Anodes = A, Lnodes = L, Ynodes = "sf12mcs_dv_1",
                                  abar = abar, SL.library = "glm", gbounds = toy_gbound,
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
  expect_error(fit_ltmle_one(mids, 1:2, "glm", "MCS", toy_regimes, toy_cv, toy_gbound),
               "imp_idx")
  # ltmle checks only the names it is given, so SuperLearner would read an
  # unnamed list by position
  expect_error(fit_ltmle_one(mids, 1, "glm", "MCS", toy_regimes,
                             list(5L, FALSE, TRUE, NULL), toy_gbound),
               "sl_cv_control")
  # gbound goes straight to ltmle's gbounds, which needs c(lower, upper)
  expect_error(fit_ltmle_one(mids, 1, "glm", "MCS", toy_regimes, toy_cv, 0.001),
               "gbounds should have length 2")
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

test_that("pool_ltmle pools within extra by keys", {
  source(here::here("R", "pool_ltmle.R"))
  tab <- \(i, est) data.table::data.table(
    imp = i, strata_id = rep(c("111", "112"), each = 2),
    estimand = rep(c("never", "now - never"), 2), type = rep(c("mean", "ate"), 2),
    estimate = est, variance = c(0.04, 0.09, 0.05, 0.10))
  tabs <- list(tab(1L, c(50, -3, 48, -5)), tab(2L, c(51, -2, 47, -6)))
  p <- pool_ltmle(tabs, by = c("strata_id", "estimand", "type"))
  expect_named(p, c("strata_id", "estimand", "type", "effect", "se", "ll", "ul"))
  expect_equal(nrow(p), 4L)
  rr <- mice::pool.scalar(Q = c(48, 47), U = c(0.05, 0.05), n = Inf, k = 1)
  expect_equal(p[strata_id == "112" & estimand == "never", effect], rr$qbar)
  expect_equal(p[strata_id == "112" & estimand == "never", se], sqrt(rr$t))
})

test_that("pool_cate pools gate and overall with Rubin and keeps diag as a check", {
  source(here::here("R", "pool_ltmle.R"))
  source(here::here("R", "pool_cate.R"))
  gate_tab <- \(i, sids, est) data.table::data.table(
    imp = i, strata_id = rep(sids, each = 2), strata_label = rep(paste0("lab", sids), each = 2),
    estimand = rep(c("never", "now - never"), length(sids)), type = rep(c("mean", "ate"), length(sids)),
    estimate = est, variance = rep(c(0.04, 0.09), length(sids)))
  overall_tab <- \(i, est) data.table::data.table(
    imp = i, estimand = c("never", "now - never"), type = c("mean", "ate"),
    estimate = est, variance = c(0.01, 0.02))
  diag_tab <- \(i, sids, n_follow, min_g) data.table::data.table(
    imp = i, strata_id = sids, strata_label = paste0("lab", sids), regime = "now",
    n_j = 100L + i, n_follow = n_follow, pct_follow = n_follow, min_g = min_g,
    median_g = min_g, share_g_at_bound = 0)
  cate_one <- list(
    list(gate = gate_tab(1L, c("111", "112"), c(50, -3, 48, -5)), overall = overall_tab(1L, c(49, -4)),
         diag = diag_tab(1L, c("111", "112"), c(10L, 0L), c(0.2, NA))),
    # stratum 112 is absent from imputation 2: imputed sex/race/hiqual can empty a cell
    list(gate = gate_tab(2L, "111", c(51, -2)), overall = overall_tab(2L, c(50, -3)),
         diag = diag_tab(2L, "111", 12L, 0.3)))

  p <- expect_no_warning(pool_cate(cate_one))
  expect_named(p, c("gate", "overall", "diag", "m"))
  expect_equal(p$m, 2L)
  expect_equal(p$gate, pool_ltmle(purrr::map(cate_one, "gate"),
                                  by = c("strata_id", "strata_label", "estimand", "type")))
  expect_equal(p$overall, pool_ltmle(purrr::map(cate_one, "overall")))
  expect_equal(nrow(p$gate), 4L)
  # one imputation only: Rubin has no between-imputation variance, so the SE is NA
  expect_true(all(is.na(p$gate[strata_id == "112", se])))
  expect_named(p$diag, c("strata_id", "strata_label", "regime", "n_j_min", "n_j_max",
                         "n_follow_min", "pct_follow", "min_g", "median_g", "share_g_at_bound"))
  d111 <- p$diag[strata_id == "111"]
  expect_equal(c(d111$n_j_min, d111$n_j_max, d111$n_follow_min, d111$min_g), c(101, 102, 10, 0.2))
  expect_true(is.na(p$diag[strata_id == "112", min_g]))
  expect_error(pool_cate(cate_one[[1]]), "single")
})

test_that("toy pipeline: estimate_cate -> pool_cate on two imputations", {
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "strata_creation.R"))
  source(here::here("R", "meta_learner.R"))
  source(here::here("R", "pool_ltmle.R"))
  source(here::here("R", "pool_cate.R"))
  df   <- make_toy_2wave_df(n_per_cell = 400, het = 1)
  mids <- make_toy_2wave_mids(df, m = 2)
  p <- purrr::map(1:2, \(i) estimate_cate(mids, i, c("SL.mean", "SL.glm"), "MCS",
                                          toy_regimes, toy_cv, gbound = 0.001)) |> pool_cate()
  expect_equal(nrow(p$gate), 84L)
  expect_equal(p$overall$estimand, ltmle_estimands)
  expect_false(anyNA(p$gate$se))
  expect_true(all(p$gate$ll < p$gate$effect & p$gate$effect < p$gate$ul))
  tr <- merge(p$gate, attr(df, "truth_by_stratum"), by = c("strata_label", "estimand"))
  expect_lt(max(abs(tr$effect - tr$truth) / tr$se), 4)
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
                  c(confounders("MCS"), "econ_dist_bin_0", "econ_dist_bin_1"))
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
    fit_ltmle_one(mids, i, "glm", "MCS", toy_regimes, toy_cv, toy_gbound)
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

## ---- two-time-point DR-learner (toy data) ------------------------------------

test_that("make_toy_2wave_df keeps the pinned toy at het = 0 and adds stratum truths", {
  d0 <- make_toy_2wave_df(n_per_cell = 50)
  # values recorded from the generator before `het` existed
  expect_equal(attr(d0, "truth"),
               c(never = 51.4187601478786, before = 49.4570934812119,
                 now = 45.4187601478786, permanent = 42.4570934812119), tolerance = 1e-10)
  expect_equal(sum(d0$sf12mcs_dv_1), 29710.7993430944, tolerance = 1e-10)
  tb <- attr(make_toy_2wave_df(n_per_cell = 50, het = 1), "truth_by_stratum")
  expect_named(tb, c("strata_label", "estimand", "truth"))
  expect_equal(nrow(tb), 84L)
  expect_setequal(unique(tb$estimand), ltmle_estimands)
  expect_gt(diff(range(tb[estimand == "permanent - never", truth])), 3)

  # het adds exactly its term to each stratum's truth (same seed, same noise)
  t0 <- attr(make_toy_2wave_df(n_per_cell = 50, het = 0), "truth_by_stratum")
  m  <- merge(tb, t0, by = c("strata_label", "estimand"), suffixes = c("_1", "_0"))
  parts <- strsplit(m$strata_label, ":", fixed = TRUE)
  male  <- purrr::map_lgl(parts, \(p) p[1] == "Male")
  low   <- purrr::map_lgl(parts, \(p) p[3] == "Low")
  term  <- list(never = 0, before = 2 * low, now = -3 * male,
                permanent = -3 * male + 2 * low,
                "before - never" = 2 * low, "now - never" = -3 * male,
                "permanent - never" = -3 * male + 2 * low)
  expect_equal(nrow(m), 84L)
  expect_equal(m$truth_1 - m$truth_0,
               purrr::map_dbl(seq_len(nrow(m)), \(i) rep_len(term[[m$estimand[i]]], nrow(m))[i]),
               tolerance = 1e-10)
})

test_that("cate_nodes splits confounders like fit_ltmle_one and picks the outcome column", {
  source(here::here("R", "confounders.R"))
  source(here::here("R", "meta_learner.R"))
  purrr::walk(c("MCS", "PCS"), \(o) {
    nd <- cate_nodes(o)
    x  <- confounders(o)
    expect_identical(nd$W, grep("_1$", x, value = TRUE, invert = TRUE))
    expect_identical(nd$L, grep("_1$", x, value = TRUE))
    expect_identical(nd$A, c("econ_dist_bin_0", "econ_dist_bin_1"))
  })
  expect_identical(cate_nodes("MCS")$Y, "sf12mcs_dv_1")
  expect_identical(cate_nodes("PCS")$Y, "sf12pcs_dv_1")
  expect_error(cate_nodes("XXX"))
})

test_that("estimate_cate refuses stacked imputations and malformed arguments", {
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "strata_creation.R"))
  source(here::here("R", "meta_learner.R"))
  mids <- make_toy_2wave_mids(make_toy_2wave_df(n_per_cell = 100), m = 2)
  sl <- c("SL.mean", "SL.glm")
  expect_error(estimate_cate(mids, 1:2, sl, "MCS", toy_regimes, toy_cv), "imp_idx")
  expect_error(estimate_cate(mids, 1, sl, "MCS", toy_regimes, toy_cv, gbound = c(0.001, 1)), "gbound")
  expect_error(estimate_cate(mids, 1, sl, "MCS", toy_regimes, toy_cv, gbound = 0.5), "gbound")
  expect_error(estimate_cate(mids, 1, sl, "MCS", unname(toy_regimes), toy_cv), "regimes")
  expect_error(estimate_cate(mids, 1, sl, "MCS", list(never = c(0L, 2L)), toy_cv), "regimes")
  expect_error(estimate_cate(mids, 1, sl, "MCS", toy_regimes, list(5L, FALSE, TRUE)), "sl_cv_control")
})

test_that("estimate_cate recovers stratum means and differences on the two-wave toy", {
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "strata_creation.R"))
  source(here::here("R", "meta_learner.R"))
  # 800 per cell: at 400 a stratum's SE (~1 point) against the truths' spread
  # (sd ~1.8) puts cor(perm) below 0.8 for 1 seed in 5
  df   <- make_toy_2wave_df(n_per_cell = 800, het = 1)
  mids <- make_toy_2wave_mids(df, m = 1)
  res  <- estimate_cate(mids, 1, c("SL.mean", "SL.glm"), "MCS", toy_regimes, toy_cv, gbound = 0.001)

  expect_named(res, c("gate", "overall", "diag"))
  expect_named(res$gate, c("imp", "strata_id", "strata_label", "estimand", "type", "estimate", "variance"))
  expect_equal(nrow(res$gate), 84L)
  expect_type(res$gate$strata_id, "character")
  expect_named(res$overall, c("imp", "estimand", "type", "estimate", "variance"))
  expect_equal(res$overall$estimand, ltmle_estimands)
  expect_equal(res$overall$type, rep(c("mean", "ate"), c(4, 3)))
  expect_named(res$diag, c("imp", "strata_id", "strata_label", "regime", "n_j", "n_follow",
                           "pct_follow", "min_g", "median_g", "share_g_at_bound"))
  expect_equal(nrow(res$diag), 48L)
  expect_true(all(res$diag$min_g >= 0.001, na.rm = TRUE))

  # the whole-sample estimate is the n_j-weighted mean of the strata, exactly
  n_j <- unique(res$diag[, .(strata_id, n_j)])
  w   <- merge(res$gate, n_j, by = "strata_id")[, .(est = sum(estimate * n_j) / sum(n_j)), by = estimand]
  expect_equal(w$est[match(res$overall$estimand, w$estimand)], res$overall$estimate, tolerance = 1e-10)

  # recovery within 4 SE of the stratum truth; heterogeneity found; reductions negative
  tr <- merge(res$gate, attr(df, "truth_by_stratum"), by = c("strata_label", "estimand"))
  expect_equal(nrow(tr), 84L)
  expect_lt(max(abs(tr$estimate - tr$truth) / sqrt(tr$variance)), 4)
  perm <- tr[estimand == "permanent - never"]
  expect_gt(stats::cor(perm$estimate, perm$truth), 0.8)
  expect_true(all(perm$estimate < 0))

  # variance is calibrated: mean z^2 is ~1.1 over 20 toy seeds (0.85-1.66); a
  # variance off by 4x either way falls outside
  expect_gt(mean(((tr$estimate - tr$truth)^2) / tr$variance), 0.5)
  expect_lt(mean(((tr$estimate - tr$truth)^2) / tr$variance), 2)
})

test_that("estimate_cate overall agrees with fit_ltmle_one on the same imputation", {
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "strata_creation.R"))
  source(here::here("R", "meta_learner.R"))
  source(here::here("R", "fit_ltmle_one.R"))
  mids <- make_toy_2wave_mids(make_toy_2wave_df(n_per_cell = 400), m = 1)
  dr <- estimate_cate(mids, 1, c("SL.mean", "SL.glm"), "MCS", toy_regimes, toy_cv, gbound = 0.001)$overall
  lt <- fit_ltmle_one(mids, 1, "glm", "MCS", toy_regimes, toy_cv, toy_gbound)
  expect_equal(dr$estimand, lt$estimand)
  expect_equal(dr$type, lt$type)
  # same estimands, same sign: one-step DR and TMLE differ by well under one SE
  expect_lt(max(abs(dr$estimate - lt$estimate) / sqrt(lt$variance)), 1)
  # same influence-curve variance up to estimator differences (0.94-1.21 over 20 seeds)
  ratio <- dr$variance / lt$variance
  expect_gt(min(ratio), 0.5)
  expect_lt(max(ratio), 2)
})

test_that("estimate_cate is doubly robust: either nuisance alone recovers the truth", {
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "strata_creation.R"))
  source(here::here("R", "meta_learner.R"))
  # glm for one nuisance, the mean for the other: recovery rests on the first one's
  # wiring alone, so double robustness cannot mask a bug in it. SuperLearner finds
  # learners from globalenv(), not from the test environment.
  right_for <- \(binary) \(Y, X, newX, family, ...)
    (if (all(Y %in% 0:1) == binary) SuperLearner::SL.glm else SuperLearner::SL.mean)(Y, X, newX, family, ...)
  assign("SL.g_right", right_for(TRUE),  envir = globalenv())
  assign("SL.q_right", right_for(FALSE), envir = globalenv())
  withr::defer(rm("SL.g_right", "SL.q_right", envir = globalenv()))

  df    <- make_toy_2wave_df(n_per_cell = 800, het = 1)
  mids  <- make_toy_2wave_mids(df, m = 1)
  truth <- stats::setNames(c(attr(df, "truth"), attr(df, "truth")[-1] - attr(df, "truth")[1]),
                           ltmle_estimands)
  # max |z| cut-offs (overall, gate), set from a 20-seed sweep. Correct code: g right
  # <= 2.04 / 3.25, Q right <= 3.43 / 3.61. Mutants: pi1 not cumulative and g0 arm
  # swap (g right) >= 3.06 / 4.58; Q2 fitted on A0 = a0 rows (Q right) >= 10.1 / 5.04.
  cut <- list(SL.g_right = c(ov = 2.6, gate = 3.9), SL.q_right = c(ov = 6, gate = 4.3))
  purrr::walk(names(cut), \(lib) {
    r  <- estimate_cate(mids, 1, lib, "MCS", toy_regimes, toy_cv, gbound = 0.001)
    ov <- r$overall
    expect_lt(max(abs(ov$estimate - truth[ov$estimand]) / sqrt(ov$variance)), cut[[lib]][["ov"]])
    tr <- merge(r$gate, attr(df, "truth_by_stratum"), by = c("strata_label", "estimand"))
    expect_lt(max(abs(tr$estimate - tr$truth) / sqrt(tr$variance)), cut[[lib]][["gate"]])
  })
})

test_that("estimate_cate enforces gbound on the cumulative g", {
  source(here::here("R", "run_mice.R"))
  source(here::here("R", "sl_wrappers.R"))
  source(here::here("R", "confounders.R"))
  source(here::here("R", "strata_creation.R"))
  source(here::here("R", "meta_learner.R"))
  mids <- make_toy_2wave_mids(make_toy_2wave_df(n_per_cell = 400), m = 1)
  res  <- estimate_cate(mids, 1, c("SL.mean", "SL.glm"), "MCS", toy_regimes, toy_cv, gbound = 0.2)
  expect_true(all(res$diag$min_g >= 0.2, na.rm = TRUE))
  # P(permanent) is ~0.06 here, so the bound binds
  expect_true(any(res$diag$share_g_at_bound > 0, na.rm = TRUE))
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
