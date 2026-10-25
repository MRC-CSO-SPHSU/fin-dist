# Pure-mathematics tests for R/rma_reml.R. Data-free and fast by design:
# metafor is the ORACLE here, never a runtime dependency of the pipeline.
here::i_am("tests/test_rma_reml.R")
pacman::p_load(testthat, here)

# metafor's REML Fisher scoring stops at threshold = 1e-5 by default. Without
# tightening it, an oracle test measures metafor's convergence, not our maths.
ORACLE_CTRL <- list(threshold = 1e-12, maxiter = 10000, tol = 1e-14)

# A 2 x 2 x 3 design mirroring the real sex x race x education strata.
maihda_design <- function() {
  g <- expand.grid(sex  = c("Female", "Male"),
                   race = c("Non.white", "White"),
                   educ = c("High", "Medium", "Low"))
  model.matrix(~ sex + race + educ, g)      # 12 x 5
}
BETA_TRUE <- c(2, -0.3, 2.0, -0.2, -0.6)

test_that("rma_reml reproduces metafor's tau2 and Knapp-Hartung inference", {
  skip_if_not_installed("metafor")
  source(here::here("R", "rma_reml.R"))

  X <- maihda_design()
  set.seed(20260731)
  vi <- runif(12, 0.05, 5)
  yi <- as.vector(X %*% BETA_TRUE) + rnorm(12, 0, sqrt(1.2)) + rnorm(12, 0, sqrt(vi))

  hd <- rma_reml(yi, vi, X, knha = TRUE)
  mf <- metafor::rma(yi = yi, vi = vi, mods = X[, -1, drop = FALSE],
                     method = "REML", test = "knha", control = ORACLE_CTRL)

  expect_equal(hd$tau2,  mf$tau2,             tolerance = 1e-9)
  expect_equal(hd$beta,  as.vector(mf$beta),  tolerance = 1e-9)
  expect_equal(hd$se,    as.vector(mf$se),    tolerance = 1e-9)
  expect_equal(hd$tval,  as.vector(mf$zval),  tolerance = 1e-9)
  expect_equal(hd$pval,  as.vector(mf$pval),  tolerance = 1e-9)
  expect_equal(hd$ci.lb, as.vector(mf$ci.lb), tolerance = 1e-9)
  expect_equal(hd$ci.ub, as.vector(mf$ci.ub), tolerance = 1e-9)
  expect_equal(hd$dfs,   12L - 5L)
})

test_that("rma_reml returns exactly zero at the tau2 boundary", {
  skip_if_not_installed("metafor")
  source(here::here("R", "rma_reml.R"))

  # Cell estimates that scatter far LESS than their sampling SEs imply, so REML
  # has no between-cell variance left to attribute. The design doc flags this
  # branch as never exercised by the development probe -- it is constructed here
  # on purpose, not hoped for.
  X <- maihda_design()
  set.seed(11)
  vi <- rep(4, 12)
  yi <- as.vector(X %*% BETA_TRUE) + rnorm(12, 0, 0.2)

  hd <- rma_reml(yi, vi, X, knha = TRUE)
  mf <- metafor::rma(yi = yi, vi = vi, mods = X[, -1, drop = FALSE],
                     method = "REML", test = "knha", control = ORACLE_CTRL)

  expect_identical(hd$tau2, 0)              # exactly zero, not 1e-9
  expect_equal(mf$tau2, 0, tolerance = 1e-9)   # oracle agrees it is a boundary case
  expect_true(all(is.finite(hd$beta)))         # inference must survive tau2 = 0
  expect_true(all(is.finite(hd$se)))
  expect_equal(hd$beta, as.vector(mf$beta), tolerance = 1e-9)
  expect_equal(hd$se,   as.vector(mf$se),   tolerance = 1e-9)
})

test_that("rma_reml reproduces metafor's heterogeneity and moderator statistics", {
  skip_if_not_installed("metafor")
  source(here::here("R", "rma_reml.R"))

  X <- maihda_design()
  set.seed(20260731)
  vi <- runif(12, 0.05, 5)
  yi <- as.vector(X %*% BETA_TRUE) + rnorm(12, 0, sqrt(1.2)) + rnorm(12, 0, sqrt(vi))

  hd <- rma_reml(yi, vi, X, knha = TRUE)
  mf <- metafor::rma(yi = yi, vi = vi, mods = X[, -1, drop = FALSE],
                     method = "REML", test = "knha", control = ORACLE_CTRL)

  expect_equal(hd$se_tau2, mf$se.tau2, tolerance = 1e-9)
  expect_equal(hd$QE,      mf$QE,      tolerance = 1e-9)
  expect_equal(hd$QEp,     mf$QEp,     tolerance = 1e-9)
  expect_equal(hd$I2,      mf$I2,      tolerance = 1e-9)
  expect_equal(hd$H2,      mf$H2,      tolerance = 1e-9)
  expect_equal(hd$QM,      mf$QM,      tolerance = 1e-9)
  expect_equal(hd$QMp,     mf$QMp,     tolerance = 1e-9)
  expect_equal(hd$m_mods,  4L)

  # null model: no moderators to test
  X0  <- model.matrix(~ 1, data.frame(i = 1:12))
  hd0 <- rma_reml(yi, vi, X0, knha = TRUE)
  mf0 <- metafor::rma(yi = yi, vi = vi, method = "REML", test = "knha",
                      control = ORACLE_CTRL)
  expect_equal(hd0$tau2, mf0$tau2, tolerance = 1e-9)
  expect_equal(hd0$QE,   mf0$QE,   tolerance = 1e-9)
  expect_equal(hd0$I2,   mf0$I2,   tolerance = 1e-9)
  expect_true(is.na(hd0$QM))
  expect_equal(hd0$m_mods, 0L)
})

test_that("blup_reml reproduces metafor::blup", {
  skip_if_not_installed("metafor")
  source(here::here("R", "rma_reml.R"))

  X <- maihda_design()
  set.seed(20260731)
  vi <- runif(12, 0.05, 5)
  yi <- as.vector(X %*% BETA_TRUE) + rnorm(12, 0, sqrt(1.2)) + rnorm(12, 0, sqrt(vi))

  hd <- rma_reml(yi, vi, X, knha = TRUE)
  mf <- metafor::rma(yi = yi, vi = vi, mods = X[, -1, drop = FALSE],
                     method = "REML", test = "knha", control = ORACLE_CTRL)

  b_hd <- blup_reml(hd); b_mf <- metafor::blup(mf)
  expect_equal(b_hd$pred, as.vector(b_mf$pred), tolerance = 1e-9)
  expect_equal(b_hd$se,   as.vector(b_mf$se),   tolerance = 1e-9)

  # BLUPs shrink toward the additive surface: each sits between the raw
  # estimate and the fitted value, and is never more extreme than the raw one.
  fitted_add <- as.vector(X %*% hd$beta)
  expect_true(all(abs(b_hd$pred - fitted_add) <= abs(yi - fitted_add) + 1e-12))
})

test_that("rstudent_reml reproduces metafor, including on an outlier cell", {
  skip_if_not_installed("metafor")
  source(here::here("R", "rma_reml.R"))

  X <- maihda_design()
  set.seed(20260731)
  vi <- runif(12, 0.05, 5)
  yi <- as.vector(X %*% BETA_TRUE) + rnorm(12, 0, sqrt(1.2)) + rnorm(12, 0, sqrt(vi))

  # Inject an outlier, mimicking Male:Non.white:Low in the real data.
  yi[9] <- yi[9] - 8

  hd <- rstudent_reml(yi, vi, X, knha = TRUE)
  mf <- metafor::rma(yi = yi, vi = vi, mods = X[, -1, drop = FALSE],
                     method = "REML", test = "knha", control = ORACLE_CTRL)
  rs <- rstudent(mf)

  expect_equal(hd$resid, as.vector(rs$resid), tolerance = 1e-9)
  expect_equal(hd$se,    as.vector(rs$se),    tolerance = 1e-9)
  expect_equal(hd$z,     as.vector(rs$z),     tolerance = 1e-9)

  # Guard the specific mistake this formula invites: the textbook form
  # sqrt(v_j + tau2_del + vpred) omits the deleted fit's Knapp-Hartung factor
  # and is wrong by ~7% on the outlier cell. Assert we are NOT computing it.
  f9    <- rma_reml(yi[-9], vi[-9], X[-9, , drop = FALSE], knha = TRUE)
  X9    <- X[9, , drop = FALSE]
  naive <- sqrt(vi[9] + f9$tau2 + as.numeric(X9 %*% f9$vb %*% t(X9)))
  expect_gt(abs(naive - hd$se[9]) / hd$se[9], 0.03)
})

test_that("rma_reml refuses the unimplemented knha = FALSE variant", {
  source(here::here("R", "rma_reml.R"))
  X <- maihda_design()
  set.seed(20260731)
  vi <- runif(12, 0.05, 5)
  yi <- as.vector(X %*% BETA_TRUE) + rnorm(12, 0, sqrt(1.2)) + rnorm(12, 0, sqrt(vi))
  expect_error(rma_reml(yi, vi, X, knha = FALSE))
  expect_error(rstudent_reml(yi, vi, X, knha = FALSE))
})
