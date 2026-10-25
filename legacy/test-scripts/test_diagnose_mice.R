# unit tests for tests/diagnose_mice.R on small synthetic data
here::i_am("tests/test_diagnose_mice.R")

pacman::p_load(testthat, mice, here)

source(here::here("tests", "diagnose_mice.R"))

## helpers ------------------------------------------------------------------
# data engineered so mice logs events at initialization:
#   - `dup` duplicates `a`  -> "collinear" event for one member of the pair
#                              (which one is a mice implementation detail)
#   - `const` is constant   -> "constant"  event, dep == "" (setup)
make_messy <- function() {
  set.seed(42)
  n  <- 60
  df <- data.frame(
    a = rnorm(n),
    b = rnorm(n),
    f = factor(sample(c("lo", "hi"), n, replace = TRUE))
  )
  df$dup   <- df$a   # duplicated BEFORE adding NAs, so dup stays complete
  df$const <- 1
  df$a[sample(n, 12)] <- NA
  df$b[sample(n, 12)] <- NA
  df$f[sample(n, 12)] <- NA
  list(data = df,
       mids = mice::mice(df, m = 2, maxit = 2, seed = 1, printFlag = FALSE))
}

make_clean <- function() {
  set.seed(42)
  n  <- 60
  df <- data.frame(
    a = rnorm(n),
    b = rnorm(n),
    f = factor(sample(c("lo", "hi"), n, replace = TRUE))
  )
  df$a[sample(n, 12)] <- NA
  df$b[sample(n, 12)] <- NA
  df$f[sample(n, 12)] <- NA
  list(data = df,
       mids = mice::mice(df, m = 2, maxit = 2, seed = 1, printFlag = FALSE))
}

## tests --------------------------------------------------------------------
test_that("diagnose_mice breaks down logged events and flags culprits", {
  l <- make_messy()
  expect_false(is.null(l$mids$loggedEvents))  # precondition of this test

  res <- NULL
  expect_output(res <- diagnose_mice(l$mids, l$data), "logged events")
  expect_output(diagnose_mice(l$mids, l$data), "setup")

  expect_type(res, "list")
  expect_named(res, c("events", "summary"))
  expect_s3_class(res$summary, "data.frame")
  expect_setequal(
    names(res$summary),
    c("variable", "n_missing", "pct_missing", "method",
      "n_logged_as_dep", "n_logged_as_out")
  )
  expect_setequal(res$summary$variable, names(l$data))

  # the engineered culprits are flagged as removed-from-other-models:
  # mice drops ONE member of the collinear pair (a, dup) - which member is
  # a mice implementation detail - plus the constant column
  s <- res$summary
  expect_gte(s$n_logged_as_out[s$variable == "a"] +
             s$n_logged_as_out[s$variable == "dup"], 1L)
  expect_gte(s$n_logged_as_out[s$variable == "const"], 1L)

  # missingness is computed from the input data
  expect_equal(s$n_missing[s$variable == "a"], 12L)
  expect_equal(s$n_missing[s$variable == "const"], 0L)
})

test_that("diagnose_mice handles a clean mids", {
  l <- make_clean()
  expect_null(l$mids$loggedEvents)  # precondition of this test

  res <- NULL
  expect_output(res <- diagnose_mice(l$mids, l$data), "No logged events")
  expect_equal(sum(res$summary$n_logged_as_dep), 0L)
  expect_equal(sum(res$summary$n_logged_as_out), 0L)
})

test_that("diagnose_mice writes plots when plot_dir is given, none otherwise", {
  l  <- make_clean()
  td <- file.path(tempdir(), "mice_diag_test")
  unlink(td, recursive = TRUE)

  expect_output(diagnose_mice(l$mids, l$data, plot_dir = td), "Plots written")

  files <- list.files(td)
  expect_true("trace.png" %in% files)
  expect_true("density_a.png" %in% files)   # continuous -> densityplot
  expect_true("density_b.png" %in% files)
  expect_true("props_f.png" %in% files)     # factor -> proportion barchart

  # plot_dir = NULL must not create the directory
  unlink(td, recursive = TRUE)
  diagnose_mice(l$mids, l$data, plot_dir = NULL) |> capture.output() -> ignored
  expect_false(dir.exists(td))
})
