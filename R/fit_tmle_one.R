# Single-point TMLE (via the tmle package) for one imputed database.
fit_tmle_one <- function(wide_mids, imp_idx, sl_libs, outcome, gbound = NULL) {
  # One imputation per call. mice::complete() accepts a vector for `action` and
  # silently returns the imputations *stacked*, so an unbranched target would
  # fit a single TMLE on m * n rows instead of erroring. Fail fast instead.
  stopifnot(
    "imp_idx must be a single imputation index -- is the target missing pattern = map(tmle_imp_idx)?" =
      length(imp_idx) == 1L
  )
  # Assign custom SL wrappers to globalenv() so SuperLearner resolves them on each batchtools worker.
  register_sl_wrappers()

  # Confounders (W), all measured pre-exposure. Shared with estimate_cate()
  # via R/confounders.R -- do not inline a copy here.
  W <- confounders(outcome)

  dat <- mice::complete(wide_mids, action = imp_idx)
  A    <- as.integer(as.character(dat$econ_dist_bin_0))
  Y    <- dat[[if (outcome == "MCS") "sf12mcs_dv_0" else "sf12pcs_dv_0"]]
  Wmat <- model.matrix(~ ., dat[W])[, -1]

  tmle::tmle(
    Y            = Y,
    A            = A,
    W            = Wmat,
    Q.SL.library = sl_libs,
    g.SL.library = sl_libs,
    cvQinit      = TRUE,
    V.Q          = 5,
    V.g          = 5,
    gbound       = gbound, 
    family       = "gaussian", # for continuous outcomes
    fluctuation = "logistic" # to bound the outcome to [0,1] for the logistic fluctuation step
  )
}
