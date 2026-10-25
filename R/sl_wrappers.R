# Custom SuperLearner learners for the TMLE step.

# SuperLearner resolves SL.library names by string against globalenv(), and S3
# dispatch for the predict methods falls back to globalenv() as well, so on a
# fresh batchtools worker neither would be visible. Call this before any
# SuperLearner()/tmle() call. Keep the list in sync with `sl_libs` in _targets.R.
register_sl_wrappers <- function() {
  nms <- c("SL.xgboost.tmle",
           "SL.rf.tmle",   "predict.SL.rf.tmle",
           "SL.poly.tmle", "predict.SL.poly.tmle", ".poly_design",
           "SL.poly2",     "SL.poly3",
           "SL.nnet.tmle", "predict.SL.nnet.tmle")
  src <- environment(register_sl_wrappers)
  for (nm in nms) assign(nm, get(nm, envir = src), envir = globalenv())
  invisible(nms)
}

## XGboost covers non-linearities and interactions, even in contexts with n>p
SL.xgboost.tmle <- function(Y, X, newX, family, ...) {
  Y <- as.numeric(Y)
  if (!all(Y %in% c(0, 1))) family <- gaussian()
  SuperLearner::SL.xgboost(Y = Y, X = X, newX = newX, family = family, ...)
}

# Ridge and lasso only useful when p>>n, which is not the case here.
# SL.glmnet.tmle <- function(Y, X, newX, family, ...) {
#  if (!all(Y %in% c(0, 1))) family <- gaussian()
#  SuperLearner::SL.glmnet(Y = Y, X = X, newX = newX, family = family, ...)
#}

# Random forest: 500 trees, minimum 20 observations per terminal node.
#
# Always grown as a regression forest, whatever `family` says. Leaf values are
# averages of training outcomes, so predictions stay inside the convex hull of Y
# by construction -- in [0,1] for the Q-model's Ystar and for the {0,1} g-model
# alike. That avoids the out-of-range predictions an unbounded learner produces
# on the rescaled outcome, so no family branching is needed here.
#
# num.threads = 1 on purpose: targets already runs 8 workers (or one SLURM job
# per branch), so letting ranger claim every core would oversubscribe the node.
# `id` is absorbed rather than forwarded: SuperLearner passes it to every
# learner, and ranger would warn "Unused arguments: id" on each of the thousands
# of fits in the branch grid.
SL.rf.tmle <- function(Y, X, newX, family,
                       obsWeights = rep(1, length(Y)), id = NULL,
                       num.trees = 500, min.node.size = 20, ...) {
  if (!requireNamespace("ranger", quietly = TRUE)) stop("package 'ranger' required")

  Xm  <- data.matrix(X)
  NXm <- data.matrix(newX)
  colnames(Xm) <- colnames(NXm) <- make.names(colnames(Xm), unique = TRUE)

  fit <- ranger::ranger(
    x = Xm, y = Y,
    num.trees     = num.trees,
    min.node.size = min.node.size,
    case.weights  = obsWeights,
    num.threads   = 1,
    ...
  )
  pred <- stats::predict(fit, data = NXm)$predictions
  out  <- list(object = fit, xnames = colnames(Xm))
  class(out) <- "SL.rf.tmle"
  list(pred = pred, fit = out)
}

predict.SL.rf.tmle <- function(object, newdata, ...) {
  NXm <- data.matrix(newdata)
  colnames(NXm) <- object$xnames
  stats::predict(object$object, data = NXm)$predictions
}

# Polynomial main-effects regression of a given degree.
#
# Only genuinely continuous columns are expanded: any power of a 0/1 dummy is
# the dummy itself, so expanding those would make the design singular. Columns
# are standardised before the powers are taken, because raw x^3 on the SF-12 or
# age scale is badly conditioned.
#
# quasibinomial for a [0,1] outcome (Ystar under fluctuation = "logistic", and
# the binary g-model): the logit link keeps fitted values inside (0,1), and the
# quasi- family avoids glm's non-integer-successes warning. Raw MCS/PCS, as fit
# by estimate_cate(), falls through to gaussian.
SL.poly.tmle <- function(Y, X, newX, family,
                         obsWeights = rep(1, length(Y)),
                         degree = 2, ...) {
  Xm  <- data.matrix(X)
  NXm <- data.matrix(newX)
  colnames(Xm) <- colnames(NXm) <- make.names(colnames(Xm), unique = TRUE)

  cont <- apply(Xm, 2, function(x) length(unique(x)) > 2)
  ctr  <- colMeans(Xm[, cont, drop = FALSE])
  scl  <- apply(Xm[, cont, drop = FALSE], 2, stats::sd); scl[scl < 1e-8] <- 1

  Xd  <- .poly_design(Xm,  cont, ctr, scl, degree)
  NXd <- .poly_design(NXm, cont, ctr, scl, degree)

  fam <- if (all(Y >= 0 & Y <= 1)) stats::quasibinomial() else stats::gaussian()
  fit <- stats::glm(Y ~ ., data = data.frame(Y = Y, Xd),
                    family = fam, weights = obsWeights)

  pred <- as.numeric(stats::predict(fit, newdata = NXd, type = "response"))
  out  <- list(object = fit, cont = cont, ctr = ctr, scl = scl, degree = degree)
  class(out) <- "SL.poly.tmle"
  list(pred = pred, fit = out)
}

.poly_design <- function(M, cont, ctr, scl, degree) {
  Cm <- scale(M[, cont, drop = FALSE], center = ctr, scale = scl)
  pw <- do.call(cbind, lapply(seq(2, degree), function(d) {
    P <- Cm^d
    colnames(P) <- paste0(colnames(Cm), "_p", d)
    P
  }))
  as.data.frame(cbind(M[, !cont, drop = FALSE], Cm, pw))
}

predict.SL.poly.tmle <- function(object, newdata, ...) {
  NXm <- data.matrix(newdata)
  colnames(NXm) <- make.names(colnames(NXm), unique = TRUE)
  NXd <- .poly_design(NXm, object$cont, object$ctr, object$scl, object$degree)
  as.numeric(stats::predict(object$object, newdata = NXd, type = "response"))
}

SL.poly2 <- function(...) SL.poly.tmle(..., degree = 2)
SL.poly3 <- function(...) SL.poly.tmle(..., degree = 3)

# Neural networks are helpful for complex non-linearities.
#
# Branches on range(Y), NOT on `family`. tmle() rescales the outcome to [0,1]
# in .initStage1() when fluctuation = "logistic", but still forwards
# family = "gaussian" to estimateQ(), so `family` misdescribes the outcome the
# Q-model actually sees. Stock SL.nnet trusts it, takes the gaussian branch and
# fits an unbounded linear output unit to a [0,1] outcome; predictions escape
# the unit interval, get clipped by .bound() before qlogis(), and the clipping
# shifts the marginal mean. Call sites hitting each branch:
#   bounded   -- fit_tmle_one(): Q on Ystar in [0,1], g on A in {0,1}
#   unbounded -- estimate_cate(): mu0/mu1 on the raw MCS/PCS scale
#
# nnet applies one global `decay`/`rang` to every weight, so unstandardised
# inputs leave the optimisation ill-conditioned. That, not capacity, is why
# stock SL.nnet oscillates between coef 0 and a runaway fit across imputations:
# standardise the design matrix and keep the best of several random restarts.
#
# One hidden layer of 50 units. `nnet` implements a single hidden layer only, so
# a 50/50/50/50 stack is not expressible here; this is the widest single-layer
# equivalent (27 -> 50 -> 1, ~1,450 weights).
#
# maxit = 200 and n.restarts = 2 are set from measurements on wave-3 PCS
# (n = 21,233, p = 27), where one fit costs 27s at maxit = 200 against 107s at
# maxit = 1000. cvQinit = TRUE makes tmle() fit the Q-model 36 times per call
# (one full SuperLearner plus a 5-fold outer loop each running its own 5-fold
# SuperLearner) plus 6 for the g-model, so per-fit cost is multiplied by ~42.
#
# The cut is safe for the marginal estimands and only those: nnet reports
# convergence = 1 (iteration cap hit) at maxit = 1000 as well as at 200, and the
# two disagree by 1.7 PCS points on average per individual (max 14.4) -- yet the
# marginal mean moves by only 0.011 points. Individual predictions feed
# estimate_cate()'s mu0/mu1, so raise maxit if the DR-learner arm looks unstable.
#
# n.restarts is low because standardisation plus decay already removed the
# instability the restarts were insurance against: across 6 free random inits the
# marginal mean has sd 0.011 PCS points, against the 0.287 that stock SL.nnet
# used to contribute to the between-imputation variance of EY0.
SL.nnet.tmle <- function(Y, X, newX, family,
                         obsWeights = rep(1, length(Y)), id = NULL,
                         size = 50, decay = 0.02, maxit = 200,
                         n.restarts = 2, ...) {
  if (!requireNamespace("nnet", quietly = TRUE)) stop("package 'nnet' required")

  Xm  <- data.matrix(X)
  ctr <- colMeans(Xm)
  scl <- apply(Xm, 2, stats::sd); scl[scl < 1e-8] <- 1
  Xs  <- scale(Xm, center = ctr, scale = scl)
  NXs <- scale(data.matrix(newX), center = ctr, scale = scl)

  bounded <- all(Y >= 0 & Y <= 1)
  # Continuous Y needs scaling too: `decay` penalises the output weights, so an
  # outcome on the ~50-point SF-12 scale would be shrunk into a near-flat fit.
  y.ctr <- if (bounded) 0 else mean(Y)
  y.scl <- if (bounded) 1 else stats::sd(Y)
  if (y.scl < 1e-8) y.scl <- 1
  Yt <- (Y - y.ctr) / y.scl

  best <- NULL
  for (i in seq_len(n.restarts)) {
    fit <- try(nnet::nnet(
      x = Xs, y = Yt, weights = obsWeights,
      size = size, decay = decay, maxit = maxit,
      linout = !bounded, entropy = bounded, softmax = FALSE,
      rang = 0.7, trace = FALSE,
      MaxNWts = size * (ncol(Xs) + 2) + 100
    ), silent = TRUE)
    if (inherits(fit, "try-error")) next
    # Same data and decay across restarts, so $value is directly comparable.
    if (is.null(best) || fit$value < best$value) best <- fit
  }
  if (is.null(best)) stop("SL.nnet.tmle: all nnet restarts failed")

  pred <- as.numeric(stats::predict(best, newdata = NXs, type = "raw")) * y.scl + y.ctr
  out  <- list(object = best, ctr = ctr, scl = scl, y.ctr = y.ctr, y.scl = y.scl)
  class(out) <- "SL.nnet.tmle"
  list(pred = pred, fit = out)
}

predict.SL.nnet.tmle <- function(object, newdata, ...) {
  NXs <- scale(data.matrix(newdata), center = object$ctr, scale = object$scl)
  as.numeric(stats::predict(object$object, newdata = NXs, type = "raw")) *
    object$y.scl + object$y.ctr
}