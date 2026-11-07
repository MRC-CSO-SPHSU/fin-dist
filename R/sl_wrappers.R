# Custom SuperLearner learners for the ltmle and DR-learner fits.

# SuperLearner resolves SL.library names by string against globalenv(), and S3
# dispatch for the predict methods falls back to globalenv() too, so a fresh
# batchtools worker sees neither. Call before any SuperLearner()/ltmle() call.
# Keep the list in sync with `sl_libs` in _targets.R.
register_sl_wrappers <- function() {
  nms <- c("SL.xgboost.ltmle",
           "SL.rf.ltmle",   "predict.SL.rf.ltmle",
           "SL.poly.ltmle", "predict.SL.poly.ltmle", ".poly_design",
           "SL.poly2",     "SL.poly3",
           "SL.nnet.ltmle", "predict.SL.nnet.ltmle")
  src <- environment(register_sl_wrappers)
  purrr::walk(nms, \(nm) assign(nm, get(nm, envir = src), envir = globalenv()))
  invisible(nms)
}

# XGBoost: non-linearities and interactions, even when n > p.
# With xgboost > 3.0 (3.2.1.1 installed), stock SuperLearner::SL.xgboost fits
# y = as.factor(Y) under family = binomial(). ltmle() passes binomial() for
# every Q regression too, so the continuous [0,1] Q outcomes error ("Got
# non-binary factor 'y'"). Non-binary Y is therefore fit as gaussian.
SL.xgboost.ltmle <- function(Y, X, newX, family, ...) {
  Y <- as.numeric(Y)
  if (!all(Y %in% c(0, 1))) family <- gaussian()
  SuperLearner::SL.xgboost(Y = Y, X = X, newX = newX, family = family, ...)
}

# No glmnet learner: ridge and lasso only help when p >> n, not the case here.

# Random forest: 500 trees, at least 20 observations per terminal node.
# Always a regression forest, whatever `family` says: leaf values are averages
# of training outcomes, so predictions cannot leave the range of Y, in [0,1] for
# the Q-model's Ystar and the {0,1} g-model alike. No family branching needed.
# num.threads = 1: targets supplies the parallelism (4 future.callr workers
# locally, one SLURM job per branch on the cluster), so ranger claiming every
# core would oversubscribe the node.
# `id` is absorbed, not forwarded: SuperLearner passes it to every learner, and
# ranger would warn "Unused arguments: id" on each of thousands of fits.
SL.rf.ltmle <- function(Y, X, newX, family,
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
  class(out) <- "SL.rf.ltmle"
  list(pred = pred, fit = out)
}

predict.SL.rf.ltmle <- function(object, newdata, ...) {
  NXm <- data.matrix(newdata)
  colnames(NXm) <- object$xnames
  stats::predict(object$object, data = NXm)$predictions
}

# Polynomial main-effects regression of a given degree.
# Only genuinely continuous columns are expanded: any power of a 0/1 dummy is
# the dummy itself, so expanding it would make the design singular. Columns are
# standardised before the powers, because raw x^3 on the SF-12 or age scale is
# badly conditioned.
# quasibinomial for a [0,1] outcome (ltmle's Q outcomes, i.e. the rescaled Y or
# the previous step's Q, and the binary g-model): the logit link keeps fitted
# values inside (0,1), and the quasi- family avoids glm's non-integer-successes
# warning. Raw MCS/PCS, as fit by estimate_cate(), falls through to gaussian.
SL.poly.ltmle <- function(Y, X, newX, family,
                         obsWeights = rep(1, length(Y)),
                         degree = 2, ...) {
  Xm  <- data.matrix(X)
  NXm <- data.matrix(newX)
  colnames(Xm) <- colnames(NXm) <- make.names(colnames(Xm), unique = TRUE)

  cont <- purrr::map_lgl(asplit(Xm, 2), \(x) length(unique(x)) > 2)
  ctr  <- colMeans(Xm[, cont, drop = FALSE])
  scl  <- purrr::map_dbl(asplit(Xm[, cont, drop = FALSE], 2), stats::sd)
  scl[scl < 1e-8] <- 1

  Xd  <- .poly_design(Xm,  cont, ctr, scl, degree)
  NXd <- .poly_design(NXm, cont, ctr, scl, degree)

  fam <- if (all(Y >= 0 & Y <= 1)) stats::quasibinomial() else stats::gaussian()
  fit <- stats::glm(Y ~ ., data = data.frame(Y = Y, Xd),
                    family = fam, weights = obsWeights)

  pred <- as.numeric(stats::predict(fit, newdata = NXd, type = "response"))
  out  <- list(object = fit, cont = cont, ctr = ctr, scl = scl, degree = degree)
  class(out) <- "SL.poly.ltmle"
  list(pred = pred, fit = out)
}

.poly_design <- function(M, cont, ctr, scl, degree) {
  Cm <- scale(M[, cont, drop = FALSE], center = ctr, scale = scl)
  pw <- purrr::map(seq(2, degree), \(d) {
    P <- Cm^d
    colnames(P) <- paste0(colnames(Cm), "_p", d)
    P
  }) |>
    purrr::reduce(cbind)
  as.data.frame(cbind(M[, !cont, drop = FALSE], Cm, pw))
}

predict.SL.poly.ltmle <- function(object, newdata, ...) {
  NXm <- data.matrix(newdata)
  colnames(NXm) <- make.names(colnames(NXm), unique = TRUE)
  NXd <- .poly_design(NXm, object$cont, object$ctr, object$scl, object$degree)
  as.numeric(stats::predict(object$object, newdata = NXd, type = "response"))
}

SL.poly2 <- function(...) SL.poly.ltmle(..., degree = 2)
SL.poly3 <- function(...) SL.poly.ltmle(..., degree = 3)

# Neural network: complex non-linearities.
#
# Branches on range(Y), NOT on `family`, which misdescribes the outcome. ltmle()
# rescales Y to [0,1] and calls SuperLearner with family = binomial() for the g
# and Q regressions alike (Estimate() turns quasibinomial into binomial first),
# yet each Q outcome (the rescaled Y, or the previous step's Q) is continuous in
# [0,1]. tmle() passed family = "gaussian" for the same [0,1] outcome. ltmle
# clips ensemble predictions outside [0,1] back into [0,1] with a warning, and
# the clipping shifts the marginal mean, so the bounded branch keeps a logistic
# output unit.
# Bounded branch: fit_ltmle_one(), Q on Ystar in [0,1] and g on A in {0,1}.
# Unbounded branch: estimate_cate(), mu0/mu1 on the raw MCS/PCS scale.
#
# nnet applies one global `decay`/`rang` to every weight, so unstandardised
# inputs leave the optimisation ill-conditioned. That, not capacity, is why
# stock SL.nnet oscillates between coef 0 and a runaway fit across imputations.
# Fix: standardise the design matrix, keep the best of several random restarts.
#
# One hidden layer of 50 units. nnet implements a single hidden layer only, so a
# 50/50/50/50 stack is not expressible; this is the widest single-layer
# equivalent (27 inputs, 50 hidden units, 1 output, ~1,450 weights).
#
# maxit = 200 and n.restarts = 2 come from measurements on wave-3 PCS
# (n = 21,233, p = 27): one fit costs 27s at maxit = 200 against 107s at
# maxit = 1000. Per imputation, ltmleMSM() runs 7 SuperLearner fits, each V-fold
# (V = 5) plus a full refit: g at A0 and at A1 once each, Q at the Y node once,
# and Q at the L1 block once per regime (4).
#
# The cut is safe for the marginal estimands only. nnet reports convergence = 1
# (iteration cap hit) at maxit = 1000 as well as at 200, and the two disagree by
# 1.7 PCS points on average per individual (max 14.4), yet the marginal mean
# moves by only 0.011 points. Individual predictions feed estimate_cate()'s
# mu0/mu1, so raise maxit if the DR-learner arm looks unstable.
#
# n.restarts is low because standardisation plus decay already removed the
# instability the restarts insured against: across 6 free random inits the
# marginal mean has sd 0.011 PCS points, against the 0.287 that stock SL.nnet
# contributed to the between-imputation variance of EY0.
SL.nnet.ltmle <- function(Y, X, newX, family,
                         obsWeights = rep(1, length(Y)), id = NULL,
                         size = 50, decay = 0.02, maxit = 200,
                         n.restarts = 2, ...) {
  if (!requireNamespace("nnet", quietly = TRUE)) stop("package 'nnet' required")

  Xm  <- data.matrix(X)
  ctr <- colMeans(Xm)
  scl <- purrr::map_dbl(asplit(Xm, 2), stats::sd); scl[scl < 1e-8] <- 1
  Xs  <- scale(Xm, center = ctr, scale = scl)
  NXs <- scale(data.matrix(newX), center = ctr, scale = scl)

  bounded <- all(Y >= 0 & Y <= 1)
  # Continuous Y needs scaling too: `decay` penalises the output weights, so an
  # outcome on the ~50-point SF-12 scale would be shrunk into a near-flat fit.
  y.ctr <- if (bounded) 0 else mean(Y)
  y.scl <- if (bounded) 1 else stats::sd(Y)
  if (y.scl < 1e-8) y.scl <- 1
  Yt <- (Y - y.ctr) / y.scl

  fits <- purrr::map(seq_len(n.restarts), \(i) try(nnet::nnet(
    x = Xs, y = Yt, weights = obsWeights,
    size = size, decay = decay, maxit = maxit,
    linout = !bounded, entropy = bounded, softmax = FALSE,
    rang = 0.7, trace = FALSE,
    MaxNWts = size * (ncol(Xs) + 2) + 100
  ), silent = TRUE)) |>
    purrr::discard(\(f) inherits(f, "try-error"))
  if (!length(fits)) stop("SL.nnet.ltmle: all nnet restarts failed")
  # Same data and decay across restarts, so $value is directly comparable.
  # Strict `<` keeps the first best fit on ties.
  best <- purrr::reduce(fits, \(b, f) if (f$value < b$value) f else b)

  pred <- as.numeric(stats::predict(best, newdata = NXs, type = "raw")) * y.scl + y.ctr
  out  <- list(object = best, ctr = ctr, scl = scl, y.ctr = y.ctr, y.scl = y.scl)
  class(out) <- "SL.nnet.ltmle"
  list(pred = pred, fit = out)
}

predict.SL.nnet.ltmle <- function(object, newdata, ...) {
  NXs <- scale(data.matrix(newdata), center = object$ctr, scale = object$scl)
  as.numeric(stats::predict(object$object, newdata = NXs, type = "raw")) *
    object$y.scl + object$y.ctr
}
