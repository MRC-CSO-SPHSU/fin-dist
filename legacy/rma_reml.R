# Random-effects meta-regression by REML, written from the mathematics.
# No package dependency: the pipeline must not carry metafor to the cluster.
# Verified against metafor 4.8.0 to ~1e-12 (design 2026-07-31, 4.4).
#
# Model:  y_j = x_j' beta + u_j + e_j,   u_j ~ N(0, tau2),  e_j ~ N(0, v_j)
#   =>    y   ~ N(X beta, diag(v) + tau2 I),   v KNOWN (squared GATE SEs).
#
# Given tau2 this is weighted least squares with w_j = 1/(v_j + tau2). tau2 is
# the ONLY quantity needing numerical work, and it is a scalar.

# Weighted-LS pieces at a given tau2.
.wls <- function(yi, vi, X, tau2) {
  wi   <- 1 / (vi + tau2)
  XtW  <- crossprod(X, diag(wi, nrow = length(yi)))
  XtWX <- XtW %*% X
  vb0  <- solve(XtWX)
  beta <- vb0 %*% (XtW %*% yi)
  r    <- as.vector(yi - X %*% beta)
  list(wi = wi, vb0 = vb0, beta = beta, r = r, RSS = sum(wi * r^2),
       logdet = as.numeric(determinant(XtWX, logarithm = TRUE)$modulus))
}

# REML profile log-likelihood, up to an additive constant free of tau2:
#   l(tau2) = -1/2 sum log(v_j + tau2) - 1/2 log|X'WX| - 1/2 sum w_j r_j^2
.ll_reml <- function(tau2, yi, vi, X) {
  f <- .wls(yi, vi, X, tau2)
  -0.5 * sum(log(vi + tau2)) - 0.5 * f$logdet - 0.5 * f$RSS
}

# REML score: dl/dtau2 = -1/2 tr(P) + 1/2 y'PPy, with the REML projection
#   P = W - WX (X'WX)^-1 X'W.   Its root is the interior REML estimate.
.score_reml <- function(tau2, yi, vi, X) {
  P  <- .P_matrix(yi, vi, X, tau2)
  0.5 * (sum((P %*% yi)^2) - sum(diag(P)))
}

# REML projection matrix at a given tau2 (shared by .score_reml and se(tau2)).
.P_matrix <- function(yi, vi, X, tau2) {
  W  <- diag(1 / (vi + tau2), nrow = length(yi))
  WX <- W %*% X
  W - WX %*% solve(crossprod(X, WX)) %*% t(WX)
}

rma_reml <- function(yi, vi, X, knha = TRUE, level = 0.95) {
  k <- length(yi); p <- ncol(X)
  stopifnot(k > p, all(vi > 0), nrow(X) == k, !anyNA(yi), !anyNA(vi))
  # Only the Knapp-Hartung formulation is implemented and verified. The z/chi-square
  # variant would need different reference distributions for pval, ci and QM, and is
  # unused, so it is refused rather than silently approximated.
  stopifnot(isTRUE(knha))

  ## tau2: bounded scalar maximisation, then polished on the score root.
  tau2_max <- 100 * max(vi)
  opt  <- optimize(.ll_reml, c(0, tau2_max), yi = yi, vi = vi, X = X,
                   maximum = TRUE, tol = .Machine$double.eps^0.5)
  tau2 <- opt$maximum
  if (.ll_reml(0, yi, vi, X) >= opt$objective) {
    tau2 <- 0                                   # REML boundary solution, exactly 0
  } else {
    # optimize()'s golden section only reaches sqrt(eps); the score root reaches ~1e-12.
    lo <- max(0, tau2 * 0.5); hi <- min(tau2_max, tau2 * 2 + 1e-8)
    s_lo <- .score_reml(lo, yi, vi, X); s_hi <- .score_reml(hi, yi, vi, X)
    if (is.finite(s_lo) && is.finite(s_hi) && s_lo * s_hi < 0)
      tau2 <- uniroot(.score_reml, c(lo, hi), yi = yi, vi = vi, X = X,
                      tol = .Machine$double.eps^0.75)$root
  }

  ## coefficients at the REML tau2
  f    <- .wls(yi, vi, X, tau2)
  dfs  <- k - p
  s2w  <- f$RSS / dfs                           # Knapp-Hartung scale factor
  vb   <- s2w * f$vb0                           # KH-scaled vb (knha asserted TRUE above)
  se   <- sqrt(diag(vb))
  beta <- as.vector(f$beta)
  crit <- qt(1 - (1 - level) / 2, dfs)

  # Strip names to match metafor output format
  beta <- unname(beta)
  se   <- unname(se)

  ## se(tau2) from the inverse REML Fisher information, I(tau2) = 1/2 tr(PP)
  P  <- .P_matrix(yi, vi, X, tau2)
  se_tau2 <- sqrt(2 / sum(P * t(P)))            # tr(PP) = sum(P * t(P)); P symmetric

  ## Heterogeneity, all on FIXED-effect weights w0 = 1/v (not the REML weights)
  P0  <- .P_matrix(yi, vi, X, 0)
  QE  <- as.numeric(crossprod(yi, P0 %*% yi))
  vt  <- dfs / sum(diag(P0))                    # typical within-cell variance
  I2  <- 100 * tau2 / (vt + tau2)
  H2  <- tau2 / vt + 1

  ## Omnibus moderator test (F, Knapp-Hartung), excluding the intercept
  QM <- QMp <- NA_real_; m_mods <- p - 1L
  if (m_mods > 0) {
    idx <- 2:p
    QM  <- as.numeric(crossprod(beta[idx],
                                solve(vb[idx, idx, drop = FALSE], beta[idx]))) / m_mods
    QMp <- pf(QM, m_mods, dfs, lower.tail = FALSE)
  }

  list(k = k, p = p, dfs = dfs, tau2 = tau2, s2w = s2w,
       beta = beta, se = se, tval = beta / se,
       pval = 2 * pt(-abs(beta / se), dfs),
       ci.lb = beta - crit * se, ci.ub = beta + crit * se,
       vb = vb, vb0 = f$vb0, yi = yi, vi = vi, X = X,
       se_tau2 = se_tau2, QE = QE, QEp = pchisq(QE, dfs, lower.tail = FALSE),
       vt = vt, I2 = I2, H2 = H2, QM = QM, QMp = QMp, m_mods = m_mods)
}

# BLUP: shrink each raw cell estimate toward the fitted additive surface.
#   lambda_j = tau2/(tau2 + v_j);  blup_j = lambda_j*y_j + (1-lambda_j)*x_j'beta
# Variance follows metafor::blup: lambda_j*v_j (NOT lambda_j^2*v_j) plus the
# uncertainty in the fitted surface.
blup_reml <- function(fit) {
  li    <- fit$tau2 / (fit$tau2 + fit$vi)
  Xb    <- as.vector(fit$X %*% fit$beta)
  pred  <- li * fit$yi + (1 - li) * Xb
  vpred <- li * fit$vi + (1 - li)^2 * rowSums((fit$X %*% fit$vb) * fit$X)
  data.frame(pred = pred, se = sqrt(vpred))
}

# Deleted (leave-one-out) residuals. Refit without cell j, predict it, standardise.
#
#   resid_j = y_j - x_j' beta_(-j)
#   se_j    = sqrt( s2w_(-j) * (v_j + tau2_(-j)) + x_j' vb_(-j) x_j )
#
# The deleted fit's Knapp-Hartung factor s2w_(-j) rescales BOTH v_j and
# tau2_(-j). This is metafor's influence.rma.uni and it is NOT the textbook
# form: when cell j is an outlier, dropping it makes the remaining cells fit
# better than their SEs predict, so s2w_(-j) < 1 and se_j can legitimately fall
# BELOW sqrt(v_j). Using sqrt(v_j + tau2_(-j) + vpred) instead disagrees by
# ~7-17% on exactly the outlier cell the decomposition hinges on. Do not
# "simplify" this.
rstudent_reml <- function(yi, vi, X, knha = TRUE) {
  # Only the Knapp-Hartung formulation is implemented and verified (see rma_reml()):
  # the deleted-fit se below combines fj$s2w with fj$vb, which is KH-scaled by
  # construction. A z/chi-square fj$vb would carry no such scale and this line
  # would be silently wrong. Refused rather than silently approximated.
  stopifnot(isTRUE(knha))
  k <- length(yi)
  as.data.frame(t(vapply(seq_len(k), function(j) {
    fj    <- rma_reml(yi[-j], vi[-j], X[-j, , drop = FALSE], knha = knha)
    Xj    <- X[j, , drop = FALSE]
    pred  <- as.numeric(Xj %*% fj$beta)
    vpred <- as.numeric(Xj %*% fj$vb %*% t(Xj))
    resid <- yi[j] - pred
    sr    <- sqrt(fj$s2w * (vi[j] + fj$tau2) + vpred)
    c(resid = resid, se = sr, z = resid / sr,
      tau2_del = fj$tau2, s2w_del = fj$s2w)
  }, numeric(5))))
}
