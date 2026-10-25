# MAIHDA decomposition of the DR-learner GATEs, one imputed dataset.
# Design: .claude/plans/2026-07-31-maihda-decomposition-design.md
#
# Consumes cate_one[[m]]$gate (12 rows) and returns the random-effects
# decomposition: tau2 with and without additive main effects, PCV, per-cell
# BLUPs and deleted residuals, and leave-one-cell-out PCV.

# Recover sex / race / education by SPLITTING THE LABEL, not by parsing the
# strata_id digits. strata_label comes from interaction(sex, race, hiqual,
# sep = ":"), and sanitize_factor_levels() in R/helpers.R (called from
# build_data()) has already stripped any ":" from the underlying factor
# levels via make.names(), so the split is exact. Reversing 100*sex + 10*race
# + hiqual would instead depend on as.numeric(factor) level ordering and, if
# that drifts, would silently mislabel cells rather than error.
.gate_factors <- function(gate_tbl) {
  parts <- strsplit(as.character(gate_tbl$strata_label), ":", fixed = TRUE)
  if (any(lengths(parts) != 3L)) {
    stop("fit_gate_meta: every strata_label must split into exactly three ",
         "parts on ':'; got lengths ",
         paste(sort(unique(lengths(parts))), collapse = "/"), ".")
  }
  m <- do.call(rbind, parts)
  out <- gate_tbl
  out$sex  <- factor(m[, 1])
  out$race <- factor(m[, 2])
  out$educ <- factor(m[, 3])

  n_cells <- nrow(unique(out[c("sex", "race", "educ")]))
  if (nrow(out) != 12L || n_cells != 12L) {
    stop("fit_gate_meta: expected a full 12-cell sex x race x education ",
         "design; got ", nrow(out), " rows spanning ", n_cells, " distinct cells.")
  }
  out
}

# PCV: share of between-stratum effect variance explained by ADDITIVE main
# effects. The remainder is the intersectional interaction.
.pcv <- function(tau2_null, tau2_main) (tau2_null - tau2_main) / tau2_null

fit_gate_meta <- function(gate_tbl) {
  imp_idx <- as.integer(gate_tbl$imp[1])
  g  <- .gate_factors(gate_tbl)
  yi <- g$estimate
  vi <- g$se^2
  X0 <- model.matrix(~ 1,                 data = g)
  X1 <- model.matrix(~ sex + race + educ, data = g)

  f0 <- rma_reml(yi, vi, X0, knha = TRUE)   # null: mu and tau2_null
  f1 <- rma_reml(yi, vi, X1, knha = TRUE)   # additive main effects

  scalars <- tibble::tibble(
    imp          = imp_idx,
    tau2_null    = f0$tau2, se_tau2_null = f0$se_tau2,
    tau2_main    = f1$tau2, se_tau2_main = f1$se_tau2,
    pcv          = .pcv(f0$tau2, f1$tau2),
    I2_null      = f0$I2,   I2_main      = f1$I2,
    vt_null      = f0$vt,   vt_main      = f1$vt,
    mu           = f0$beta[1], se_mu     = f0$se[1],
    QE_null      = f0$QE,   QE_main      = f1$QE,
    QM           = f1$QM,   QMp          = f1$QMp,
    s2w_main     = f1$s2w
  )

  coefs <- tibble::tibble(imp = imp_idx, term = colnames(X1),
                          estimate = f1$beta, se = f1$se)

  bl <- blup_reml(f1)
  rs <- rstudent_reml(yi, vi, X1, knha = TRUE)
  cells <- tibble::tibble(
    imp           = imp_idx,
    strata_id     = as.character(g$strata_id),
    strata_label  = as.character(g$strata_label),
    estimate      = yi,
    se            = g$se,
    additive_pred = as.vector(X1 %*% f1$beta),
    se_pred       = sqrt(rowSums((X1 %*% f1$vb) * X1)),
    blup          = bl$pred,   se_blup      = bl$se,
    del_resid     = rs$resid,  se_del_resid = rs$se,
    z             = rs$z,      tau2_del     = rs$tau2_del,
    n_j              = g$n_j,
    pct_exposed      = g$pct_exposed,
    min_g            = g$min_g,
    median_g         = g$median_g,
    share_g_at_bound = g$share_g_at_bound
  )

  # Leave-one-cell-out: how much of PCV rests on any single stratum? At J = 12
  # with a thin cell this is the difference between a finding and an artefact,
  # so it is a first-class output rather than a diagnostic.
  loco <- purrr::map_dfr(seq_len(nrow(g)), function(d) {
    d0 <- rma_reml(yi[-d], vi[-d], X0[-d, , drop = FALSE], knha = TRUE)
    d1 <- rma_reml(yi[-d], vi[-d], X1[-d, , drop = FALSE], knha = TRUE)
    tibble::tibble(
      imp = imp_idx,
      dropped_id    = as.character(g$strata_id[d]),
      dropped_label = as.character(g$strata_label[d]),
      tau2_null_d = d0$tau2, se_tau2_null_d = d0$se_tau2,
      tau2_main_d = d1$tau2, se_tau2_main_d = d1$se_tau2,
      pcv_d       = .pcv(d0$tau2, d1$tau2)
    )
  })

  list(scalars = scalars, coefs = coefs, cells = cells, loco = loco)
}
