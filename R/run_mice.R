# Congeniality guard for the interaction terms run_mice() adds. If mice drops one
# for collinearity, it records the fact in loggedEvents and carries on silently:
# the very attenuation those terms exist to prevent.
assert_congenial_terms <- function(mids, a_col, y_col) {
  keys <- c(a_col, y_col)  # before any early return, so a missing argument errors

  le <- mids$loggedEvents
  if (!is.data.frame(le) || nrow(le) == 0L || !"out" %in% names(le)) return(invisible(TRUE))

  # Fixed substring match: factor dummies carry the level ("econ_dist_bin_01"), and
  # "econ_dist_bin_lagged_0" does not contain "econ_dist_bin_0".
  dropped <- purrr::map(strsplit(le$out, ",\\s*"), \(tm) {
    purrr::keep(trimws(tm), \(term)
      grepl(":", term, fixed = TRUE) && purrr::some(keys, \(k) grepl(k, term, fixed = TRUE)))
  })

  wipeout <- grepl("All predictors are constant", le$out, fixed = TRUE)

  bad <- lengths(dropped) > 0L | wipeout
  if (!any(bad)) return(invisible(TRUE))

  stop("run_mice: congeniality broken -- mice dropped ", sum(lengths(dropped)),
       " exposure-interaction term(s) across ", sum(bad), " logged event(s)",
       if (any(wipeout)) " (including a full-predictor wipeout)" else "", ":\n",
       paste(utils::capture.output(print(le[bad, , drop = FALSE])), collapse = "\n"),
       "\n\nFall back to the additive forms 'A0 * A1 * (sex + race + hiqual)' for the ",
       "outcome and 'y * A_other * (sex + race + hiqual)' for each exposure ",
       "(design 8.4 rung 2) and state in the paper that the imputation model is ",
       "congenial for the additive part of heterogeneity only.")
}

run_mice <- function(wide_data, m = 10, maxit = 10, seed = 20260522,
                     outcome = c("MCS", "PCS")) {
  outcome <- match.arg(outcome)

  method_list <- mice::make.method(wide_data)

  # Keyed by group: column name without its _<wave> or _base suffix.
  methods_by_group <- c(
    # Time-invariant confounders
    sex_dv_fact = "logreg",
    hiqual_dv_fact = "polr",
    race = "logreg",
    # Lagged confounders
    pcs_lagged = "pmm",
    mcs_lagged = "pmm",
    econ_dist_bin_lagged = "logreg",
    econ_benefits_lagged = "logreg",
    home_owner_lagged = "logreg",
    mastat_lagged = "logreg",
    dnc_lagged = "polr",
    econ_emp_bin_fact = "logreg",
    log_income = "pmm",       # continuous
    # Other confounders
    gor_dv_fact = "polyreg",
    age_dv = "norm",
    # Exposure
    econ_dist_bin = "logreg",
    # Outcome
    sf12mcs_dv = "pmm",
    sf12pcs_dv = "pmm"
  )

  group   <- sub("_(\\d+|base)$", "", names(method_list))
  matched <- group %in% names(methods_by_group)
  method_list[matched] <- methods_by_group[group[matched]]

  # ---- Congenial imputation model (design 8) -------------------------------
  # Downstream analyses let the effect differ by levels of S = sex x race x hiqual,
  # so the imputation model must include S too. Otherwise imputed values are drawn
  # under effect homogeneity, biasing GATE estimates towards the null.

  # Final-wave outcome; earlier-wave SF-12 columns are auxiliaries.
  y_all <- grep(paste0("^sf12", tolower(outcome), "_dv_\\d+$"), names(wide_data), value = TRUE)
  y_col <- y_all[which.max(as.integer(sub(".*_", "", y_all)))]
  a_col <- grep("^econ_dist_bin_\\d+$", names(wide_data), value = TRUE)
  S     <- "sex_dv_fact_base * race_base * hiqual_dv_fact_base"

  stopifnot(length(y_col) == 1L, length(a_col) > 0L)

  # One formula per variable, every other column as predictor. formulas mode, not
  # predictorMatrix: supplying both to mice() leaves predictorMatrix silently inert.
  form_list <- mice::make.formulas(wide_data)

  # Y model: saturated A0 x A1 x S -- the outcome that carries the interaction estimand.
  # Separate A_t x S terms would force a zero A0:A1 interaction.
  form_list[[y_col]] <- stats::update.formula(
    form_list[[y_col]], stats::reformulate(c(".", paste(c(a_col, S), collapse = " * ")))
  )

  # A models ("full" congeniality): y x A_other x S, so each exposure is imputed
  # under the same regime x S heterogeneity the Y model allows. y enters centred
  # (same model): raw SF-12 (~50) makes the y:S columns near-collinear with S, and
  # mice's remove.lindep() then drops them.
  y_c <- sprintf("I(%s - %.2f)", y_col, mean(wide_data[[y_col]], na.rm = TRUE))
  form_list[a_col] <- purrr::map(a_col, \(a) stats::update.formula(
    form_list[[a]],
    paste("~ . -", y_col, "+", paste(c(y_c, setdiff(a_col, a), S), collapse = " * "))
  ))

  mids <-  mice::mice(
           data            = wide_data,
           m               = m,
           maxit           = maxit,
           seed            = seed,
           method          = method_list,
           formulas        = form_list,
           # remove.lindep() at its default eps (1e-4) drops nested A0 x A1 x S terms
           # from well-posed models; 1e-8 still drops empty cells and exact duplicates.
           eps             = 1e-8
          )

  le <- mids$loggedEvents
  if (is.null(le) || nrow(le) == 0) {
    message("run_mice: no logged events.")
  } else {
    message("run_mice: ", nrow(le), " logged event(s) during imputation:")
    message(paste(utils::capture.output(print(le)), collapse = "\n"))
  }

  # stop if mice dropped any congenial interaction term
  assert_congenial_terms(mids, a_col = a_col, y_col = y_col)

  mids
}
