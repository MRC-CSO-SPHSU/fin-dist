# Congeniality guard for the imputation model.
# The reason why this function exist is related to the interaction term introduced
# to preserve congeniality in the imputation model. This is explained in the next function.
# If mice drops one of the interaction terms for collinearity it records the fact in
# loggedEvents and carries on silently, which is precisely the attenuation this
# change exists to prevent.

assert_congenial_terms <- function(mids, a_col) {
  le <- mids$loggedEvents
  if (is.null(le) || !is.data.frame(le) || nrow(le) == 0L) return(invisible(TRUE))
  if (!"out" %in% names(le)) return(invisible(TRUE))

  keys <- c(a_col)

  dropped <- purrr::map(strsplit(le$out, ",\\s*"), \(tm) {
    purrr::keep(trimws(tm), \(term)
      grepl(":", term, fixed = TRUE) &&
        any(purrr::map_lgl(keys, \(k) grepl(k, term, fixed = TRUE))))
  })

  wipeout <- grepl("All predictors are constant", le$out, fixed = TRUE)

  bad <- lengths(dropped) > 0L | wipeout
  if (!any(bad)) return(invisible(TRUE))

  stop("run_mice: congeniality broken -- mice dropped ", sum(lengths(dropped)),
       " exposure-interaction term(s) across ", sum(bad), " logged event(s)",
       if (any(wipeout)) " (including a full-predictor wipeout)" else "", ":\n",
       paste(utils::capture.output(print(le[bad, , drop = FALSE])), collapse = "\n"),
       "\n\nFall back to the additive form 'A * (sex + race + hiqual)' ",
       "(design 8.4 rung 2) and state in the paper that the imputation model is ",
       "congenial for the additive part of heterogeneity only.")
}

run_mice <- function(wide_data, m = 10, maxit = 10, seed = 20260522,
                     outcome = c("MCS", "PCS")) {

  outcome <- match.arg(outcome)

  method_list <- mice::make.method(wide_data)

  methods_by_group <- c(
    # Time-invariant confounders
    sex_dv_fact_base = "logreg",
    hiqual_dv_fact_base = "polr",
    race_base = "logreg",
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
  # so the imputation model must include S too, otherwise imputed values are drawn 
  # under effect homogeneity and bias towards the null GATE estimates.

  y_col <- if (outcome == "MCS") "sf12mcs_dv_0" else "sf12pcs_dv_0"  # wave-3 outcome only
  a_col <- "econ_dist_bin_0"
  S     <- "sex_dv_fact_base * race_base * hiqual_dv_fact_base"

  stopifnot(all(c(y_col, a_col) %in% names(wide_data)))

  # pred_mat exist only as input of make.formulas(), which creates formulas for imputing all DT vars
  # Supplying both would leave predictorMatrix silently inert.
  pred_mat  <- mice::make.predictorMatrix(wide_data)
  form_list <- mice::make.formulas(wide_data, predictorMatrix = pred_mat) # Creates the formula for all variables in the DT

  # Y model: saturated A x S -- the outcome that carries the interaction estimand.
  form_list[[y_col]] <- stats::update.formula(
    form_list[[y_col]], stats::reformulate(c(".", paste(a_col, "*", S)))
  )

  mids <-  mice::mice(
           data            = wide_data,
           m               = m,
           maxit           = maxit,
           seed            = seed,
           method          = method_list,
           formulas        = form_list
          )

  le <- mids$loggedEvents
  if (is.null(le) || nrow(le) == 0) {
    message("run_mice: no logged events.")
  } else {
    message("run_mice: ", nrow(le), " logged event(s) during imputation:")
    message(paste(utils::capture.output(print(le)), collapse = "\n"))
  }

  # Confirming that the interaction terms are still present in the imputation model
  assert_congenial_terms(mids, a_col = a_col)

  mids
}
