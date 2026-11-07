# Diagnostics for a mice run:
#   1. mids$loggedEvents (collinearity or other problems), per variable.
#   2. A per-variable missingness + method + event-count table.
#   3. Optional PNGs: observed-vs-imputed plots (density for continuous,
#      category proportions for factors) and convergence traces.
#
# Standalone: source this file, then call diagnose_mice() on any mids object,
# with one plot folder per dataset:
#   diagnose_mice(mids, data_given_to_mice, plot_dir = "outputs/figs/mice_diag/<dataset>")

diagnose_mice <- function(mids, raw_data, plot_dir = NULL, plot_vars = NULL) {
  stopifnot(inherits(mids, "mids"), is.data.frame(raw_data))
  if (!all(names(raw_data) %in% names(mids$method))) {
    warning("raw_data columns do not all match the data mice saw; ",
            "method / event counts may be incomplete.", call. = FALSE)
  }

  le <- mids$loggedEvents
  le_raw <- le          # kept for the unrecognized-format fallback below
  le_unrecognized <- !is.null(le) &&
    (!is.data.frame(le) || !all(c("dep", "out") %in% names(le)))
  if (le_unrecognized) {
    le <- NULL          # skip the structured breakdown and event counts
  }
  if (!is.null(le)) {
    # initialization events (constant/collinear removals) have dep == ""
    le$dep[le$dep == ""] <- "(setup)"
  }

  ## ---- Section 1: loggedEvents breakdown ---------------------------------
  cat("\n---- mice diagnosis: logged events ----\n")
  if (le_unrecognized) {
    cat("loggedEvents present but in an unrecognized format; raw object:\n")
    print(le_raw)
  } else if (is.null(le)) {
    cat("No logged events \n")
  } else {
    cat(nrow(le), "logged events.\n")
    cat("\nEvents by variable (dep = model that logged the event,",
        "out = term removed):\n")
    print(table(dep = le$dep, out = le$out))
    cat("\nFirst", min(10L, nrow(le)), "raw rows:\n")
    print(utils::head(le, 10L))
  }

  ## ---- Section 2: missingness + method table ------------------------------
  # 'out' entries can be comma-separated lists and, for factors, term labels
  out_tokens <- character(0L)
  if (!is.null(le)) {
    out_tokens <- trimws(unlist(strsplit(as.character(le$out), ",")))
  }
  n_as_dep <- purrr::map_int(names(raw_data),
                             \(v) if (is.null(le)) 0L else sum(le$dep == v))
  n_as_out <- purrr::map_int(names(raw_data), \(v) sum(startsWith(out_tokens, v)))

  summary_table <- data.frame(
    variable        = names(raw_data),
    n_missing       = purrr::map_int(raw_data, \(x) sum(is.na(x))),
    pct_missing     = round(100 * purrr::map_dbl(raw_data, \(x) mean(is.na(x))), 1L),
    method          = unname(mids$method[names(raw_data)]),
    n_logged_as_dep = n_as_dep,
    n_logged_as_out = n_as_out,
    row.names       = NULL
  )
  summary_table <- summary_table[
    order(-summary_table$n_logged_as_dep,
          -summary_table$n_logged_as_out,
          -summary_table$pct_missing), ]

  cat("\n---- mice diagnosis: missingness and methods ----\n")
  cat("(n_logged_as_dep = events in this variable's own imputation model;\n",
      " n_logged_as_out = times it was removed from other models - high\n",
      " values flag collinear variables)\n", sep = "")
  print(summary_table, row.names = FALSE)

  ## ---- Section 3: plots ----------------------------------------------------
  if (!is.null(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    imp_counts <- purrr::map_int(mids$imp, \(x) if (is.null(x)) 0L else nrow(x))
    if (is.null(plot_vars)) {
      plot_vars <- names(imp_counts)[imp_counts > 0L]
    }

    # every plot opens its own png device and closes exactly that one on exit
    save_png <- function(file, make_plot) {
      tryCatch({
        grDevices::png(file.path(plot_dir, file),
                       width = 1600, height = 1200, res = 150)
        dev_id <- grDevices::dev.cur()
        on.exit(grDevices::dev.off(dev_id), add = TRUE)
        p <- make_plot()
        # lattice objects only render under Rscript when printed explicitly
        if (inherits(p, "trellis")) print(p)
      }, error = function(e) {
        warning("Plot failed (", file, "): ", conditionMessage(e),
                call. = FALSE)
      })
    }

    save_png("trace.png", function() plot(mids))

    purrr::walk(plot_vars, \(v) {
      col <- raw_data[[v]]
      if (is.numeric(col)) {
        fml <- stats::as.formula(paste("~", v))
        save_png(paste0("density_", v, ".png"),
                 function() mice::densityplot(mids, fml))
      } else if (is.factor(col)) {
        save_png(paste0("props_", v, ".png"), function() {
          lv  <- levels(col)
          obs <- factor(as.character(col[!is.na(col)]), levels = lv)
          imp <- factor(unlist(purrr::map(mids$imp[[v]], as.character)),
                        levels = lv)
          tab <- rbind(Observed = prop.table(table(obs)),
                       Imputed  = prop.table(table(imp)))
          graphics::barplot(tab, beside = TRUE, las = 2,
                            legend.text = rownames(tab),
                            main = v, ylab = "Proportion")
        })
      }
    })
    cat("\nPlots written to ", plot_dir, "\n", sep = "")
  }

  cat("\nNote: m = ", mids$m, ", maxit = ", mids$iteration,
      " - small test runs are indicative, not confirmatory.\n", sep = "")

  invisible(list(events = le, summary = summary_table))
}
