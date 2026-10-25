# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

A causal-inference study estimating the effect of **financial distress (`econ_dist_bin`) on mental health (`sf12mcs_dv`)** using UKHLS (Understanding Society) data. The repository currently carries **two methodological tracks on two branches** (see "Branch state" below):

- **Longitudinal track (`main`)** — waves 3–5, time-varying treatment. Headline estimator is **LTMLE** (longitudinal TMLE), with **g-formula via multiple imputation (`gFormulaMI`)** and **IPTW** as sensitivity analyses.
- **Single-point track (`exp/single-treat`, the current branch)** — reframes the question as a **point treatment**: the effect of financial distress at one wave on mental health at a later wave, estimated with **TMLE** (`tmle::tmle()`), with **g-formula via `gfoRmula`** as the sensitivity analysis. Widens the data window to **waves 3–10** to gain person-time.

Missing data is handled with **`mice`** throughout, and everything is orchestrated as a **`targets`** pipeline. Language is R throughout; there is no build/lint toolchain — the unit of work is the targets DAG (a `testthat` suite lives in `tests/`, see below).

## Branch state — which analysis is live

> **`CLAUDE.md` is gitignored**, so this file does **not** change when you `git checkout` between branches. The same file describes both tracks; read this section to know which one applies to your working tree. Check the branch with `git branch --show-current`.

- **`main` = "as usual" longitudinal LTMLE.** The pipeline documented in "Longitudinal pipeline" below is the source of truth there. Do not port single-point changes onto it unless asked.
- **`exp/single-treat` = the transition to single-point TMLE.** This is where the LTMLE→TMLE and `gFormulaMI`→`gfoRmula` migration happens. **It is mid-flight**: new functions exist (`R/fit_tmle_one.R`, `R/build_data.R`, `R/run_mice_long.R`) and `07_single_treatment.qmd` is the new report, but `_targets.R` has **not** been fully rewired (it still builds the LTMLE DAG, renders `06_full_models.qmd`, and its `tar_option_set(packages=...)` still lists `gFormulaMI`/`ltmle` but not `gfoRmula`/`tmle`/`miceadds`). Treat the single-point DAG as **intended, not yet wired** — see "Single-point pipeline (in progress)".

## Project plans & decisions (single source of truth)

Implementation plans and design decisions for ongoing work live in **`.claude/plans/`** (project-local, gitignored). That folder is the **single source of truth** for what we are building and why: consult it before starting work, and keep it current as plans are approved or revised.

**Codebase context** (what the code *actually does today*, as opposed to what the plans intend) lives in **`.claude/context/`**:

- [`.claude/context/codebase-map.md`](.claude/context/codebase-map.md) — full function inventory, DAG wiring, variable definitions, and a list of the signature mismatches that currently stop `tar_make()` from running end-to-end.
- [`.claude/context/hierarchical-tmle-balzer.md`](.claude/context/hierarchical-tmle-balzer.md) — API and hard data requirements of the vendored Balzer cluster-TMLE code in `R/helpers_cluster_tmle.R`, plus the gap analysis for wiring it into the MAIHDA target.

- [`.claude/plans/point-treatment-dr.md`](.claude/plans/point-treatment-dr.md) — the original point-treatment doubly-robust plan. **Note it is now partially superseded**: the live implementation on `exp/single-treat` uses the **`tmle` package** (not a one-A-node `ltmle::ltmle()` call), drops the IPTW sensitivity arm, moves to **long-format multilevel imputation**, and uses **`gfoRmula`** for the g-formula arm. Update this plan to match the realized design before relying on it.

## Running the pipeline

```r
# From the project root, in R:
targets::tar_make_future(workers = 8L)   # primary entrypoint (parallel)
targets::tar_make()                       # serial fallback

targets::tar_visnetwork()                 # inspect the DAG / what's stale
targets::tar_read(comparison)             # pull a computed target into the session
targets::tar_read(ltmle_results)          # (main) longitudinal estimates
targets::tar_read(tmle_results)           # (exp/single-treat) point estimates, once wired
targets::tar_outdated()                   # which targets would rerun
```

On the Glasgow SLURM cluster, submit the whole thing as a batch job:

```bash
sbatch run_quarto05.sbatch    # activates conda `quarto` env, runs tar_make_future(workers=8L)
```

The final `report` target renders the current report via `tar_quarto()`: **`06_full_models.qmd`** on `main`, transitioning to **`07_single_treatment.qmd`** on `exp/single-treat`.

### Tests

```r
testthat::test_file("tests/test_code.R")   # structural smoke tests for the data + mice functions
```

`tests/test_code.R` checks that the data builders and `mice` calls run and produce the expected columns / a valid `mids`; it also drives `tests/diagnose_mice.R` (convergence plots). On `exp/single-treat` some tests still `source()` deleted files (`R/build_wide_data.R`) and stale `fnct/` paths — update them alongside the migration.

## Longitudinal pipeline (`_targets.R` on `main`)

The DAG, in dependency order:

1. `pop_data` — `import_data() |> clean_data() |> preproc_data()` (functions in `R/import_cleaning.R`). `import_data()` reads raw UKHLS `.dta` files for **all 10 waves** (`letters[1:10]`) and **caches to `data/cache/raw_data.fst`**; pass `force = TRUE` to rebuild from source.
2. `wide_data` — reshape long → wide with one column per (variable × wave), e.g. `econ_dist_bin_0/1/2`. Uses `make_wide()` in `R/helpers.R`.
3. `wide_mids` — `run_mice()`, the single wide-format multiple-imputation step (one cached `mids` object).
4. `ltmle_data_list` → `ltmle_one` → `ltmle_results` — the main estimator. `ltmle_one` **dynamically branches over the (regime × imputation) grid** (`work_grid`); each branch is one `fit_ltmle_one()` call. `pool_ltmle()` combines per-imputation fits with Rubin's rules (`mice::pool.scalar`).
5. `mi_results` (`run_gformula`) and the IPTW arm — sensitivity analyses; depend only on `wide_mids` and run in parallel with the LTMLE branch.
6. `comparison` — `assemble_comparison()` joins the estimators on `intervention`.
7. `report` — renders the report `.qmd`.

### Key design facts (read before editing the longitudinal track)

- **`R/` holds every pipeline function** (data helpers included — there is **no `fnct/` directory** anymore; older docs and some `tests/` `source()` calls still reference `fnct/import_cleaning.R` / `fnct/helpers.R`, which are now `R/import_cleaning.R` / `R/helpers.R`). The statistical code in `R/` was extracted *verbatim* from the original `.qmd` chunks — preserve behavior exactly unless explicitly changing methodology. The report reads results with `tar_read()` rather than recomputing.
- **The three treatment timepoints are `_0/_1/_2`** = UKHLS waves 3/4/5. Node roles in LTMLE: `econ_dist_bin_*` = A (treatment), `sf12mcs_dv_*` = Y (outcome), `pcs_lagged_*` = L (time-varying confounder), `*_base` = baseline covariates. `Qform`/`gform`/`regimes` in `_targets.R` encode these.
- **`fit_ltmle_one()` does `assign("SL.xgboost.ltmle", ..., envir = globalenv())`** on purpose: `ltmle`/`tmle` resolve `SL.library` learner names by string against the global env, and on a fresh batchtools worker the custom wrapper (`R/sl_wrappers.R`) wouldn't otherwise be present. Don't "clean this up." The same shim is required in `fit_tmle_one()`.
- **`SL.xgboost.ltmle`** wraps `SL.xgboost` to handle quirks: factor A/L nodes coded 1/2 (coerced to 0/1) and bounded-continuous Q-models that need `gaussian()` instead of `binomial`.
- **Config lives at the top of `_targets.R`** as dev-vs-final settings: `mice_m` (5 dev / 35 final), `mice_maxit` (10 / 15), `gformula_M` (20 / 50). Flip these for production runs. Seed is `20260522`.

## Single-point pipeline (in progress, `exp/single-treat`)

The reframing collapses the longitudinal estimand to a **point exposure**, so there is no treatment–confounder feedback and standard conditional-ignorability TMLE applies. The estimator API changes accordingly — `tmle::tmle()` takes separate `A` / `W` / `Y` arguments (not a node-ordered frame with `abar`), and operates on **long-format** columns with *no wave suffix* (`econ_dist_bin`, `sf12mcs_dv`) rather than the `_0/_1/_2` wide columns.

New / changed building blocks (functions exist; wiring into `_targets.R` is the remaining work):

- **`R/build_data.R`** (`build_data()`, replaces the deleted `build_wide_data.R`) — returns a `list(wide_data, long_data)` from `pop_data`. The `datasets` target wraps it. `wide_data` keeps the wide LTMLE shape; `long_data` is one row per person-wave for the multilevel imputation + point-treatment fit.
- **`R/run_mice_long.R`** (`run_mice_long()`) — long-format multilevel MI via **`miceadds`** `2l.pmm` samplers, clustered on `pidp` (`predictorMatrix[, "pidp"] <- -2`). Factor levels are `make.names()`-sanitized first because the `2l.*` samplers paste design-matrix column names into an `lme4` formula, so spaces/hyphens/`+` in labels would parse as operators and crash. Intended to back a `long_mids` target.
- **`R/fit_tmle_one.R`** (`fit_tmle_one()`) — single-point TMLE via `tmle::tmle(A = econ_dist_bin, W = c(...lagged covariates...), Ynodes = "sf12mcs_dv", ...)`. Intended to back a `tmle_results` target.
- **`R/run_gformula.R`** — g-formula arm migrating from `gFormulaMI` to the **`gfoRmula`** package (both are currently loaded in `07_single_treatment.qmd` during the transition; `gfoRmula` is the target).
- **Removed on this branch:** the IPTW sensitivity arm (`R/run_iptw.R`, `R/extract_iptw.R`) and `R/build_wide_data.R`.

When wiring this up, expect to add `datasets`, `long_mids`, and `tmle_results` targets, add `gfoRmula`/`tmle`/`miceadds` to `tar_option_set(packages=...)`, point `run_mice()`'s input at `datasets$wide_data`, and switch the `tar_quarto()` report to `07_single_treatment.qmd`.

## Parallel backend (both tracks)

`_targets.R` selects the backend at runtime: if `SLURM_JOB_ID` is set **and** `sbatch` is on PATH, it uses `future.batchtools::batchtools_slurm` with `slurm.tmpl` (each worker target → its own SLURM job, resources set via `future::plan(resources=...)`). Otherwise it falls back to `future.callr` with 4 local workers. Per-target resource overrides go through `tar_resources(future = tar_resources_future(resources = list(...)))`. Note `import_data()` runs with `force = FALSE` (use the `.fst` cache) on SLURM and `force = TRUE` (rebuild from raw `.dta`) locally.

## Repo layout notes

- **`01_*`–`04_*.qmd`** are the exploratory/development history (simulation tests, simple → full model build-up).
- **`05_imputation.qmd` / `06_full_models.qmd`** drove earlier stages of the pipeline; **`07_single_treatment.qmd`** is the new report for the single-point track. When changing methods, the source of truth is `R/` plus the highest-numbered report `.qmd` for the branch you are on — not the earlier numbered files.
- `_targets/` is the cache (gitignored). `data/` (raw UKHLS `.dta`, simulated `.dta`, fst cache) is gitignored.
- `code/`, `legacy/` — older scripts and the Stata data simulator; not part of the live pipeline.
- `figs/dag/` — DAGs documenting the assumed causal structure (the single-point DAG is also drawn inline in `07_single_treatment.qmd`).
- `.claude/plans/` — **project plans & decisions, the single source of truth for ongoing work** (gitignored; see the section above).
- `.claude/methods_papers/` — reference PDFs (Thomson et al. 2023 is the methodological template).
- `CLAUDE.md` itself is gitignored (and, as noted, shared across branches).
