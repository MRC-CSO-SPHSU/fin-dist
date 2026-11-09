# targets pipeline for the UKHLS financial distress analysis.

library(targets)
library(tarchetypes)
library(future)
library(future.batchtools)
library(future.callr)

# Parallel backend picked at runtime: future.batchtools on SLURM, else
# future.callr (each target in its own R process).
on_slurm <- nzchar(Sys.getenv("SLURM_JOB_ID")) && nzchar(Sys.which("sbatch"))

if (on_slurm) {
  future::plan(
    future.batchtools::batchtools_slurm,
    template  = "slurm.tmpl",
    resources = list(
      ncpus    = 2,
      memory   = 24 * 1024,      # MB per CPU (slurm.tmpl: --mem-per-cpu), so 48 GB per worker
      walltime = 6 * 60 * 60,    # seconds (6 h)
      account  = "none"
    )
  )
} else {
  future::plan(future.callr::callr, workers = 4)
}

message("--- TARGETS FUTURE PLAN ---")
message("hostname:        ", Sys.info()[["nodename"]])
message("SLURM_JOB_ID:    '", Sys.getenv("SLURM_JOB_ID"), "'")
message("Sys.which sbatch: '", Sys.which("sbatch"), "'")
message("on_slurm:        ", on_slurm)
message("plan:            ", paste(class(future::plan()), collapse = "/"))
message("---------------------------")

# ---- Packages attached to every target's evaluation environment ----
tar_option_set(
  packages = c(
    "data.table",
    "bit64",
    "dplyr",
    "tidyr",
    "tibble",
    "purrr",
    "magrittr",
    "rlang",
    "here",
    "mice",
    "gFormulaMI",
    "estimatr",
    "car",
    "lmtest",
    "sandwich",
    "ltmle",
    "miceadds",
    "SuperLearner",
    "xgboost",
    "gam",
    "arm",
    "quarto"
  ),
  format = "rds",
  repository = "local",
  storage = "worker",
  retrieval = "worker",
  seed   = 20261022,
  # To avoid excess of calls to the cluster for unresolved workers from the controller.
  backoff = tar_backoff(min = 10, max = 60, rate = 2)
)

# ---- Sourcing functions (R/) ----
purrr::walk(list.files(here::here("R"), "\\.R$", full.names = TRUE), source)

# --------------------- Configuration -------------------------
## mice configs
mice_m      <- 50
mice_maxit  <- 20
seed_random <- 20261022
## gFormulaMI: gFormulaImpute() uses mids$m (= mice_m) when given a mids, so M is ignored.
gform_M <- 50
cate_gbound <- 0.001 # non-aggressive, symmetric bounding
## ltmle lower bound (non-aggressive) on each regime's cumulative g,
ltmle_bound <- c(0.001, 1)
sen_ltmle_bound <- c(0.025, 1)

## Analytic samples: one row per window and outcome. tar_map turns "w3-4" into
## "w3.4" in target names.
specs <- tidyr::expand_grid(
  tibble::tibble(
    label      = c("w3-4", "w8-9", "w12-13"),
    pre_wave   = c(1L, 6L, 10L),
    base_wave  = c(2L, 7L, 11L),
    start_wave = c(3L, 8L, 12L),
    final_wave = c(4L, 9L, 13L)
  ),
  outcome = c("MCS", "PCS")
)

## c(distress at start_wave, distress at final_wave). The first regime is the
## reference for contrasts in ltmle and for gFormulaMI's factor(regime).
ltmle_regimes <- list(never      = c(0L, 0L),
                      before     = c(1L, 0L),
                      now        = c(0L, 1L),
                      permanent  = c(1L, 1L))

## SuperLearner's internal V-fold CV for the ensemble weights of every ltmle Q/g
## fit, not TMLE cross-fitting (ltmle has none). ltmle accepts only
## V/stratifyCV/shuffle; stratifyCV must stay FALSE as the Q outcomes are continuous.
sl_cv_control <- list(V = 5L,
                      stratifyCV = FALSE,
                      shuffle = TRUE)

## SuperLearner library for the ltmle Q- and g-models. SL.poly2/3 and SL.*.ltmle
## are the custom wrappers in R/sl_wrappers.R.
sl_libs <- c("SL.mean", "SL.glm", "SL.poly2", "SL.poly3", "SL.gam",
             "SL.nnet.ltmle", "SL.rf.ltmle", "SL.xgboost.ltmle")

# ------------------------------- DAG -----------------------------------
# ---- Static branching: one copy of these targets per specs row ----
map_pipe <- tar_map(
  values = specs,
  names  = c(label, outcome),
  tar_target(wide_data,
             build_data(pop_data,
                        outcome    = outcome,
                        pre_wave   = pre_wave,
                        base_wave  = base_wave,
                        start_wave = start_wave,
                        final_wave = final_wave)),
  tar_target(wide_mids,
             run_mice(wide_data,
                      m       = mice_m,
                      maxit   = mice_maxit,
                      seed    = seed_random,
                      outcome = outcome)),
  tar_target(
    ltmle_one,
    fit_ltmle_one(
      wide_mids     = wide_mids,
      imp_idx       = ltmle_imp_idx,
      sl_libs       = sl_libs,
      outcome       = outcome,
      regimes       = ltmle_regimes,
      sl_cv_control = sl_cv_control,
      gbound        = ltmle_bound
    ),
    pattern   = map(ltmle_imp_idx), # one branch per imputation
    iteration = "list"              # pool_ltmle() needs the list of per-imputation tables
  ),
  tar_target(ltmle_results,
             pool_ltmle(ltmle_one)),
  # Sensitivity 1: LTMLE with gbound = 0.025 on the cumulative g
  tar_target(ltmle_sens1,
             fit_ltmle_one(
               wide_mids     = wide_mids,
               imp_idx       = ltmle_imp_idx,
               sl_libs       = sl_libs,
               outcome       = outcome,
               regimes       = ltmle_regimes,
               sl_cv_control = sl_cv_control,
               gbound        = sen_ltmle_bound
             ),
             pattern   = map(ltmle_imp_idx),
             iteration = "list"
  ),
  tar_target(ltmle_sens1_results,
             pool_ltmle(ltmle_sens1)),
  # Sensitivity 2: G-formula regime means (factor(regime) + 0)
  tar_target(mi_results,
             run_gformula(
               wide_mids     = wide_mids,
               wide_data_mi  = wide_data,
               M             = gform_M,
               estimand      = "factor(regime) + 0",
               outcome_scale = outcome,
               regimes       = ltmle_regimes
             )),
  # Sensitivity 3: G-formula contrasts vs never (factor(regime))
  tar_target(mi_ate_results,
             run_gformula(
               wide_mids     = wide_mids,
               wide_data_mi  = wide_data,
               M             = gform_M,
               estimand      = "factor(regime)",
               outcome_scale = outcome,
               regimes       = ltmle_regimes
             )),

  tar_target(comparison,
             assemble_comparison(ltmle_results,
                                 mi_results,
                                 mi_ate_results
            )),

  # ---------- Heterogeneous Treatment Effects: DR-learner CATEs ---------------
  # Still single-exposure (econ_dist_bin_0, sf12*_0) until the CATE scripts are
  # adapted to transitions. Reuses ltmle_imp_idx so its branches line up with ltmle_one.
  tar_target(
    cate_one,
    estimate_cate(
      wide_mids = wide_mids,
      imp_idx   = ltmle_imp_idx,
      sl_libs   = sl_libs,
      outcome   = outcome,
      gbound    = cate_gbound
    ),
    pattern   = map(ltmle_imp_idx),
    iteration = "list"
  ),
  tar_target(cate_results,
             pool_cate(cate_one)),
  # Sensitivity: DR-learner CATEs with gbound = 0.025
  tar_target(
    cate_sens1_one,
    estimate_cate(
      wide_mids = wide_mids,
      imp_idx   = ltmle_imp_idx,
      sl_libs   = sl_libs,
      outcome   = outcome,
      gbound    = 0.025
    ),
    pattern   = map(ltmle_imp_idx),
    iteration = "list"
  ),
  tar_target(cate_sens1_results,
             pool_cate(cate_sens1_one))
)

# --- DAG ---
list(
  tar_target(pop_data,
    if (on_slurm) {
      import_data(force = FALSE) |> clean_data() |> preproc_data()
    } else {
      import_data(force = TRUE) |> clean_data() |> preproc_data()
    },
    cue = tar_cue(mode = "always"),
    description = "Data importing and formatting"
  ),

  # Imputation index, shared by every specs row
  tar_target(ltmle_imp_idx, seq_len(mice_m)),

  map_pipe,

  tar_combine(ltmle_results_all, map_pipe[["ltmle_results"]],
              command = dplyr::bind_rows(!!!.x, .id = "target")),

  # Single report covering both outcomes
  tarchetypes::tar_quarto(report, "07_single_treatment.qmd")
)