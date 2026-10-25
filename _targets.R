# _targets.R — pipeline orchestration for the financial distress analysis from the UKHLS data.
# Parallel backend: future + future.batchtools.
# See documentation for more details

library(targets)
library(tarchetypes)
library(future)
library(future.batchtools)
library(future.callr)

# Detect SLURM at runtime, if not in cluster, run locally with future.callr (in a separate r process)
on_slurm <- nzchar(Sys.getenv("SLURM_JOB_ID")) && nzchar(Sys.which("sbatch"))

if (on_slurm) {
  future::plan(
    future.batchtools::batchtools_slurm,
    template  = "slurm.tmpl",
    resources = list(
      ncpus    = 2,
      memory   = 24 * 1024,      # MB per CPU (= 24 GB per worker).
      walltime = 6 * 60 * 60,    # seconds (= 6 h wall)
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
    "tmle",
    "miceadds",
    "SuperLearner",
    "xgboost",
    "gam",
    "arm",
    "quarto"
  ),
  format = "rds",
  seed   = 20261022
)

# ---- Sourcing functions (R/) ----
for (f in c(list.files("R", "\\.R$", full.names = TRUE))) source(f); rm(f)

# --------------------- Configuration -------------------------
## mice configs
mice_m      <- 50
mice_maxit  <- 20
seed_random <- 20261022
## gFormulaMI configs
gform_M <- 50
cate_gbound <- 0.001 # non-aggressive, symmetric bounding
tmle_bound <- 0.001 # non-aggressive symmetric bounding

## Analytical samples for the whole pipeline and outcomes

specs <- tidyr::expand_grid(
  tibble::tibble(label = c("w3", "w8", "w13"),
  target_wave = c(3L, 8L, 13L),
  pre_wave = c(2L, 7L, 12L),
  base_wave = c(1L, 6L, 11L)),
  outcome = c("MCS", "PCS")
)

## SuperLearner library for the TMLE Q- and g-models. SL.xgboost.tmle is the
## custom wrapper in R/sl_wrappers.R (re-scaled outcome (0,1) handling).
sl_libs <- c("SL.mean", "SL.glm", "SL.poly2", "SL.poly3", "SL.gam",
             "SL.nnet.tmle", "SL.rf.tmle", "SL.xgboost.tmle")

# ------------------------------- DAG -----------------------------------
# ---- Loop for each outcome scale (MCS, PCS) using static branching ----
map_pipe <- tar_map(
  values = specs,
  names  = c(label, outcome),
  tar_target(wide_data,
             build_data(pop_data, 
                        outcome = outcome,
                        pre_wave = pre_wave,
                        target_wave = target_wave,
                        base_wave = base_wave)),
  tar_target(wide_mids,
             run_mice(wide_data,
                      m       = mice_m,
                      maxit   = mice_maxit,
                      seed    = seed_random,
                      outcome = outcome)),
  tar_target(
    tmle_one,
    fit_tmle_one(
      wide_mids = wide_mids,
      imp_idx   = tmle_imp_idx,
      sl_libs   = sl_libs,
      outcome   = outcome,
      gbound    = tmle_bound
    ),
    pattern   = map(tmle_imp_idx), # one branch per imputation
    iteration = "list"             # collect the per-imputation tmle fits in a list
  ),
  tar_target(tmle_results,
             pool_tmle(tmle_one)),
  # Sensitivity 1: TMLE with gbound = 0.025
  tar_target(tmle_sens1,
             fit_tmle_one(
               wide_mids = wide_mids,
               imp_idx   = tmle_imp_idx,
               sl_libs   = sl_libs,
               outcome   = outcome,
               gbound    = 0.025
             ),
             pattern   = map(tmle_imp_idx), # one branch per imputation
             iteration = "list"             # collect the per-imputation tmle fits in a list
  ),
  tar_target(tmle_sens1_results,
             pool_tmle(tmle_sens1)),
  # Sensitivity 2: G-formula with the same estimand as TMLE (factor(regime) + 0)
  tar_target(mi_results,
             run_gformula(
               wide_mids     = wide_mids,
               wide_data_mi  = wide_data,
               M             = gform_M,
               estimand      = "factor(regime) + 0",
               outcome_scale = outcome
             )),
  # Sensitivity 3: G-formula with the estimand of interest (factor(regime))
  tar_target(mi_ate_results,
             run_gformula(
               wide_mids     = wide_mids,
               wide_data_mi  = wide_data,
               M             = gform_M,
               estimand      = "factor(regime)",
               outcome_scale = outcome
             )),

  tar_target(comparison,
             assemble_comparison(tmle_results,
                                 mi_results,
                                 mi_ate_results
            )),

  # ---------- Heterogeneous Treatment Effects: DR-learner CATEs ---------------
  # Reuses tmle_imp_idx so branching stays aligned with tmle_one.
  tar_target(
    cate_one,
    estimate_cate(
      wide_mids = wide_mids,
      imp_idx   = tmle_imp_idx,
      sl_libs   = sl_libs,
      outcome   = outcome,
      gbound    = cate_gbound
    ),
    pattern   = map(tmle_imp_idx), # one branch per imputation
    iteration = "list"             # collect the per-imputation results in a list
  ),
  tar_target(cate_results,
             pool_cate(cate_one)),
  # Sensitivity: DR-learner CATEs with gbound = 0.025
  tar_target(
    cate_sens1_one,
    estimate_cate(
      wide_mids = wide_mids,
      imp_idx   = tmle_imp_idx,
      sl_libs   = sl_libs,
      outcome   = outcome,
      gbound    = 0.025
    ),
    pattern   = map(tmle_imp_idx), # one branch per imputation
    iteration = "list"             # collect the per-imputation results in a list
  ),
  tar_target(cate_sens1_results,
             pool_cate(cate_sens1_one))
)

# --- DAG ---
list(
  # Data prep: import/clean/preproc
  tar_target(pop_data,
    if (on_slurm) {
      import_data(force = FALSE) |> clean_data() |> preproc_data()
    } else {
      import_data(force = TRUE) |> clean_data() |> preproc_data()
    },
    cue = tar_cue(mode = "always"),
    description = "Data importing and formatting"
  ),

  # Imputation index, shared by both outcome scales
  tar_target(tmle_imp_idx, seq_len(mice_m)),

  map_pipe,

  tar_combine(tmle_results_all, map_pipe[["tmle_results"]],
              command = dplyr::bind_rows(!!!.x, .id = "target")),

  # Single report covering both outcomes
  tarchetypes::tar_quarto(report, "07_single_treatment.qmd")
)