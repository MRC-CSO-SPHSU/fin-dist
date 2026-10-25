#######################################################################################################
# PROJECT: Financial distress as a mediator of the relationship between employment and health outcomes
# DESCRIPTION: Mock ltmle analysis — empty DT skeleton with time-varying nodes (no data loaded)
# AUTHORS:  Darwin del Castillo
# LAST UPDATE: 07 May 2026
#######################################################################################################

pacman::p_load(data.table, 
               ltmle, 
               SuperLearner)

source(here::here("fnct", "import_cleaning.R"))  # data loading and cleaning functions
source(here::here("fnct", "helpers.R"))          # data manipulation and formatting helpers

# ── 0. Constants ───────────────────────────────────────────────────────────────
T_WAVES <- 3L           # waves to include (subset of 10 for feasibility)
MCS_RANGE <- c(0, 100)  # SF-12 MCS observed range — used to rescale Y to [0,1]
PCS_RANGE <- c(0, 100)  # SF-12 PCS observed range

# ── 1. Column schema ─────────────────────────────────────────────────────────
# Baseline (time-invariant) — measured at wave 1, carried forward
baseline_cols <- c(
  "age_0",        # age at baseline
  "sex_0",        # sex (0 = female, 1 = male)
  "race_0",       # ethnicity (0 = white, 1 = non-white)
  "hiqual_0",     # highest qualification (ordinal)
  "gor_0",        # government office region (1-12)
  "mcs_0",        # SF-12 MCS at baseline (rescaled [0,1])
  "pcs_0"         # SF-12 PCS at baseline (rescaled [0,1])
)

# Time-varying — per wave t: L_t → A_t → M_t → Y_t
# ltmle reads left-to-right column position to infer temporal ordering
time_cols <- unlist(lapply(seq_len(T_WAVES), function(t) {
  c(
    paste0("M_dist_",      t),   # EXPOSURE: financial distress (binary indicator) 
    paste0("Y_mcs_",      t),   # OUTCOME: SF-12 mental component (rescaled [0,1])
    paste0("Y_pcs_",      t)    # OUTCOME: SF-12 physical component (rescaled [0,1])
  )
}))

all_cols <- c(baseline_cols, time_cols)

# ── 3. Node classification ───────────────────────────────────────────────────
# Anodes: the treatment/exposure nodes ltmle intervenes on
Anodes <- grep("^A_emp_",    names(DT), value = TRUE)

# Lnodes: time-varying confounders + mediator (M in Lnodes → total effect)
Lnodes <- c(
  grep("^L_",      names(DT), value = TRUE),
  grep("^M_dist_", names(DT), value = TRUE)
)

# Ynodes: one outcome at a time (ltmle handles a single Y per call)
Ynodes_mcs <- grep("^Y_mcs_", names(DT), value = TRUE)
Ynodes_pcs <- grep("^Y_pcs_", names(DT), value = TRUE)

# ── 4. Intervention rules ────────────────────────────────────────────────────
# Static rules: assign employment status for all T waves
abar_employed   <- rep(0L, T_WAVES)  # always employed  (A_emp = 0)
abar_unemployed <- rep(1L, T_WAVES)  # always unemployed (A_emp = 1)

# ── 5. SuperLearner library ──────────────────────────────────────────────────
SL.lib <- c("SL.glm", "SL.glmnet", "SL.ranger")

# ── 6. ltmle calls ───────────────────────────────────────────────────────────
# NOTE: these calls require a populated DT — structural skeleton only

# Mental health outcome
result_mcs <- ltmle(
  data          = DT,
  Anodes        = Anodes,
  Lnodes        = Lnodes,
  Ynodes        = Ynodes_mcs,
  abar          = list(abar_employed, abar_unemployed), 
  SL.library    = NULL, # glm will be used for pipeline
  estimate.time = TRUE
)

# Physical health outcome
# result_pcs <- ltmle(
#  data          = DT,
#  Anodes        = Anodes,
#  Lnodes        = Lnodes,
#  Ynodes        = Ynodes_pcs,
#  abar          = list(abar_employed, abar_unemployed),
#  SL.library    = SL.lib,
#  estimate.time = FALSE
#)

# ── 7. Results ───────────────────────────────────────────────────────────────
ltmle::summary(result_mcs)
ltmle::summary(result_pcs)