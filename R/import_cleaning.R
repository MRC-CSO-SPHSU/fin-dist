import_data <- function(force = FALSE) {
  cache_path <- here::here("data", "cache", "raw_data.fst")
  if (!force && file.exists(cache_path)) {
    message("Loading from cache. Use force = TRUE to rebuild.")
    return(fst::read_fst(cache_path, as.data.table = TRUE))
  }

  message("Building raw_data.fst and writing to cache.")
  
  # Defining paths
  ukhls_raw <- here::here("data", "raw", "ukhls")
  
  library(data.table)
  library(haven)
  library(purrr)

  waves <- letters[1:13]

  # loading variables from indall
  ## variables to keep
  vars_indall <- c("pidp", "hidp", "hhorig", "fnspid", "mnspid", "ivfio", "ppid",
                   "pns1pid", "pns2pid", "age_dv", "sex_dv", "gor_dv",
                   "depchl_dv", "intdaty_dv", "intdatm_dv", "intdatd_dv",
                   "mastat_dv", "nchild_dv", "ethn_dv")
  indall <- rbindlist(map(seq_along(waves), function(i) {
    w <- waves[i]

    # build prefix names
    wave_vars <- setdiff(vars_indall, "pidp")
    cols <- c("pidp", paste0(w, "_", c(wave_vars)))

    dt <- read_dta(paste0(ukhls_raw, "/", w, "_indall.dta"), col_select = all_of(cols)) |> 
      zap_labels() |> 
      as.data.table()

    # strip wave prefix from variable names
    prefix <- paste0(w, "_")
    setnames(dt, \(x) ifelse(startsWith(x, prefix), substr(x, nchar(prefix) + 1, nchar(x)), x))

    dt$wave <- i
    return(dt)
  }), fill = TRUE)
  indall[, pidp := bit64::as.integer64(pidp)]

  # loading variables from indresp
  vars_indresp <- c("pidp", "hidp", "jbstat", "scsf1", "sclfsato",
                   "fimnlabgrs_dv", "fimnpen_dv", "finnow",
                   "hiqual_dv", "sf12mcs_dv", "sf12pcs_dv")
  indresp <- rbindlist(map(seq_along(waves), function(i) {
    w <- waves[i]

    ## build prefix names (scsf1 not available in wave 1)
    wave_vars <- setdiff(vars_indresp, "pidp")
    if (i == 1) wave_vars <- setdiff(wave_vars, "scsf1")
    cols <- c("pidp", paste0(w, "_", wave_vars))
    dt <- read_dta(paste0(ukhls_raw, "/", w, "_indresp.dta"), col_select = all_of(cols)) |> 
      zap_labels() |> 
      as.data.table()
    
    ## strip wave prefix from variable names
    prefix <- paste0(w, "_")
    setnames(dt, \(x) ifelse(startsWith(x, prefix), substr(x, nchar(prefix) + 1, nchar(x)), x))
    dt$wave <- i
    return(dt)
  }), fill = TRUE)
  indresp[, pidp := bit64::as.integer64(pidp)]

  # loading variables from hhrep
  vars_hhrep <- c("hidp", "fihhmnsben_dv", "hsownd", "tenure_dv")
  hhrep <- rbindlist(map(seq_along(waves), function(i) {
    w <- waves[i]
    ## build prefix names
    cols <- c(paste0(w, "_", c(vars_hhrep)))
    dt <- read_dta(paste0(ukhls_raw, "/", w, "_hhresp.dta"), col_select = all_of(cols)) |>
      zap_labels() |> 
      as.data.table()
    ## strip wave prefix from variable names
    prefix <- paste0(w, "_")
    setnames(dt, \(x) ifelse(startsWith(x, prefix), substr(x, nchar(prefix) + 1, nchar(x)), x))
    dt$wave <- i
    return(dt)
  }), fill = TRUE)

  # loading variables from income
  vars_income <- c("pidp", "hidp", "frmnthimp_dv", "ficode")
  income <- rbindlist(map(seq_along(waves), function(i) {
    w <- waves[i]
    cols <- c("pidp", paste0(w, "_", setdiff(vars_income, "pidp")))

    dt <- read_dta(paste0(ukhls_raw, "/", w, "_income.dta"), col_select = all_of(cols)) |>
      zap_labels() |>
      as.data.table()

    prefix <- paste0(w, "_")
    setnames(dt, \(x) ifelse(startsWith(x, prefix), substr(x, nchar(prefix) + 1, nchar(x)), x))

    dt$wave <- i
    dt
  }), fill = TRUE)
  income[, pidp := bit64::as.integer64(pidp)]

  ## income component splits (ficode to inc_* columns)
  income[, inc_pp  := fifelse(ficode == 4,  frmnthimp_dv, NA_real_)]
  income[, inc_tu  := fifelse(ficode == 25, frmnthimp_dv, NA_real_)]
  income[, inc_ma  := fifelse(ficode == 26, frmnthimp_dv, NA_real_)]
  income[, inc_fm  := fifelse(ficode == 27, frmnthimp_dv, NA_real_)]
  income[, inc_oth := fifelse(ficode == 38, frmnthimp_dv, NA_real_)]

  ## collapse to one row per person-wave: sum income components
  income_collapsed <- income[,
    map(.SD, \(x) sum(x, na.rm = TRUE)),
    by = .(pidp, hidp, wave),
    .SDcols = c("inc_pp", "inc_tu", "inc_ma", "inc_fm", "inc_oth")
  ]

  ## benefits receipt: UC and legacy benefits — collapse by household-wave (max)
  income[, benefits_uc := fifelse(ficode == 40, 1L, 0L)]
  income[, benefits_lb := fifelse(ficode %in% c(15L, 16L, 19L, 20L, 22L, 33L), 1L, 0L)]

  benefits_collapsed <- income[,
    .(benefits_uc = max(benefits_uc), benefits_lb = max(benefits_lb)),
    by = .(hidp, wave)
  ]

  # Merge datasets
  raw_data <- merge(indall, indresp, by = c("pidp", "hidp", "wave"), all.x = TRUE)
  raw_data <- merge(raw_data, hhrep, by = c("hidp", "wave"), all.x = TRUE)
  raw_data <- merge(raw_data, income_collapsed, by = c("pidp", "hidp", "wave"), all.x = TRUE)
  raw_data <- merge(raw_data, benefits_collapsed, by = c("hidp", "wave"), all.x = TRUE)

  
  fst::write_fst(raw_data, cache_path)

  return(raw_data)
}

clean_data <- function(DT) {
  library(data.table)

  if (!is.data.table(DT)) stop("Input must be a data.table")

  raw_data <- DT

  # Defining aliases
  age_responsible <- 18
  age_seek_employment <- 16
  age_leave_school <- 16
  age_form_partnership <- 18
  age_leave_parents <- 18
  age_retirement_min <- 50
  age_forced_retirement <- 75
  age_retirement <- 65
  age_reproductive_max <- 49
  age_max_dependent <- 17

  # drop if fre hhorig comes from ukhls iemb 2014-15 (hhorig = 8)
  raw_data <- raw_data[hhorig != 8]

  # probe for household identifier (needed for household-level variables)
  raw_data[, idhh := hidp]

  # probe for individual identifier (needed for individual-level variables)
  raw_data[, idind := pidp]

  # partner identifier probe (needed for household income)
  raw_data[, partnerid := ppid]

  # probe for mother and father id (mnspid = mother's pidp, fnspid = father's pidp)
  raw_data[, idmother := mnspid]
  raw_data[, idfather := fnspid]

  # interview date (needed for age calculation)
  raw_data[, int_date := lubridate::make_date(intdaty_dv, intdatm_dv, intdatd_dv)]

  # create age_probe variable for composite variables (e.g., employment status, partnership status, etc.)
  raw_data[, age_probe := fifelse(age_dv < 0, NA_integer_, as.integer(age_dv))]

  # create gender probe variable and recoding sex (2 -> 0 for female)
  raw_data[, gender_probe := fifelse(sex_dv == 2, 0, sex_dv)]

  # partner dummy variable (1 if has partner, 0 if not)
  raw_data[, has_partner := fifelse(partnerid > 0 & !is.na(partnerid), 1L, 0L)]

  # dependent children dummy variable (1 if has dependent children, 0 if not)
  raw_data[, depChild := fcase(
  age_probe >= age_max_dependent, 0L,                    # adults: never dependent
  is.na(age_probe) | !depchl_dv %in% c(1:2), NA_integer_,   # unknown status
  depchl_dv == 1 & (pns1pid > 0 | pns2pid > 0), 1L,
  default = 0L
  )]
  raw_data[, dnc := sum(depChild), by = .(wave, idhh)]   # NA if any member unknown

  ## drop temporary depChild
  raw_data[, depChild := NULL]

  # flag for being at or above state pension age (1 if at/above SPA, 0 if below)
  raw_data[, age_pension := fcase(
  # missing age/sex/interview year -> NA, not the default 0L ("below SPA")
  is.na(age_probe) | is.na(gender_probe) | intdaty_dv < 0, NA_integer_,
  # Men
  gender_probe == 1 & age_probe >= 66 & intdaty_dv >= 2020,             1L,
  gender_probe == 1 & age_probe >= 65 & intdaty_dv >= 2009 & intdaty_dv < 2020, 1L,
  # Women — phased SPA increase
  gender_probe == 0 & age_probe >= 66 & intdaty_dv >= 2021,              1L,
  gender_probe == 0 & age_probe >= 65 & intdaty_dv >= 2019 & intdaty_dv < 2021, 1L,
  gender_probe == 0 & age_probe >= 64 & intdaty_dv >= 2018 & intdaty_dv < 2019, 1L,
  gender_probe == 0 & age_probe >= 63 & intdaty_dv >= 2016 & intdaty_dv < 2018, 1L,
  gender_probe == 0 & age_probe >= 62 & intdaty_dv >= 2014 & intdaty_dv < 2016, 1L,
  gender_probe == 0 & age_probe >= 61 & intdaty_dv >= 2012 & intdaty_dv < 2014, 1L,
  gender_probe == 0 & age_probe >= 60 & intdaty_dv >= 2009 & intdaty_dv < 2012, 1L,

  default = 0L
)]
  
  # race into white vs non-white (1 if white, 0 if non-white)
  raw_data[, race := fcase(ethn_dv < 0, NA_integer_,
                           ethn_dv %in% c(1:4), 0L,
                           ethn_dv %in% c(5:97), 1L,
                           default = NA_integer_)]
  
  # job status recode: 1 if employed, 0 if unemployed or inactive
  raw_data[, les_c3 := fcase(
  jbstat < 0,                           NA_integer_,
  jbstat %in% c(1,2,5,12,13,14,15), 1L, # 1 = employed
  jbstat == 7,                       2L, # 2 = students
  jbstat %in% c(3,6,8,10,11,97,9,4), 3L, # 3 = unemployed or inactive
  default = NA_integer_
  )]

  ## people below age to leave home are not at risk of work, so set to not employed
  raw_data[age_probe < age_leave_parents, les_c3 := 3L]
  ## people under 16 is student
  raw_data[age_probe < age_seek_employment, les_c3 := 2L]

  # les_c4 is cloned but with a retired category
  raw_data[, les_c4 := les_c3]
  raw_data[jbstat == 4, les_c4 := 4L]

  # flag for adult children in the household (1 if has adult children, 0 if not)
  ## build parent tables
  mothers <- raw_data[gender_probe == 0, .(
    wave, hidp,
    idmother       = pidp,
    age_mother     = age_probe,
    pension_mother = age_pension,
    les_c4mother   = les_c4
  )]

  fathers <- raw_data[gender_probe == 1, .(
    wave, hidp,
    idfather       = pidp,
    age_father     = age_probe,
    pension_father = age_pension,
    les_c4father   = les_c4
  )]
  
  ## submerge parent info back into main data
  raw_data <- merge(raw_data, mothers, by = c("wave", "hidp", "idmother"), all.x = TRUE)
  raw_data <- merge(raw_data, fathers, by = c("wave", "hidp", "idfather"), all.x = TRUE)
  
  ## Phase 2 — set flag to 1
  raw_data[, adultchildflag := as.integer(
    (idmother > 0 | idfather > 0) & between(age_probe, age_max_dependent, age_forced_retirement) & partnerid <= 0
  )]

  ## Phase 3a — knock back to 0: parents retired (age_pension == 1 or les_c4 == 4)
  raw_data[pension_mother == 1  & is.na(pension_father),  adultchildflag := 0L]
  raw_data[is.na(pension_mother) & pension_father == 1,   adultchildflag := 0L]
  raw_data[pension_mother == 1  & pension_father == 1,    adultchildflag := 0L]

  raw_data[les_c4mother == 4L & is.na(les_c4father),     adultchildflag := 0L]
  raw_data[is.na(les_c4mother) & les_c4father == 4L,     adultchildflag := 0L]
  raw_data[les_c4mother == 4L & les_c4father == 4L,      adultchildflag := 0L]

  raw_data[les_c4mother == 4L & pension_father == 1,     adultchildflag := 0L]
  raw_data[les_c4father == 4L & pension_mother == 1,     adultchildflag := 0L]

  ## Phase 3b — knock back to 0: insufficient age gap (< 15 years)
  raw_data[(age_father - age_probe) <= 15 & is.na(age_mother),              adultchildflag := 0L]
  raw_data[is.na(age_father) & (age_mother - age_probe) <= 15,              adultchildflag := 0L]
  raw_data[(age_father - age_probe) <= 15 & (age_mother - age_probe) <= 15, adultchildflag := 0L]

  # household composition variables: number of adults (nadults) and number of children (nchildren)
  raw_data[, house_comp := fcase(
    has_partner == 1 & dnc == 0,                                                        1L,  # couple, no children
    has_partner == 1 & dnc > 0  & !is.na(dnc),                                         2L,  # couple, children
    has_partner == 0 & (dnc == 0 | age_probe <= age_responsible | adultchildflag == 1), 3L,  # single, no children
    has_partner == 0 & dnc > 0  & !is.na(dnc),                                         4L,  # single, children
    default = NA_integer_
  )]

  # CPI for inflation adjustment
  cpi_lookup <- c("2009" = 0.879, "2010" = 0.901, "2011" = 0.936, "2012" = 0.960,
                  "2013" = 0.982, "2014" = 0.996, "2015" = 1.000, "2016" = 1.010,
                  "2017" = 1.036, "2018" = 1.060, "2019" = 1.078, "2020" = 1.089,
                  "2021" = 1.116, "2022" = 1.205, "2023" = 1.286, "2024" = 1.329,
                  "2025" = 1.329)
  raw_data[, cpi := cpi_lookup[as.character(intdaty_dv)]]

  # ---- YPNB: Gross personal non-benefit income ----
  # recode sentinel missing values (-9, -1 -> NA) before summing
  sentinel_cols <- c("fimnlabgrs_dv", "fimnpen_dv", "inc_pp", "inc_tu", "inc_ma", "inc_fm", "inc_oth")
  raw_data[, (sentinel_cols) := map(.SD, \(x) fifelse(x < 0, NA_real_, x)), .SDcols = sentinel_cols]

  raw_data[, ypnb := rowSums(.SD, na.rm = TRUE), .SDcols = sentinel_cols]
  raw_data[, ypnb := fifelse(rowSums(!is.na(.SD)) == 0L, NA_real_, ypnb), .SDcols = sentinel_cols]
  raw_data[ypnb < 0, ypnb := 0]
  raw_data[, ypnb := ypnb / cpi]   # deflate to 2015 prices

  # ---- YPNBSP: Partner's gross personal non-benefit income ----
  # merge partner's ypnb into main dataset using partnerid and hidp as linking keys
  ypnb_lookup <- raw_data[, .(pidp, wave, hidp, ypnbsp = ypnb)]
  raw_data <- merge(raw_data, ypnb_lookup,
                    by.x = c("wave", "hidp", "partnerid"),
                    by.y = c("wave", "hidp", "pidp"),
                    all.x = TRUE)
  # ypnbsp is already in 2015 prices (copy of ypnb) - no second CPI division

  # ---- OECD modified equivalence scale (moecd_eq) ----
  raw_data[, depChild_013  := fcase(age_probe >= 0 & age_probe <= 13, 1L,
    is.na(age_probe) | !depchl_dv %in% c(1:2), NA_integer_,
    depchl_dv == 1 & (pns1pid > 0 | pns2pid > 0), 1L,
    default = 0L 
    )]
  raw_data[, depChild_1418 := fcase(age_probe >= 14 & age_probe <= 18, 1L,
    is.na(age_probe) | !depchl_dv %in% c(1:2), NA_integer_,
    depchl_dv == 1 & (pns1pid > 0 | pns2pid > 0), 1L,
    default = 0L)
  ]
  raw_data[, dnc013  := sum(depChild_013), by = .(wave, idhh)]
  raw_data[, dnc1418 := sum(depChild_1418), by = .(wave, idhh)]
  raw_data[, c("depChild_013", "depChild_1418") := NULL]

  raw_data[, moecd_eq := fcase(
    house_comp == 1L, 1.5,
    house_comp == 2L, 0.3 * dnc013 + 0.5 * dnc1418 + 1.5,
    house_comp == 3L, 1.0,
    house_comp == 4L, 0.3 * dnc013 + 0.5 * dnc1418 + 1.0,
    default = NA_real_
  )]
  raw_data[, c("dnc013", "dnc1418") := NULL]

  # ---- YHHNB: Equivalised household non-benefit income ----
  raw_data[, yhhnb := fcase(
    house_comp %in% c(1L, 2L), ypnb + fifelse(is.na(ypnbsp), 0, ypnbsp),
    house_comp %in% c(3L, 4L), ypnb,
    default = NA_real_
  )]
  raw_data[, yhhnb := yhhnb / moecd_eq]
  # ypnb is already in 2015 prices - no additional CPI division here
  # (do-file line 391 divides yhhnb by CPI again: double-deflation bug)
  raw_data[, yhhnb_asinh := asinh(yhhnb)]
  raw_data[, log_income  := yhhnb_asinh]

  # ---- HOME OWNERSHIP ----
  raw_data[, home_owner := fcase(
    tenure_dv < 0,                NA_integer_,
    tenure_dv >= 1 & tenure_dv <= 2, 1L,
    tenure_dv %in% c(3L, 4L, 5L, 6L, 7L),    0L,
    tenure_dv == 8L & hsownd == 97,              NA_integer_,
    default = NA_integer_
  )]

  # ---- UNEMPLOYMENT DUMMY ----
  raw_data[, unemp := as.integer(jbstat == 3)]
  raw_data[is.na(les_c3),                                       unemp := NA_integer_]
  raw_data[!is.na(age_probe) & age_probe < age_seek_employment, unemp := NA_integer_]
  raw_data[les_c4 == 4L & !is.na(unemp) & unemp == 1L,         unemp := 0L]

  # ---- ECONOMIC BENEFITS ----
  raw_data[, econ_benefits := fcase(
    fihhmnsben_dv > 0 & !is.na(fihhmnsben_dv), 1L,
    fihhmnsben_dv == 0,                         0L,
    default = NA_integer_
  )]
  raw_data[benefits_uc == 1L, econ_benefits := 1L]

  raw_data[, econ_benefits_nonuc := econ_benefits]
  raw_data[benefits_uc == 1L, econ_benefits_nonuc := 0L]

  raw_data[, econ_benefits_uc := econ_benefits]
  raw_data[benefits_uc == 0L, econ_benefits_uc := 0L]

  raw_data[, econ_benefits_lb := benefits_lb]
  raw_data[econ_benefits_uc == 1L, econ_benefits_lb := 0L]

  # ---- FINANCIAL DISTRESS ----
  raw_data[finnow < 0, finnow := NA_integer_]
  raw_data[, econ_dist := finnow]

  ## Main variable: economic distress dummy (1 if in financial distress, 0 if not)
  raw_data[, econ_dist_bin := fcase(
    econ_dist %in% 1:3, 0L,
    econ_dist %in% 4:5, 1L,
    default = NA_integer_
  )]

  raw_data[, econ_dist_bin_fact := as.factor(econ_dist_bin)]

  ####### Final modfications and recodes ########

  # generate long-term sick, retired and student dummies (integer 0/1)
  raw_data[, econ_ltsick  := fifelse(jbstat == 8,  1L, 0L)]
  raw_data[, econ_retire  := fifelse(les_c4 == 4L, 1L, 0L)]
  raw_data[, econ_student := fifelse(les_c4 == 2L, 1L, 0L)]

  # generating main exposure variable: employment status
  # long-term sick takes priority over les_c4 category
  raw_data[, econ_emp := fcase(
    econ_ltsick == 1L, "Long-term sick",
    les_c4 == 1L,      "Employed or self-employed",
    les_c4 == 2L,      "Student",
    les_c4 == 3L,      "Not employed (at risk of work)",
    les_c4 == 4L,      "Retired",
    default = NA_character_
  )]

  # binary employment: 0 = employed, 1 = not employed
  raw_data[, econ_emp_bin := fcase(
    econ_emp == "Employed or self-employed", 0L,
    !is.na(econ_emp),                        1L,
    default = NA_integer_
  )]

  #sensitivity exposure: binary, labour-force. Students and retired are not at risk of work.
  raw_data[, econ_emp_bin_narrow := fcase(
    econ_emp == "Employed or self-employed",                             0L,
    econ_emp %in% c("Not employed (at risk of work)", "Long-term sick"), 1L,
    default = NA_integer_
  )]

  #  convert les_c3 and les_c4 to factors (after all integer recodes are done)
  raw_data[, les_c3 := factor(les_c3, levels = 1:3,
    labels = c("Employed or self-employed", "Student", "Not employed"))]
  raw_data[, les_c4 := factor(les_c4, levels = 1:4,
    labels = c("Employed or self-employed", "Student", "Not employed", "Retired"))]

  # convert race into factor
  raw_data[, race := factor(race, levels = c(0L, 1L), labels = c("White", "Non-white"))]
  
  # convert house ownership to factor
  raw_data[, home_owner := factor(home_owner, levels = c(0L, 1L), labels = c("Renter", "Owner"))]

  # convert benefits receipt to factor (after econ_benefits_* derivations above)
  raw_data[, econ_benefits := factor(econ_benefits, levels = c(0L, 1L), labels = c("No benefits", "Benefits"))]

  # rename gender_probe to sex_dv (overwrite original sex_dv which has 1/2 coding and is less intuitive than 0/1)
  raw_data[, sex_dv := gender_probe]
  ## changing sex_dv == -9 to NA (missing) to match original coding
  raw_data[sex_dv == -9, sex_dv := NA_integer_]
  ## keeping only no missing values of sex
  raw_data <- raw_data[!is.na(sex_dv)]
  ## convert sex_dv to factor
  raw_data[, sex_dv := factor(sex_dv, levels = c(0, 1), labels = c("Female", "Male"))]

  # rename imputed age_probe to age_dv (overwrite original age_dv which has negative values for missing)
  raw_data[, age_dv := age_probe]

  # recoding gor_dv==-9 to NA and then forward/backward fill within person to impute missing region of residence (gor_dv)
  raw_data[gor_dv == -9, gor_dv := NA_integer_]
  setorder(raw_data, pidp, wave)
  raw_data[, gor_dv := nafill(gor_dv, type = "locf"), by = pidp]  # forward fill
  raw_data[, gor_dv := nafill(gor_dv, type = "nocb"), by = pidp]  # backward fill
  
  # transform gor_dv into factor
  gor_labels <- c("North East", "North West", "Yorkshire and the Humber", "East Midlands", "West Midlands",
                  "East of England", "London", "South East", "South West", "Wales", "Scotland", "Northern Ireland")
  raw_data[, gor_dv_fact := factor(gor_dv, levels = 1:12, labels = gor_labels)]

  # marital status (mastat_dv): recode to partnered (1) vs not partnered (0)
  # negatives and 0 "Child under 16" become NA (no adult marital status), then dropped.
  raw_data[mastat_dv < 0, mastat_dv := NA_integer_]
  ## partnered     = Married (2), Same-sex civil partnership (3), Living as couple (10)
  ## not partnered = Single/never married (1), Separated but legally married (4),
  ##                 Divorced (5), Widowed (6), Separated from civil partner (7),
  ##                 Former civil partner (8), Surviving civil partner (9)
  raw_data[, mastat_dv := fifelse(mastat_dv %in% c(2L, 3L, 10L), 1L,
                            fifelse(mastat_dv %in% c(1L, 4L, 5L, 6L, 7L, 8L, 9L), 0L, NA_integer_))]
  raw_data <- raw_data[!is.na(mastat_dv)]
  # convert to factor so mice's logreg method matches the variable type
  raw_data[, mastat_dv := factor(mastat_dv, levels = c(0L, 1L), labels = c("Not partnered", "Partnered"))]

  # scsf1, sf12mcs_dv and sf12pcs_dv recode missing from negative to NA
  raw_data[scsf1 < 0, scsf1 := NA_integer_]
  raw_data[sf12mcs_dv < 0, sf12mcs_dv := NA_integer_]
  raw_data[sf12pcs_dv < 0, sf12pcs_dv := NA_integer_]

  # recoding hiqual_dv missings to NA
  raw_data[hiqual_dv < 0, hiqual_dv := NA_integer_]

  # drop intermediate variables used for cleaning
  raw_data[, c("age_probe", "gender_probe", "has_partner", "age_pension",
               "idmother", "age_mother", "pension_mother", "les_c4mother",
               "idfather", "age_father", "pension_father", "les_c4father",
               "partnerid", "adultchildflag", "idhh", "idind") := NULL]
  
  # dropping other variables not needed
  raw_data[, c("jbstat", "unemp", "cpi", "econ_benefits_uc", "econ_benefits_lb", 
               "yhhnb_asinh", "yhhnb", "ypnbsp", "ypnb", "moecd_eq", "house_comp",
              "econ_benefits_nonuc", "hsownd", "tenure_dv", "depchl_dv", "nchild_dv", 
              "intdaty_dv", "intdatm_dv", "intdatd_dv", "sclfsato", "finnow", 
              "fimnpen_dv", "fimnlabgrs_dv", "scsf1", "inc_pp", "inc_tu", "inc_ma", 
              "inc_fm", "inc_oth", "benefits_uc", "benefits_lb", "fihhmnsben_dv", "ethn_dv",
              "econ_ltsick", "econ_student", "econ_retire") := NULL]
  return(raw_data)
}

preproc_data <- function(DT) {
  library(data.table)

  if (!is.data.table(DT)) stop("Input must be a data.table")

  # ---- 1. Filter to working-age and select analysis columns ----
  DT <- DT[age_dv >= 25 & age_dv <= 64]

  cols_to_keep <- c("pidp", "wave",
                    "sf12mcs_dv", "sf12pcs_dv", "log_income",
                    "econ_emp_bin", "econ_dist", "econ_dist_bin", "econ_dist_bin_fact", "econ_benefits", "gor_dv_fact",
                    "gor_dv", "mastat_dv", "home_owner", "dnc", "age_dv",
                    "race", "sex_dv", "hiqual_dv", "econ_emp", "econ_emp_bin_narrow")

  DT <- DT[, ..cols_to_keep]
  DT[, response := 1L]

  # ---- 2. Expand to balanced panel (all pidp x wave combinations) ----
  grid     <- CJ(pidp = unique(DT$pidp), wave = sort(unique(DT$wave)))
  exp_data <- merge(grid, DT, by = c("pidp", "wave"), all.x = TRUE)
  exp_data[is.na(response), response := 0L]
  setorder(exp_data, pidp, wave)

  # impute time-invariant variables for non-response rows (reponse = 0)
  exp_data[, sex_dv    := sex_dv[!is.na(sex_dv)][1L], by = pidp]          # time-invariant sex
  exp_data[, race      := race[!is.na(race)][1L],      by = pidp]          # time-invariant race
  exp_data[, hiqual_dv := hiqual_dv[!is.na(hiqual_dv)][1L], by = pidp]     # time-invariant highest qualification
  # collapse UKHLS hiqual_dv codes: 1-2 (degree/other higher degree) = High,
  # 3-5 (A-level/GCSE/other qualification) = Medium, 9 (no qualification) = Low
  exp_data[, hiqual_dv_fact := factor(fcase(
    hiqual_dv %in% 1:2, "High",
    hiqual_dv %in% 3:5, "Medium",
    hiqual_dv == 9L,    "Low"
  ), levels = c("High", "Medium", "Low"))]

  # ---- 3. Wave-1 baseline variable creation ----
#  base_cols <- c("pidp", "age_dv", "sex_dv", "gor_dv", "gor_dv_fact", "mastat_dv",
#                 "home_owner", "dnc", "hiqual_dv_fact", "race",
#                 "sf12mcs_dv", "sf12pcs_dv")

#  base_data <- exp_data[wave == 1L, ..base_cols]

  # rename all except pidp with _base suffix
#  setnames(base_data, setdiff(names(base_data), "pidp"),
#           paste0(setdiff(names(base_data), "pidp"), "_base"))

  # expand to all waves: each person gets 10 rows (person-first order)
#  n_persons <- nrow(base_data)
#  base_data  <- base_data[rep(seq_len(n_persons), each = 10L)]
#  base_data[, wave := rep(1:10, times = n_persons)]
#  base_data[, t0   := wave - 1L]

  # ---- 4. Join baseline with time-varying panel ----
#  pop_data <- merge(base_data, exp_data, by = c("pidp", "wave"), all.x = TRUE)

  # ---- 5. Recoding dnc to categorical variable ----
  pop_data <- exp_data[, dnc_fact := factor(fcase(
                         is.na(dnc), NA_character_,
                         dnc == 0L,  "Zero",
                         dnc == 1L,  "One",
                         default    = "Two or more"
                        ), levels = c("Zero", "One", "Two or more"))]
  
  #--- Transforming econ_emp_bin into factor variable ----
  pop_data[, econ_emp_bin_fact := as.factor(econ_emp_bin)]

  #--- Renaming sex_dv as a factor variable ----
  pop_data[, sex_dv_fact := sex_dv]

  # ---- 5. Final column selection ----
  final_cols <- c("pidp", "wave", "response", "race",
                  "sf12mcs_dv", "sf12pcs_dv", "log_income",
                  "econ_emp_bin", "econ_emp_bin_fact", "econ_dist", "econ_dist_bin", "econ_dist_bin_fact", "econ_benefits",
                  "gor_dv", "mastat_dv", "home_owner", "dnc", "dnc_fact", "age_dv",
                  "sex_dv_fact", "hiqual_dv_fact", "gor_dv_fact")

  pop_data[, ..final_cols]
} 