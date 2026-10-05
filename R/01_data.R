read_source <- function(path) {
  raw <- as.data.frame(readxl::read_excel(path, skip = 3, .name_repair = "minimal",
                                          col_types = "text"))
  names(raw) <- trimws(names(raw))
  raw
}

clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[toupper(x) %in% c("", "NA", "NC", "ND", "NR", "N/A", "NULL", "-")] <- NA
  x
}

to_number <- function(x) {
  x <- gsub(",", ".", clean_text(x), fixed = TRUE)
  x[!grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", x)] <- NA
  as.numeric(x)
}

to_date <- function(x) {
  x <- clean_text(x)
  serial <- grepl("^[0-9]{5}(\\.[0-9]+)?$", x)
  out <- as.Date(rep(NA_character_, length(x)))
  out[serial] <- as.Date(floor(as.numeric(x[serial])), origin = "1899-12-30")
  for (format in c("%Y-%m-%d", "%m/%d/%Y")) {
    todo <- is.na(out) & !is.na(x) & !serial
    out[todo] <- as.Date(x[todo], format = format)
  }
  out
}

ckd_epi_2009 <- function(creatinine, age, male) {
  creatinine <- creatinine / 88.4
  k <- ifelse(male == 1, 0.9, 0.7)
  alpha <- ifelse(male == 1, -0.411, -0.329)
  141 * pmin(creatinine / k, 1)^alpha * pmax(creatinine / k, 1)^(-1.209) *
    0.993^age * ifelse(male == 1, 1, 1.018)
}

event_by_horizon <- function(time, bridge_event) {
  out <- ifelse(!is.na(time) & time >= 0, as.numeric(time <= settings$horizon), NA)
  out[is.na(out) & !bridge_event] <- 0
  out
}

build_cohort <- function(raw) {
  column <- function(name) raw[[name]]
  d <- data.frame(row_id = seq_len(nrow(raw)))

  d$center <- gsub("[[:space:]]+", " ", trimws(column("Center")))

  # Phenotype variables
  d$age <- to_number(column("Age"))
  d$male <- to_number(column("Sex (male = 1)"))
  d$creatinine <- to_number(column("Creatinin (micromol/l)"))
  d$lactate <- to_number(column("Lactate (mmol/L)"))
  d$alt <- to_number(column("ALAT (UI/L)"))
  d$platelets <- to_number(column("Platelets (number/mm3)")) / 1000
  d$egfr <- ckd_epi_2009(d$creatinine, d$age, d$male)
  d$adult <- !is.na(d$age) & d$age >= 18
  d$clustering_complete <- complete.cases(d[c("egfr", "lactate", "alt", "platelets")])
  d$phenotype_population <- d$adult & d$clustering_complete

  # Devices and treatment strategies
  device <- toupper(clean_text(column("IAPB/impella")))
  device[device %in% c("5", "5.00")] <- "5.0"
  device[device == "5.50"] <- "5.5"
  device[device == "2.50"] <- "2.5"
  d$device <- device
  d$ecmo <- to_number(column("ECMO"))
  d$ecmo_unloading <- to_number(column("ECMO Unloading"))
  d$device_start <- to_date(column("IABP/impella start date (mm/dd/yyyy)"))
  d$device_stop <- to_date(column("IABP/impella withdrawal date (mm/dd/yyyy)"))
  d$ecmo_start <- to_date(column("ECMO start date (mm/dd/yyyy)"))
  d$ecmo_stop <- to_date(column("ECMO withdrawal date (mm/dd/yyyy)"))
  device_days <- to_number(column("IABP/impella total length (days)"))
  ecmo_days <- to_number(column("ECMO total length (days)"))
  d$ecmo[is.na(d$ecmo) & d$ecmo_unloading %in% 1 & !is.na(d$ecmo_start)] <- 1

  start_text <- tolower(clean_text(column("IABP/impella start date (mm/dd/yyyy)")))
  text_before <- grepl("before", start_text)
  text_after <- grepl("after", start_text)
  text_more_24h <- grepl("more than 24h", start_text)

  has_ecmo <- d$ecmo %in% 1
  mafp <- device %in% c("CP", "2.5", "5.0", "5.5", "BOTH")
  device_minus_ecmo <- as.numeric(d$device_start - d$ecmo_start)
  before_ecmo <- has_ecmo & ((!is.na(device_minus_ecmo) & device_minus_ecmo < 0) |
                               (is.na(device_minus_ecmo) & text_before))
  same_day <- has_ecmo & !is.na(device_minus_ecmo) & device_minus_ecmo == 0

  d$strategy <- NA
  d$strategy[mafp & d$ecmo %in% 0] <- 1
  d$strategy[mafp & before_ecmo] <- 1
  d$strategy[has_ecmo & d$ecmo_unloading %in% 1 & !is.na(device) & !before_ecmo] <- 2
  d$strategy[mafp & same_day] <- 2
  d$iabp_before_ecmo <- device %in% "IABP" & before_ecmo
  d$treatment <- d$strategy - 1
  d$support_group <- ifelse(d$iabp_before_ecmo, 2, d$strategy)

  d$t0 <- d$device_start
  d$t0[has_ecmo] <- pmin(d$device_start[has_ecmo], d$ecmo_start[has_ecmo], na.rm = TRUE)
  admission <- to_date(column("admission date (mm/dd/yyyy)"))
  d$year <- as.numeric(format(d$t0, "%Y"))
  d$year[is.na(d$year)] <- as.numeric(format(admission[is.na(d$year)], "%Y"))

  rescue <- d$strategy %in% 1 & has_ecmo
  d$rescue_ecmo <- as.numeric(rescue)
  d$rescue_ecmo_day <- NA
  d$rescue_ecmo_day[rescue] <- as.numeric(d$ecmo_start[rescue] - d$t0[rescue])
  d$rescue_ecmo_day[rescue & (is.na(d$rescue_ecmo_day) | d$rescue_ecmo_day < 1)] <- 1

  strategy_2 <- d$strategy %in% 2
  d$unloading_day <- NA
  d$unloading_day[strategy_2] <- pmax(device_minus_ecmo[strategy_2], 0)
  d$unloading_day[strategy_2 & is.na(device_minus_ecmo) & text_after] <- 1
  d$unloading_day[strategy_2 & is.na(device_minus_ecmo) & text_after & text_more_24h] <- 2
  d$delayed_unloading <- strategy_2 & !is.na(d$unloading_day) & d$unloading_day >= 1
  d$unloading_device <- ifelse(device %in% "IABP", "IABP", ifelse(device %in% "BOTH", "Both", "mAFP"))

  d$tte_population <- d$phenotype_population & !d$iabp_before_ecmo & !is.na(d$strategy)

  # Baseline covariates
  d$bmi <- to_number(column("BMI (kg/m2)"))
  d$ischemic_cmp <- to_number(column("Ischemic CMP"))
  d$dilated_cmp <- to_number(column("Dilated CMP"))
  d$valvular_cmp <- to_number(column("Valvulopathy CMP"))
  d$atrial_fibrillation <- to_number(column("AF"))
  d$previous_hf <- to_number(column("Previous HF"))
  d$copd <- to_number(column("COPD"))
  d$diabetes <- to_number(column("Diabetes"))
  d$ckd <- to_number(column("CKD_GFR<60"))
  d$dialysis <- to_number(column("Chronic dialysis"))
  d$cirrhosis <- to_number(column("known cirrhosis"))
  d$pad <- to_number(column("Arteriopathy"))

  d$cardiac_arrest <- to_number(column("pre assist cardiac arrest"))
  arrest <- d$cardiac_arrest %in% 1
  no_arrest <- d$cardiac_arrest %in% 0
  ohca_text <- toupper(clean_text(column("OHCA")))
  ohca <- ifelse(ohca_text %in% c("EH", "1"), 1, ifelse(ohca_text %in% c("IH", "0"), 0, NA))
  among_arrest <- function(x) ifelse(arrest, x, ifelse(no_arrest, 0, NA))
  d$ohca <- among_arrest(ohca)
  d$shockable_rhythm <- among_arrest(to_number(column("Shockable rhythm")))
  d$no_flow <- among_arrest(to_number(column("NoFlow (min)")))
  d$low_flow <- among_arrest(to_number(column("lowflow (min)")))

  d$amics <- to_number(column("AMICS"))
  d$electrical_storm <- to_number(column("Electrical storm"))
  d$postcardiotomy <- to_number(column("Post cardiotomy"))
  d$intoxication <- to_number(column("Intoxication"))
  d$myocarditis <- to_number(column("Fulminant myocarditis"))
  d$acute_on_chronic_hf <- to_number(column("Acute on CHF"))
  d$acute_valvular <- to_number(column("Acute valvulopathy"))
  d$septic_shock <- to_number(column("Sepsis shock"))

  d$scai <- factor(toupper(clean_text(column("SCAI (lactate/ALAT/ACR)"))), levels = c("C", "D", "E"))
  d$sofa <- to_number(column("SOFA score"))
  d$lvef <- to_number(column("LVEF (%)"))
  d$norepinephrine <- to_number(column("Norepinephrine"))
  d$dobutamine <- to_number(column("Dobutamine"))
  d$epinephrine <- to_number(column("Epinephrine"))
  d$rrt_shock <- to_number(column("RRT_shock"))
  ventilation <- toupper(clean_text(column("Ventilation")))
  d$ventilation <- ifelse(ventilation == "MV", 1, ifelse(ventilation %in% c("VNI", "O2"), 0, NA))

  # In-hospital outcomes
  d$hospital_death <- to_number(column("Hospit"))
  d$d28 <- to_number(column("D28"))
  death_date <- to_date(column("Date of death (mm/dd/yyyy)"))
  discharge_date <- to_date(column("Date of discharge from hospital (mm/dd/yyyy)"))
  d$death_time <- as.numeric(death_date - d$t0)
  d$death_time[d$death_time < 0] <- NA
  discharge_time <- as.numeric(discharge_date - d$t0)
  discharge_time[discharge_time < 0] <- NA
  d$event_time <- ifelse(d$hospital_death == 1,
                         ifelse(is.na(d$death_time), discharge_time, d$death_time),
                         discharge_time)
  late_death <- d$hospital_death == 1 & !is.na(d$event_time) & d$event_time > settings$horizon
  d$in_hospital_death <- as.numeric(d$hospital_death == 1 & !late_death)

  # One-year recovery
  bridge <- clean_text(column("Bridge to (recovery, Transplant, LLVAD)"))
  lvad_time <- as.numeric(to_date(column("date_lvad")) - d$t0)
  transplant_time <- as.numeric(to_date(column("date_of_transplantation")) - d$t0)

  d$dead_365 <- ifelse(is.na(d$death_time), NA, as.numeric(d$death_time <= settings$horizon))
  no_time <- is.na(d$death_time)
  d$dead_365[no_time & (d$hospital_death == 1 | d$d28 %in% 1)] <- 1
  d$dead_365[no_time & is.na(death_date) & d$hospital_death == 0 & d$d28 %in% 0] <- 0

  lvad_365 <- event_by_horizon(lvad_time, bridge %in% c("LVAD", "LVAD+Transplant"))
  transplant_365 <- event_by_horizon(transplant_time, bridge %in% c("Transplant", "LVAD+Transplant"))
  d$advanced_365 <- ifelse(lvad_365 %in% 1 | transplant_365 %in% 1, 1,
                           ifelse(lvad_365 %in% 0 & transplant_365 %in% 0, 0, NA))
  lvad_time[lvad_time < 0] <- NA
  transplant_time[transplant_time < 0] <- NA
  d$advanced_time <- pmin(lvad_time, transplant_time, na.rm = TRUE)

  device_end <- as.numeric(d$device_stop - d$t0)
  device_offset <- as.numeric(d$device_start - d$t0)
  device_end[is.na(device_end)] <- (device_offset + device_days)[is.na(device_end)]
  ecmo_end <- as.numeric(d$ecmo_stop - d$t0)
  ecmo_offset <- as.numeric(d$ecmo_start - d$t0)
  ecmo_end[is.na(ecmo_end)] <- (ecmo_offset + ecmo_days)[is.na(ecmo_end)]
  d$support_end <- ifelse(has_ecmo, pmax(device_end, ecmo_end), device_end)
  free_365 <- ifelse(is.na(d$support_end), NA, as.numeric(d$support_end <= settings$horizon))

  d$recovery_365 <- NA
  d$recovery_365[d$dead_365 %in% 1 | d$advanced_365 %in% 1 | free_365 %in% 0] <- 0
  d$recovery_365[d$dead_365 %in% 0 & d$advanced_365 %in% 0 & free_365 %in% 1] <- 1

  d
}

participant_flow <- function(cohort) {
  phenotype <- cohort$phenotype_population
  data.frame(
    stage = c(
      "Records assessed for phenotyping", "Excluded: age below 18 years",
      "Adults eligible for phenotyping", "Excluded: missing clustering variables",
      "Phenotype population", "Excluded: IABP initiated before V-A ECMO",
      "Excluded: device sequence not classifiable", "TTE population",
      "mAFP-first", "V-A ECMO with LV unloading"
    ),
    n = c(
      nrow(cohort), sum(!cohort$adult), sum(cohort$adult),
      sum(cohort$adult & !cohort$clustering_complete), sum(phenotype),
      sum(phenotype & cohort$iabp_before_ecmo),
      sum(phenotype & !cohort$iabp_before_ecmo & is.na(cohort$strategy)),
      sum(cohort$tte_population),
      sum(cohort$tte_population & cohort$strategy == 1),
      sum(cohort$tte_population & cohort$strategy == 2)
    )
  )
}
