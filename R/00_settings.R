settings <- list(
  source_file = file.path("..", "data.xlsx"),
  seed = 20260913,
  horizon = 365,
  imputations = 15,
  imputation_iterations = 10,
  bootstrap_replicates = 1000,
  workers = 12,
  figure_dpi = 600,
  bootstrap_folder = file.path("private", "bootstrap")
)

palette <- c(
  strategy_1 = "#0072B2",
  strategy_2 = "#D55E00",
  non_congestive = "#009E73",
  cardiorenal = "#E69F00",
  cardiometabolic = "#CC79A7",
  neutral = "#4D4D4D",
  light = "#D9E2E8"
)

strategy_labels <- c("mAFP-first", "V-A ECMO with LV unloading")
phenotype_labels <- c("Non-congestive", "Cardiorenal", "Cardiometabolic")
populations <- c("Overall", phenotype_labels)
phenotype_colors <- palette[c("non_congestive", "cardiorenal", "cardiometabolic")]
names(phenotype_colors) <- phenotype_labels

continuous_covariates <- c(
  "age", "bmi", "lvef", "lactate", "egfr", "alt", "platelets", "year"
)

binary_covariates <- c(
  "male", "ischemic_cmp", "dilated_cmp", "valvular_cmp", "atrial_fibrillation",
  "previous_hf", "diabetes", "ckd", "cirrhosis", "pad",
  "cardiac_arrest", "ohca", "shockable_rhythm", "amics", "intoxication", "myocarditis", "acute_on_chronic_hf",
  "acute_valvular", "septic_shock", "norepinephrine", "dobutamine", "rrt_shock",
  "electrical_storm", "postcardiotomy"
)

descriptive_continuous <- c(
  "age", "bmi", "lvef", "lactate", "egfr", "alt", "platelets", "low_flow", "year", "no_flow", "sofa"
)

descriptive_binary <- c(
  "male", "ischemic_cmp", "dilated_cmp", "valvular_cmp", "atrial_fibrillation",
  "previous_hf", "diabetes", "ckd", "cirrhosis", "pad", "copd", "dialysis",
  "cardiac_arrest", "ohca", "shockable_rhythm", "amics", "intoxication", "myocarditis", "acute_on_chronic_hf",
  "acute_valvular", "septic_shock", "norepinephrine", "dobutamine", "rrt_shock",
  "electrical_storm", "postcardiotomy"
)

etiology_variables <- c(
  "amics", "electrical_storm", "postcardiotomy", "intoxication", "myocarditis",
  "acute_on_chronic_hf", "acute_valvular", "septic_shock"
)

variable_labels <- c(
  age = "Age, years", male = "Male sex", bmi = "Body mass index, kg/m2",
  year = "Calendar year", ischemic_cmp = "Ischemic cardiomyopathy",
  dilated_cmp = "Dilated cardiomyopathy", valvular_cmp = "Valvular cardiomyopathy",
  atrial_fibrillation = "Atrial fibrillation", previous_hf = "Previous heart failure",
  copd = "Chronic obstructive pulmonary disease", diabetes = "Diabetes",
  ckd = "Chronic kidney disease", dialysis = "Chronic dialysis",
  cirrhosis = "Cirrhosis", pad = "Peripheral arterial disease",
  cardiac_arrest = "Cardiac arrest before support",
  ohca = "Out-of-hospital cardiac arrest", shockable_rhythm = "Shockable rhythm",
  no_flow = "No-flow time, minutes", low_flow = "Low-flow time, minutes",
  amics = "AMI-related cardiogenic shock", electrical_storm = "Electrical storm",
  postcardiotomy = "Postcardiotomy shock", intoxication = "Intoxication",
  myocarditis = "Fulminant myocarditis",
  acute_on_chronic_hf = "Acute-on-chronic heart failure",
  acute_valvular = "Acute valvular disease", septic_shock = "Septic shock",
  rare_etiology = "Shock etiologies below 10%, grouped",
  scai = "SCAI stage", sofa = "SOFA score",
  lvef = "Left ventricular ejection fraction, percent",
  norepinephrine = "Norepinephrine", dobutamine = "Dobutamine",
  epinephrine = "Epinephrine", ventilation = "Invasive ventilation",
  rrt_shock = "Renal replacement therapy at shock",
  lactate = "Lactate, mmol/L", egfr = "eGFR, mL/min/1.73 m2",
  alt = "Alanine aminotransferase, U/L", platelets = "Platelets, 10^3/mm3"
)
