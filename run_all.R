library(mgcv)
library(mice)
library(splines)

for (file in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(file)
for (folder in c("private", file.path("results", c("tables", "figures", "audit")))) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

message("Reading data")
cohort <- build_cohort(read_source(settings$source_file))

message("Deriving phenotypes")
phenotypes <- derive_phenotypes(cohort)
cohort$phenotype <- phenotypes$assignments$phenotype[match(cohort$row_id, phenotypes$assignments$row_id)]
missing_clustering <- clustering_missingness(cohort)
flow <- participant_flow(cohort)
phenotype_data <- cohort[cohort$phenotype_population, ]
tte <- cohort[cohort$tte_population, ]
saveRDS(cohort, file.path("private", "cohort.rds"))

message("Multiple imputation")
phenotype_imputed <- impute_cohort(phenotype_data, settings$imputations, settings$seed)
tte_imputed <- impute_cohort(tte, settings$imputations, settings$seed + 1)

message("Phenotype description")
phenotype_results <- describe_phenotypes(phenotype_data, phenotype_imputed)

message("Positivity and calibration")
propensity <- fit_propensity_model(tte_imputed[[1]])
calibration_data <- model_data(tte_imputed[[1]])
calibration <- lapply(c(in_hospital_death = "in_hospital_death", recovery_365 = "recovery_365"), function(outcome) {
  data.frame(observed = calibration_data[[outcome]],
             predicted = fitted(fit_binary_model(calibration_data, outcome)))
})
states <- state_occupation(tte_imputed[[1]])

message("Target trial emulation: point estimates")
cluster <- start_cluster()
point <- point_estimates(cluster, tte_imputed)

message("Target trial emulation: bootstrap")
bootstrap <- run_bootstrap(cluster, tte)
parallel::stopCluster(cluster)

results <- add_intervals(point$estimates, bootstrap)
hte <- all_hte_tests(point$estimates, bootstrap)
evalues <- evalue_table(results, bootstrap)

message("Writing tables and figures")
output <- list(
  tte = tte, phenotype_data = phenotype_data, flow = flow, missing_clustering = missing_clustering,
  phenotypes = phenotypes, phenotype_results = phenotype_results, propensity = propensity,
  calibration = calibration, states = states, point = point,
  bootstrap = bootstrap, results = results, hte = hte, evalues = evalues
)
write_all_tables(output)
save_all_figures(output)

checks <- numerical_checks(results, hte, evalues, bootstrap)
write_audit(checks, bootstrap)
message("Done")
