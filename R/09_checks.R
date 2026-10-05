numerical_checks <- function(results, hte, evalues, bootstrap) {
  by_metric <- function(metric) {
    rows <- results$metric == metric
    setNames(results$value[rows], paste(results$analysis, results$outcome, results$population)[rows])
  }
  keys <- names(by_metric("risk_1"))
  risks <- results[results$metric %in% c("risk_1", "risk_2"), ]
  rmtl <- results[results$metric %in% c("rmtl_1", "rmtl_2"), ]
  data.frame(
    check = c(
      "Standardized risks between 0 and 1",
      "Lower confidence limits not above upper limits",
      "Risk differences equal risk 2 minus risk 1",
      "Risk ratios equal risk 2 divided by risk 1",
      "RMTL between 0 and 365 days",
      "Heterogeneity p-values between 0 and 1",
      "E-values at least 1",
      "All bootstrap replicates completed"
    ),
    passed = c(
      all(risks$value >= 0 & risks$value <= 1),
      all(results$lower <= results$upper),
      all(abs(by_metric("difference")[keys] - (by_metric("risk_2")[keys] - by_metric("risk_1")[keys])) < 1e-8),
      all(abs(by_metric("ratio")[keys] - by_metric("risk_2")[keys] / by_metric("risk_1")[keys]) < 1e-8),
      all(rmtl$value >= 0 & rmtl$value <= settings$horizon),
      all(hte$p_value >= 0 & hte$p_value <= 1),
      all(c(evalues$e_value_estimate, evalues$e_value_confidence_limit) >= 1),
      length(unique(bootstrap$replicate)) == settings$bootstrap_replicates
    )
  )
}

write_audit <- function(checks, bootstrap) {
  folder <- file.path("results", "audit")
  write.csv(checks, file.path(folder, "numerical_checks.csv"), row.names = FALSE)

  code_files <- c("run_all.R", list.files("R", full.names = TRUE))
  write.csv(data.frame(file = code_files, md5 = unname(tools::md5sum(code_files))),
            file.path(folder, "code_checksums.csv"), row.names = FALSE)

  table_files <- list.files(file.path("results", "tables"), pattern = "\\.csv$", full.names = TRUE)
  patient_columns <- c("row_id", "patient", "t0", "death_time", "event_time")
  aggregated <- sapply(table_files, function(f) !any(names(read.csv(f, check.names = FALSE)) %in% patient_columns))
  write.csv(data.frame(table = basename(table_files), aggregated_only = unname(aggregated)),
            file.path(folder, "confidentiality_checks.csv"), row.names = FALSE)

  metadata <- data.frame(
    item = c("run_time", "r_version", "source_file", "source_md5", "seed", "imputations",
             "imputation_iterations", "bootstrap_replicates", "bootstrap_completed", "workers"),
    value = c(format(Sys.time()), R.version.string, basename(settings$source_file),
              unname(tools::md5sum(settings$source_file)), settings$seed, settings$imputations,
              settings$imputation_iterations, settings$bootstrap_replicates,
              length(unique(bootstrap$replicate)), settings$workers)
  )
  write.csv(metadata, file.path(folder, "run_metadata.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(folder, "session_info.txt"))

  stopifnot(all(checks$passed), all(aggregated))
}
