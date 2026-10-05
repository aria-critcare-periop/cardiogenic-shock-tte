imputation_variables <- c(
  "center", "support_group", "phenotype", "hospital_death", "event_time",
  "recovery_365", descriptive_continuous, descriptive_binary, "scai",
  "ventilation"
)

impute_cohort <- function(data, m, seed) {
  work <- data[imputation_variables]
  categorical <- c("center", "support_group", "phenotype", "hospital_death", "recovery_365",
                   descriptive_binary, "ventilation")
  for (v in categorical) work[[v]] <- factor(work[[v]])

  incomplete <- names(work)[colSums(is.na(work)) > 0]
  method <- mice::make.method(work)
  method[c("hospital_death", "event_time")] <- ""
  imputation <- mice::mice(work, m = m, maxit = settings$imputation_iterations,
                           method = method, seed = seed, printFlag = FALSE)

  lapply(seq_len(m), function(i) {
    completed <- mice::complete(imputation, i)
    out <- data
    for (v in incomplete) {
      out[[v]] <- if (v == "scai") completed[[v]] else as.numeric(as.character(completed[[v]]))
    }
    out
  })
}

missingness_table <- function(data) {
  variables <- c(descriptive_continuous, descriptive_binary, "scai", "event_time", "recovery_365")
  labels <- c(variable_labels, event_time = "Time to in-hospital death or discharge",
              recovery_365 = "Recovery at one year")
  data.frame(
    variable = variables,
    label = unname(labels[variables]),
    n_missing = sapply(data[variables], function(x) sum(is.na(x))),
    percent_missing = sapply(data[variables], function(x) 100 * mean(is.na(x)))
  )
}
