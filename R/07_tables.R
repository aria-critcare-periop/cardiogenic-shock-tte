write_table <- function(x, name) {
  write.csv(x, file.path("results", "tables", paste0(name, ".csv")), row.names = FALSE, na = "")
}

write_html_table <- function(x, name, title) {
  cell <- function(v) {
    if (is.numeric(v)) v <- ifelse(is.na(v), "", formatC(v, digits = 3, format = "f"))
    gsub("<", "&lt;", as.character(v), fixed = TRUE)
  }
  header <- paste0("<th>", gsub("_", " ", names(x)), "</th>", collapse = "")
  body <- sapply(seq_len(nrow(x)), function(i) {
    paste0("<tr>", paste0("<td>", sapply(x[i, ], cell), "</td>", collapse = ""), "</tr>")
  })
  html <- c(
    "<!doctype html><html><head><meta charset='utf-8'><style>",
    "body{font-family:Arial,sans-serif;color:#222;margin:32px}h1{font-size:18px}",
    "table{border-collapse:collapse;width:100%;font-size:11px}",
    "th{text-align:left;border-top:1.5px solid #222;border-bottom:1px solid #666;padding:7px}",
    "td{padding:6px;border-bottom:0.5px solid #d8d8d8}tr:last-child td{border-bottom:1.5px solid #222}",
    "</style></head><body>",
    paste0("<h1>", title, "</h1><table><thead><tr>", header, "</tr></thead><tbody>"),
    body, "</tbody></table></body></html>"
  )
  writeLines(html, file.path("results", "tables", paste0(name, ".html")), useBytes = TRUE)
}

standardized_mean_difference <- function(x, treatment) {
  a <- x[treatment == 0]
  b <- x[treatment == 1]
  (mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
}

baseline_by_strategy <- function(data) {
  row <- function(label, x, binary) {
    summarize <- if (binary) summary_binary else summary_continuous
    data.frame(characteristic = label,
               mafp_first = summarize(x[data$treatment == 0]),
               ecmo_lv_unloading = summarize(x[data$treatment == 1]),
               missing = sum(is.na(x)),
               standardized_mean_difference = standardized_mean_difference(x, data$treatment))
  }
  rows <- list(data.frame(characteristic = "Patients, n", mafp_first = as.character(sum(data$treatment == 0)),
                          ecmo_lv_unloading = as.character(sum(data$treatment == 1)),
                          missing = 0, standardized_mean_difference = NA))
  for (p in phenotype_labels) rows[[length(rows) + 1]] <- row(p, as.numeric(data$phenotype == p), TRUE)
  for (stage in c("C", "D", "E")) {
    rows[[length(rows) + 1]] <- row(paste("SCAI stage", stage), as.numeric(data$scai == stage), TRUE)
  }
  for (v in descriptive_continuous) rows[[length(rows) + 1]] <- row(variable_labels[[v]], data[[v]], FALSE)
  for (v in descriptive_binary) rows[[length(rows) + 1]] <- row(variable_labels[[v]], data[[v]], TRUE)
  do.call(rbind, rows)
}

effect_table <- function(results, hte, analysis, outcome, metrics) {
  rows <- results[results$analysis == analysis & results$outcome == outcome & results$metric %in% metrics, ]
  wide <- reshape(rows[c("population", "metric", "value", "lower", "upper")],
                  idvar = "population", timevar = "metric", direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  names(wide) <- sub("^(lower|upper)\\.(.*)$", "\\2_\\1", names(wide))
  wide <- wide[match(populations, wide$population), c("population", as.vector(sapply(metrics, function(m) paste0(m, c("", "_lower", "_upper")))))]
  tests <- hte[hte$analysis == analysis & hte$outcome == outcome, ]
  for (i in seq_len(nrow(tests))) wide[[paste0("hte_p_value_", tests$metric[i])]] <- tests$p_value[i]
  cbind(analysis = analysis, outcome = outcome, wide, row.names = NULL)
}

sensitivity_table <- function(results, hte) {
  do.call(rbind, lapply(sensitivity_analyses, function(a) {
    rbind(effect_table(results, hte, a, "In-hospital death", c("risk_1", "risk_2", "difference", "ratio")),
          effect_table(results, hte, a, "Recovery", c("risk_1", "risk_2", "difference", "ratio")))
  }))
}

all_hte_tests <- function(point, bootstrap) {
  specs <- unique(point[point$metric %in% c("difference", "rmtl_difference"), c("analysis", "outcome", "metric")])
  do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
    hte_test(point, bootstrap, specs$analysis[i], specs$outcome[i], specs$metric[i])
  }))
}

evalue_table <- function(results, bootstrap) {
  outcomes <- c("In-hospital death", "Recovery")
  ratios <- results[results$analysis == "Primary" & results$outcome %in% outcomes & results$metric == "ratio", ]
  effects <- data.frame(outcome = ratios$outcome, contrast = paste(ratios$population, "risk ratio"),
                        estimate = ratios$value, lower = ratios$lower, upper = ratios$upper)
  contrasts <- do.call(rbind, lapply(outcomes, function(outcome) {
    point <- ratios[ratios$outcome == outcome, ]
    rr <- setNames(point$value, point$population)
    boot <- bootstrap[bootstrap$analysis == "Primary" & bootstrap$outcome == outcome & bootstrap$metric == "ratio", ]
    boot_rr <- tapply(boot$value, list(boot$replicate, boot$population), mean)
    do.call(rbind, lapply(c("Cardiorenal", "Cardiometabolic"), function(p) {
      interval <- percentile_interval(boot_rr[, p] / boot_rr[, "Non-congestive"])
      data.frame(outcome = outcome, contrast = paste(p, "versus Non-congestive, ratio of risk ratios"),
                 estimate = rr[[p]] / rr[["Non-congestive"]], lower = interval[1], upper = interval[2])
    }))
  }))
  table <- rbind(effects, contrasts)
  table$e_value_estimate <- e_value(table$estimate)
  table$e_value_confidence_limit <- e_value_limit(table$lower, table$upper)
  table
}

propensity_summary <- function(propensity) {
  data <- propensity$data
  groups <- expand.grid(population = populations, strategy = strategy_labels, stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    keep <- in_population(data, groups$population[i]) &
      data$treatment == match(groups$strategy[i], strategy_labels) - 1
    q <- quantile(propensity$propensity[keep], c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1))
    data.frame(groups[i, ], n = sum(keep), minimum = q[1], p01 = q[2], q1 = q[3], median = q[4],
               q3 = q[5], p99 = q[6], maximum = q[7], row.names = NULL)
  }))
}

state_names <- c("Alive, support-free", "Ongoing temporary MCS", "LVAD or transplant", "Death")

state_occupation <- function(data) {
  death_time <- ifelse(is.na(data$death_time) & data$in_hospital_death == 1, data$event_time, data$death_time)
  do.call(rbind, lapply(0:1, function(treatment) {
    keep <- data$treatment == treatment
    d <- data[keep, ]
    do.call(rbind, lapply(seq(0, settings$horizon, by = 5), function(day) {
      state <- ifelse(!is.na(death_time[keep]) & death_time[keep] <= day, "Death",
               ifelse(!is.na(d$advanced_time) & d$advanced_time <= day, "LVAD or transplant",
               ifelse((!is.na(d$support_end) & d$support_end <= day) |
                        (day == settings$horizon & d$recovery_365 == 1),
                      "Alive, support-free", "Ongoing temporary MCS")))
      data.frame(strategy = strategy_labels[treatment + 1], day = day, state = state_names,
                 probability = sapply(state_names, function(s) mean(state == s)), row.names = NULL)
    }))
  }))
}

write_all_tables <- function(output) {
  with(output, {
    table_1 <- baseline_by_strategy(tte)
    write_table(table_1, "table_1_baseline_by_strategy")
    write_html_table(table_1, "table_1_baseline_by_strategy", "Table 1. Baseline characteristics by treatment strategy")

    binary_metrics <- c("risk_1", "risk_2", "difference", "ratio")
    tables <- list(
      table_2_in_hospital_mortality = effect_table(results, hte, "Primary", "In-hospital death", binary_metrics),
      table_3_cumulative_incidence_rmtl = effect_table(results, hte, "Primary", "In-hospital death, competing risks",
                                                       c(binary_metrics, "rmtl_1", "rmtl_2", "rmtl_difference")),
      table_4_recovery = effect_table(results, hte, "Primary", "Recovery", binary_metrics)
    )
    titles <- c("Table 2. In-hospital mortality by day 365", "Table 3. Cumulative incidence of in-hospital death and RMTL",
                "Table 4. Recovery at one year")
    for (i in seq_along(tables)) {
      write_table(tables[[i]], names(tables)[i])
      write_html_table(tables[[i]], names(tables)[i], titles[i])
    }

    by_phenotype <- do.call(rbind, lapply(phenotype_labels, function(p) {
      cbind(phenotype = p, baseline_by_strategy(tte[tte$phenotype == p, ]))
    }))
    write_table(by_phenotype, "table_s27_baseline_by_strategy_within_phenotype")

    write_table(flow, "table_s1_participant_flow")
    write_table(missing_clustering$by_variable, "table_s2_clustering_missingness")
    write_table(missing_clustering$pattern, "table_s3_clustering_missing_patterns")
    write_table(merge(phenotypes$quality, data.frame(k = 3, selected = TRUE), all.x = TRUE), "table_s4_cluster_quality")
    write_table(phenotypes$transformation, "table_s5_clustering_transformations")
    write_table(phenotype_results$profile, "table_s6_phenotype_profiles")
    write_table(phenotype_results$scai_cells, "table_s7_scai_by_phenotype")
    write_table(phenotype_results$scai_test, "table_s8_scai_phenotype_test")
    write_table(phenotype_results$outcomes, "table_s9_outcomes_by_phenotype_scai")
    write_table(phenotype_results$prognostic_models, "table_s10_phenotype_prognostic_models")
    write_table(phenotype_results$mortality_models, "table_s11_phenotype_mortality_models")
    write_table(phenotype_results$cumulative_incidence$gray_test, "table_s12_phenotype_gray_test")
    write_table(missingness_table(tte), "table_s13_covariate_missingness")
    write_table(as.data.frame(table(center = tte$center, strategy = strategy_labels[tte$treatment + 1]), responseName = "n"),
                "table_s14_strategy_by_center")
    write_table(as.data.frame(table(year = tte$year, strategy = strategy_labels[tte$treatment + 1], useNA = "ifany"),
                              responseName = "n"), "table_s15_strategy_by_year")
    write_table(as.data.frame(table(phenotype = tte$phenotype, strategy = strategy_labels[tte$treatment + 1]),
                              responseName = "n"), "table_s16_strategy_by_phenotype")
    write_table(propensity_summary(propensity), "table_s17_propensity_distribution")
    write_table(hte, "table_s18_heterogeneity_tests")
    write_table(sensitivity_table(results, hte), "table_s19_clinical_sensitivity_analyses")
    write_table(evalues, "table_s20_e_values")
    write_table(point$diagnostics, "table_s21_outcome_model_diagnostics")
    write_table(point$curves, "table_s22_standardized_cumulative_incidence")
    write_table(point$aalen_johansen, "table_s23_unadjusted_cumulative_incidence")
    write_table(point$rescue, "table_s24_rescue_ecmo_cumulative_incidence")
    write_table(states, "table_s25_state_occupation")
    write_table(data.frame(requested = settings$bootstrap_replicates,
                           completed = length(unique(bootstrap$replicate))), "table_s26_bootstrap_completion")
  })
}
