summary_continuous <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.1f (%.1f to %.1f)", q[2], q[1], q[3])
}

summary_binary <- function(x) {
  sprintf("%d (%.1f%%)", sum(x == 1, na.rm = TRUE), 100 * mean(x == 1, na.rm = TRUE))
}

pool_rubin <- function(estimates, variances) {
  m <- length(estimates)
  within <- mean(variances)
  between <- var(estimates)
  total <- within + (1 + 1 / m) * between
  df <- (m - 1) * (1 + within / ((1 + 1 / m) * between))^2
  estimate <- mean(estimates)
  se <- sqrt(total)
  c(estimate = estimate, lower = estimate - qt(0.975, df) * se,
    upper = estimate + qt(0.975, df) * se,
    p_value = 2 * pt(abs(estimate / se), df, lower.tail = FALSE))
}

profile_variables <- c(
  "lactate", "egfr", "alt", "platelets", "age", "male", "bmi", "ischemic_cmp",
  "dilated_cmp", "valvular_cmp", "atrial_fibrillation", "previous_hf", "copd",
  "diabetes", "ckd", "dialysis", "cirrhosis", "pad", etiology_variables,
  "cardiac_arrest", "ohca", "shockable_rhythm", "no_flow", "low_flow", "sofa",
  "lvef", "norepinephrine", "dobutamine", "epinephrine", "ventilation", "rrt_shock"
)

phenotype_profile_table <- function(data) {
  groups <- c(list(Overall = rep(TRUE, nrow(data))),
              lapply(setNames(phenotype_labels, phenotype_labels), function(p) data$phenotype == p))
  rows <- lapply(profile_variables, function(v) {
    binary <- all(data[[v]] %in% c(0, 1, NA))
    summarize <- if (binary) summary_binary else summary_continuous
    values <- sapply(groups, function(keep) summarize(data[[v]][keep]))
    data.frame(characteristic = variable_labels[[v]], t(values), missing = sum(is.na(data[[v]])),
               check.names = FALSE)
  })
  scai_rows <- lapply(c("C", "D", "E"), function(stage) {
    values <- sapply(groups, function(keep) summary_binary(as.numeric(data$scai[keep] == stage)))
    data.frame(characteristic = paste("SCAI stage", stage), t(values), missing = 0, check.names = FALSE)
  })
  size <- data.frame(characteristic = "Patients, n", t(sapply(groups, sum)), missing = 0, check.names = FALSE)
  rbind(size, do.call(rbind, rows[1:4]), do.call(rbind, scai_rows), do.call(rbind, rows[-(1:4)]))
}

scai_by_phenotype <- function(data) {
  counts <- table(scai = data$scai, phenotype = factor(data$phenotype, levels = phenotype_labels))
  test <- chisq.test(counts)
  cells <- as.data.frame(counts, responseName = "n")
  cells$percent_within_phenotype <- 100 * cells$n / ave(cells$n, cells$phenotype, FUN = sum)
  cells$percent_within_scai <- 100 * cells$n / ave(cells$n, cells$scai, FUN = sum)
  list(
    cells = cells,
    test = data.frame(test = "Pearson chi-square test", statistic = unname(test$statistic),
                      df = unname(test$parameter), p_value = test$p.value,
                      cramer_v = sqrt(unname(test$statistic) / (sum(counts) * 2)))
  )
}

outcomes_by_phenotype_scai <- function(data, imputed) {
  groups <- expand.grid(phenotype = c("All", phenotype_labels), scai = c("All", "C", "D", "E"),
                        stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    keep <- (groups$phenotype[i] == "All" | data$phenotype == groups$phenotype[i]) &
      (groups$scai[i] == "All" | data$scai == groups$scai[i])
    recovery <- mean(sapply(imputed, function(d) mean(d$recovery_365[keep])))
    data.frame(groups[i, ], n = sum(keep), in_hospital_mortality = mean(data$in_hospital_death[keep]),
               recovery = recovery, row.names = NULL)
  }))
}

pooled_terms <- function(fits, terms, labels) {
  do.call(rbind, lapply(seq_along(terms), function(i) {
    pooled <- pool_rubin(sapply(fits, function(f) coef(f)[[terms[i]]]),
                         sapply(fits, function(f) vcov(f)[terms[i], terms[i]]))
    data.frame(contrast = labels[i], odds_ratio = exp(pooled[["estimate"]]),
               lower = exp(pooled[["lower"]]), upper = exp(pooled[["upper"]]),
               p_value = pooled[["p_value"]])
  }))
}

phenotype_prognostic_models <- function(imputed) {
  terms <- c("phenotypeCardiorenal", "phenotypeCardiometabolic", "scaiD", "scaiE")
  labels <- c("Cardiorenal versus Non-congestive", "Cardiometabolic versus Non-congestive",
              "SCAI stage D versus C", "SCAI stage E versus C")
  outcomes <- c(in_hospital_death = "In-hospital death", recovery_365 = "Recovery")
  do.call(rbind, lapply(names(outcomes), function(outcome) {
    fits <- lapply(imputed, function(d) {
      d$phenotype <- factor(d$phenotype, levels = phenotype_labels)
      d$center <- factor(d$center)
      d$support_group <- factor(d$support_group)
      gam(as.formula(paste(outcome, "~ phenotype + scai + support_group + age + male + year + s(center, bs = 're')")),
          data = d, family = binomial, method = "REML")
    })
    cbind(outcome = outcomes[[outcome]], pooled_terms(fits, terms, labels))
  }))
}

phenotype_mortality_models <- function(imputed) {
  terms <- c("phenotypeCardiorenal", "phenotypeCardiometabolic")
  labels <- c("Cardiorenal versus Non-congestive", "Cardiometabolic versus Non-congestive")
  formulas <- list(
    "Unadjusted" = in_hospital_death ~ phenotype,
    "Adjusted for SOFA score, age, and sex" = in_hospital_death ~ phenotype + sofa + age + male
  )
  do.call(rbind, lapply(names(formulas), function(model) {
    fits <- lapply(imputed, function(d) {
      d$phenotype <- factor(d$phenotype, levels = phenotype_labels)
      glm(formulas[[model]], data = d, family = binomial)
    })
    cbind(model = model, pooled_terms(fits, terms, labels))
  }))
}

phenotype_cumulative_incidence <- function(imputed) {
  runs <- lapply(imputed, function(d) {
    status <- competing_status(d)
    fit <- aalen_johansen_curves(status$time, status$status, d$phenotype)
    list(curves = fit$curves, p_value = fit$tests[1, "pv"])
  })
  curves <- do.call(rbind, lapply(runs, `[[`, "curves"))
  list(
    curves = aggregate(curves["cumulative_incidence"], curves[c("series", "day")], mean),
    gray_test = data.frame(test = "Gray test, in-hospital death by phenotype",
                           median_p_value = median(sapply(runs, `[[`, "p_value")),
                           imputations = length(imputed))
  )
}

describe_phenotypes <- function(data, imputed) {
  scai <- scai_by_phenotype(data)
  list(
    profile = phenotype_profile_table(data),
    scai_cells = scai$cells,
    scai_test = scai$test,
    outcomes = outcomes_by_phenotype_scai(data, imputed),
    prognostic_models = phenotype_prognostic_models(imputed),
    mortality_models = phenotype_mortality_models(imputed),
    cumulative_incidence = phenotype_cumulative_incidence(imputed)
  )
}
