adjustment_terms <- c(continuous_covariates, binary_covariates, "scai")
adjustment_formula <- paste(adjustment_terms, collapse = " + ")

model_data <- function(data) {
  data$phenotype <- factor(data$phenotype, levels = phenotype_labels)
  data$center <- factor(data$center)
  for (v in continuous_covariates) data[[v]] <- (data[[v]] - mean(data[[v]])) / sd(data[[v]])
  data
}

in_population <- function(data, population) {
  if (population == "Overall") rep(TRUE, nrow(data)) else data$phenotype == population
}

# Binary outcomes

fit_binary_model <- function(data, outcome) {
  formula <- as.formula(paste(outcome, "~ treatment * phenotype +", adjustment_formula,
                              "+ s(center, bs = 're')"))
  gam(formula, data = data, family = binomial, method = "REML")
}

standardized_risks <- function(model, data) {
  risk_1 <- predict(model, newdata = transform(data, treatment = 0), type = "response")
  risk_2 <- predict(model, newdata = transform(data, treatment = 1), type = "response")
  do.call(rbind, lapply(populations, function(population) {
    keep <- in_population(data, population)
    r1 <- mean(risk_1[keep])
    r2 <- mean(risk_2[keep])
    data.frame(population = population,
               metric = c("risk_1", "risk_2", "difference", "ratio"),
               value = c(r1, r2, r2 - r1, r2 / r1))
  }))
}

binary_estimates <- function(data, analysis) {
  death <- fit_binary_model(data, "in_hospital_death")
  recovery <- fit_binary_model(data, "recovery_365")
  rbind(
    cbind(analysis = analysis, outcome = "In-hospital death", standardized_risks(death, data)),
    cbind(analysis = analysis, outcome = "Recovery", standardized_risks(recovery, data))
  )
}

# Competing risks

patient_days <- function(data) {
  last_day <- pmin(data$event_time, settings$horizon)
  rows <- rep(seq_len(nrow(data)), last_day + 1)
  days <- data[rows, ]
  days$day <- sequence(last_day + 1) - 1
  event_day <- days$day == days$event_time
  days$death <- as.numeric(event_day & days$hospital_death == 1)
  days$discharge <- as.numeric(event_day & days$hospital_death == 0)
  days
}

fit_hazard_model <- function(days, event) {
  formula <- as.formula(paste(event, "~ ns(day, df = 4) + treatment * phenotype +",
                              adjustment_formula, "+ s(center, bs = 're')"))
  bam(formula, data = days, family = quasibinomial(), method = "fREML")
}

daily_hazards <- function(model, newdata, days) {
  patient_part <- predict(model, newdata = newdata)
  reference <- newdata[rep(1, length(days)), ]
  reference$day <- days
  time_part <- predict(model, newdata = reference)
  plogis(outer(patient_part, time_part - time_part[1], "+"))
}

cumulative_incidence <- function(death_model, discharge_model, data) {
  days <- 0:settings$horizon
  curves <- list()
  estimates <- list()
  cif <- list()
  for (treatment in 0:1) {
    newdata <- data
    newdata$treatment <- treatment
    newdata$day <- 0
    hazard_death <- daily_hazards(death_model, newdata, days)
    hazard_discharge <- daily_hazards(discharge_model, newdata, days)
    remain <- pmax(1 - hazard_death - hazard_discharge, 0)
    in_hospital <- cbind(1, t(apply(remain, 1, cumprod))[, -length(days)])
    cif_death <- t(apply(in_hospital * hazard_death, 1, cumsum))
    cif_discharge <- t(apply(in_hospital * hazard_discharge, 1, cumsum))
    for (population in populations) {
      keep <- in_population(data, population)
      death_curve <- unname(colMeans(cif_death[keep, , drop = FALSE]))
      curves[[length(curves) + 1]] <- data.frame(
        population = population, treatment = treatment, day = days,
        death = death_curve, discharge = unname(colMeans(cif_discharge[keep, , drop = FALSE]))
      )
      cif[[paste(population, treatment)]] <- c(
        risk = death_curve[length(days)],
        rmtl = sum(death_curve[days < settings$horizon])
      )
    }
  }
  for (population in populations) {
    a <- cif[[paste(population, 0)]]
    b <- cif[[paste(population, 1)]]
    estimates[[population]] <- data.frame(
      population = population,
      metric = c("risk_1", "risk_2", "difference", "ratio", "rmtl_1", "rmtl_2", "rmtl_difference"),
      value = c(a["risk"], b["risk"], b["risk"] - a["risk"], b["risk"] / a["risk"],
                a["rmtl"], b["rmtl"], b["rmtl"] - a["rmtl"])
    )
  }
  list(estimates = do.call(rbind, estimates), curves = do.call(rbind, curves))
}

# Unadjusted cumulative incidence

aalen_johansen_curves <- function(time, status, group, days = 0:settings$horizon) {
  fit <- cmprsk::cuminc(time, status, group, cencode = 0)
  grid <- cmprsk::timepoints(fit, days)$est
  for (i in seq_len(nrow(grid))) {
    for (j in seq_along(days)) {
      if (is.na(grid[i, j])) grid[i, j] <- if (j == 1) 0 else grid[i, j - 1]
    }
  }
  curves <- data.frame(
    series = rep(rownames(grid), each = length(days)),
    day = rep(days, nrow(grid)),
    cumulative_incidence = as.vector(t(grid))
  )
  list(curves = curves, tests = fit$Tests)
}

competing_status <- function(data) {
  status <- ifelse(data$event_time > settings$horizon, 0, ifelse(data$hospital_death == 1, 1, 2))
  list(time = pmin(data$event_time, settings$horizon), status = status)
}

# Main estimation for one dataset

sensitivity_analyses <- c(
  "AMICS", "2017 or later", "V-A ECMO with mAFP versus mAFP",
  "V-A ECMO with IABP versus mAFP", "Strategy 2 excluding delayed LV unloading"
)

sensitivity_subsets <- function(data) {
  setNames(list(
    data$amics == 1,
    data$year >= 2017,
    data$treatment == 0 | data$unloading_device == "mAFP",
    data$treatment == 0 | data$unloading_device == "IABP",
    !data$delayed_unloading
  ), sensitivity_analyses)
}

estimate_all <- function(imputed, details = FALSE) {
  data <- model_data(imputed)

  estimates <- list(binary_estimates(data, "Primary"))

  days <- patient_days(data)
  primary <- cumulative_incidence(fit_hazard_model(days, "death"),
                                  fit_hazard_model(days, "discharge"), data)
  estimates[[2]] <- cbind(analysis = "Primary", outcome = "In-hospital death, competing risks",
                          primary$estimates)

  subsets <- sensitivity_subsets(imputed)
  for (name in names(subsets)) {
    estimates[[length(estimates) + 1]] <- binary_estimates(model_data(imputed[subsets[[name]], ]), name)
  }

  result <- list(estimates = do.call(rbind, estimates))
  if (details) {
    result$curves <- cbind(analysis = "Primary", primary$curves)
    result$diagnostics <- rbind(
      binary_model_diagnostics(fit_binary_model(data, "in_hospital_death"), data, "in_hospital_death"),
      binary_model_diagnostics(fit_binary_model(data, "recovery_365"), data, "recovery_365")
    )
    competing <- competing_status(imputed)
    group <- paste(strategy_labels[imputed$treatment + 1], imputed$phenotype, sep = " | ")
    result$aalen_johansen <- rbind(
      cbind(scope = "Overall", aalen_johansen_curves(competing$time, competing$status,
                                                     strategy_labels[imputed$treatment + 1])$curves),
      cbind(scope = "Phenotype", aalen_johansen_curves(competing$time, competing$status, group)$curves)
    )
    s1 <- imputed$treatment == 0
    rescue_time <- ifelse(imputed$rescue_ecmo == 1 & imputed$rescue_ecmo_day <= imputed$event_time,
                          imputed$rescue_ecmo_day, competing$time)
    rescue_status <- ifelse(imputed$rescue_ecmo == 1 & imputed$rescue_ecmo_day <= imputed$event_time,
                            1, competing$status + (competing$status > 0))
    result$rescue <- aalen_johansen_curves(rescue_time[s1], rescue_status[s1],
                                           imputed$phenotype[s1])$curves
    overall_rescue <- aalen_johansen_curves(rescue_time[s1], rescue_status[s1], rep("Overall", sum(s1)))$curves
    result$rescue <- rbind(overall_rescue, result$rescue)
  }
  result
}

# Propensity

fit_propensity_model <- function(imputed) {
  data <- model_data(imputed)
  etiologies <- intersect(etiology_variables, binary_covariates)
  rare <- etiologies[colMeans(data[etiologies]) < 0.10]
  data$rare_etiology <- as.numeric(rowSums(data[rare]) > 0)
  terms <- c(setdiff(adjustment_terms, rare), "rare_etiology")
  model <- glm(as.formula(paste("treatment ~ phenotype +", paste(terms, collapse = " + "))),
               family = binomial, data = data)
  list(model = model, propensity = fitted(model), rare = rare, data = data)
}

# Inference

c_statistic <- function(observed, predicted) {
  ranks <- rank(predicted)
  n1 <- sum(observed == 1)
  n0 <- sum(observed == 0)
  (sum(ranks[observed == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

binary_model_diagnostics <- function(model, data, outcome) {
  observed <- data[[outcome]]
  predicted <- fitted(model)
  linear <- qlogis(predicted)
  pearson <- residuals(model, type = "pearson")
  leverage <- model$hat
  cooks <- pearson^2 * leverage / (sum(model$edf) * (1 - leverage)^2)
  data.frame(
    outcome = outcome,
    converged = model$converged & model$outer.info$conv == "full convergence",
    center_sd = 1 / sqrt(model$sp[["s(center)"]]),
    observed_risk = mean(observed),
    predicted_risk = mean(predicted),
    brier_score = mean((observed - predicted)^2),
    c_statistic = c_statistic(observed, predicted),
    calibration_intercept = unname(coef(glm(observed ~ offset(linear), family = binomial))[1]),
    calibration_slope = unname(coef(glm(observed ~ linear, family = binomial))[2]),
    max_abs_pearson_residual = max(abs(pearson)),
    pearson_residuals_above_3 = sum(abs(pearson) > 3),
    max_leverage = max(leverage),
    max_cooks_distance = max(cooks)
  )
}

percentile_interval <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

add_intervals <- function(point, bootstrap) {
  key_point <- paste(point$analysis, point$outcome, point$population, point$metric)
  key_boot <- paste(bootstrap$analysis, bootstrap$outcome, bootstrap$population, bootstrap$metric)
  intervals <- t(sapply(split(bootstrap$value, key_boot), percentile_interval))
  point$lower <- intervals[key_point, 1]
  point$upper <- intervals[key_point, 2]
  point
}

hte_test <- function(point, bootstrap, analysis, outcome, metric) {
  estimates <- sapply(phenotype_labels, function(p) {
    point$value[point$analysis == analysis & point$outcome == outcome &
                   point$population == p & point$metric == metric]
  })
  rows <- bootstrap[bootstrap$analysis == analysis & bootstrap$outcome == outcome &
                      bootstrap$metric == metric & bootstrap$population %in% phenotype_labels, ]
  draws <- reshape(rows[c("replicate", "population", "value")], idvar = "replicate",
                   timevar = "population", direction = "wide")
  draws <- as.matrix(draws[paste0("value.", phenotype_labels)])
  contrast <- rbind(c(-1, 1, 0), c(-1, 0, 1))
  difference <- contrast %*% estimates
  variance <- contrast %*% cov(draws, use = "complete.obs") %*% t(contrast)
  statistic <- as.numeric(t(difference) %*% solve(variance) %*% difference)
  data.frame(analysis = analysis, outcome = outcome, metric = metric, df = 2,
             statistic = statistic, p_value = pchisq(statistic, 2, lower.tail = FALSE))
}

population_p_values <- function(bootstrap, analysis, outcome, metric) {
  null <- if (metric == "ratio") 1 else 0
  sapply(populations, function(p) {
    draws <- bootstrap$value[bootstrap$analysis == analysis & bootstrap$outcome == outcome &
                               bootstrap$population == p & bootstrap$metric == metric]
    2 * min(mean(draws <= null, na.rm = TRUE), mean(draws >= null, na.rm = TRUE))
  })
}

e_value <- function(ratio) {
  ratio <- ifelse(ratio < 1, 1 / ratio, ratio)
  ratio + sqrt(ratio * (ratio - 1))
}

e_value_limit <- function(lower, upper) {
  ifelse(lower <= 1 & upper >= 1, 1, e_value(ifelse(lower > 1, lower, upper)))
}
