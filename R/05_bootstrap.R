start_cluster <- function() {
  cluster <- parallel::makeCluster(settings$workers)
  parallel::clusterCall(cluster, setwd, getwd())
  parallel::clusterEvalQ(cluster, {
    library(mgcv)
    library(splines)
    source("R/00_settings.R")
    source("R/03_imputation.R")
    source("R/04_models.R")
    source("R/05_bootstrap.R")
    NULL
  })
  cluster
}

point_estimates <- function(cluster, imputed) {
  runs <- parallel::parLapply(cluster, imputed, estimate_all, details = TRUE)
  average <- function(name, keys) {
    stacked <- do.call(rbind, lapply(runs, `[[`, name))
    values <- setdiff(names(stacked), keys)
    aggregate(stacked[values], stacked[keys], mean)
  }
  estimates <- average("estimates", c("analysis", "outcome", "population", "metric"))
  list(
    estimates = recompute_ratios(estimates),
    curves = average("curves", c("analysis", "population", "treatment", "day")),
    diagnostics = average("diagnostics", "outcome"),
    aalen_johansen = average("aalen_johansen", c("scope", "series", "day")),
    rescue = average("rescue", c("series", "day"))
  )
}

recompute_ratios <- function(estimates) {
  key <- paste(estimates$analysis, estimates$outcome, estimates$population)
  risk_1 <- setNames(estimates$value[estimates$metric == "risk_1"], key[estimates$metric == "risk_1"])
  risk_2 <- setNames(estimates$value[estimates$metric == "risk_2"], key[estimates$metric == "risk_2"])
  ratio <- estimates$metric == "ratio"
  estimates$value[ratio] <- risk_2[key[ratio]] / risk_1[key[ratio]]
  estimates
}

run_replicate <- function(b, tte) {
  file <- file.path(settings$bootstrap_folder, sprintf("replicate_%04d.rds", b))
  if (file.exists(file)) return(file)
  set.seed(settings$seed + b)
  rows <- unlist(lapply(split(seq_len(nrow(tte)), tte$center),
                        function(i) i[sample.int(length(i), replace = TRUE)]))
  imputed <- impute_cohort(tte[rows, ], m = settings$imputations, seed = settings$seed + b)
  stacked <- do.call(rbind, lapply(imputed, function(d) estimate_all(d)$estimates))
  estimates <- recompute_ratios(aggregate(stacked["value"], stacked[c("analysis", "outcome", "population", "metric")], mean))
  estimates$replicate <- b
  saveRDS(estimates, file)
  file
}

run_bootstrap <- function(cluster, tte) {
  dir.create(settings$bootstrap_folder, recursive = TRUE, showWarnings = FALSE)
  files <- parallel::parLapplyLB(cluster, seq_len(settings$bootstrap_replicates),
                                 run_replicate, tte = tte)
  do.call(rbind, lapply(unlist(files), readRDS))
}
