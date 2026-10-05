clustering_variables <- c("egfr", "lactate", "alt", "platelets")

sample_skewness <- function(x) mean((x - mean(x))^3) / sd(x)^3

rescale_01 <- function(x) (x - min(x)) / (max(x) - min(x))

calinski_harabasz <- function(x, cluster) {
  overall <- colMeans(x)
  between <- 0
  within <- 0
  for (g in unique(cluster)) {
    members <- x[cluster == g, , drop = FALSE]
    center <- colMeans(members)
    between <- between + nrow(members) * sum((center - overall)^2)
    within <- within + sum(sweep(members, 2, center)^2)
  }
  k <- length(unique(cluster))
  (between / (k - 1)) / (within / (nrow(x) - k))
}

davies_bouldin <- function(x, cluster) {
  groups <- sort(unique(cluster))
  centers <- t(sapply(groups, function(g) colMeans(x[cluster == g, , drop = FALSE])))
  spread <- sapply(seq_along(groups), function(i) {
    members <- x[cluster == groups[i], , drop = FALSE]
    sqrt(mean(rowSums(sweep(members, 2, centers[i, ])^2)))
  })
  distance <- as.matrix(dist(centers))
  ratios <- outer(spread, spread, "+") / distance
  diag(ratios) <- NA
  mean(apply(ratios, 1, max, na.rm = TRUE))
}

ambiguous_proportion <- function(consensus_matrix) {
  values <- consensus_matrix[upper.tri(consensus_matrix)]
  mean(values > 0.1 & values < 0.9)
}

clustering_missingness <- function(cohort) {
  adults <- cohort[cohort$adult, clustering_variables]
  by_variable <- data.frame(
    variable = clustering_variables,
    n_missing = colSums(is.na(adults)),
    percent_missing = 100 * colMeans(is.na(adults))
  )
  patterns <- apply(is.na(adults), 1, function(row) paste(clustering_variables[row], collapse = " + "))
  patterns <- patterns[patterns != ""]
  pattern_table <- as.data.frame(table(pattern = patterns))
  names(pattern_table)[2] <- "n_patients"
  list(by_variable = by_variable, pattern = pattern_table[order(-pattern_table$n_patients), ])
}

derive_phenotypes <- function(cohort) {
  input <- cohort[cohort$phenotype_population, c("row_id", clustering_variables)]
  cache <- file.path("private", "phenotypes.rds")
  if (file.exists(cache)) {
    previous <- readRDS(cache)
    if (identical(previous$input, input)) return(previous)
  }

  transformed <- input[clustering_variables]
  transformation <- setNames(rep("none", 4), clustering_variables)
  for (v in clustering_variables) {
    if (abs(sample_skewness(transformed[[v]])) > 1) {
      transformed[[v]] <- log1p(transformed[[v]] - min(transformed[[v]]))
      transformation[v] <- "shifted log1p"
    }
  }
  scaled <- as.matrix(as.data.frame(lapply(transformed, rescale_01)))

  consensus <- ConsensusClusterPlus::ConsensusClusterPlus(
    t(scaled), maxK = 7, reps = 1000, pItem = 0.8, pFeature = 1,
    clusterAlg = "km", distance = "euclidean", seed = 42,
    title = file.path("private", "consensus"), plot = "png", verbose = FALSE
  )

  distances <- dist(scaled)
  quality <- do.call(rbind, lapply(2:7, function(k) {
    set.seed(1)
    fit <- kmeans(scaled, centers = k, nstart = 50, iter.max = 100)
    data.frame(
      k = k,
      within_ss = fit$tot.withinss,
      silhouette = mean(cluster::silhouette(fit$cluster, distances)[, 3]),
      calinski_harabasz = calinski_harabasz(scaled, fit$cluster),
      davies_bouldin = davies_bouldin(scaled, fit$cluster),
      ambiguous_proportion = ambiguous_proportion(consensus[[k]]$consensusMatrix)
    )
  }))

  set.seed(1)
  final <- kmeans(scaled, centers = 3, nstart = 50, iter.max = 100)
  profile <- aggregate(as.data.frame(scale(transformed)), list(cluster = final$cluster), mean)
  cardiometabolic <- profile$cluster[which.max(profile$lactate)]
  others <- profile[profile$cluster != cardiometabolic, ]
  cardiorenal <- others$cluster[which.min(others$egfr)]
  phenotype <- ifelse(final$cluster == cardiometabolic, "Cardiometabolic",
                      ifelse(final$cluster == cardiorenal, "Cardiorenal", "Non-congestive"))

  set.seed(1)
  tsne <- Rtsne::Rtsne(scaled, perplexity = 30, pca = TRUE, check_duplicates = FALSE)

  result <- list(
    input = input,
    assignments = data.frame(row_id = input$row_id, phenotype = phenotype),
    quality = quality,
    transformation = data.frame(variable = clustering_variables, transformation = unname(transformation)),
    tsne = data.frame(dimension_1 = tsne$Y[, 1], dimension_2 = tsne$Y[, 2], phenotype = phenotype)
  )
  saveRDS(result, cache)
  result
}
