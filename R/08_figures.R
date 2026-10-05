save_figure <- function(name, draw, width, height, folder = file.path("results", "figures")) {
  for (type in c("pdf", "svg", "png")) {
    path <- file.path(folder, paste0(name, ".", type))
    if (type == "pdf") cairo_pdf(path, width = width, height = height, family = "Arial")
    if (type == "svg") svg(path, width = width, height = height, family = "Arial")
    if (type == "png") png(path, width = width, height = height, units = "in",
                           res = settings$figure_dpi, type = "cairo", family = "Arial")
    par(las = 1, bty = "n", col.axis = "#2B2B2B", col.lab = "#2B2B2B")
    draw()
    dev.off()
  }
}

format_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

strategy_colors <- unname(palette[c("strategy_1", "strategy_2")])

plot_participant_flow <- function(flow) {
  n <- setNames(flow$n, flow$stage)
  box <- function(x, y, w, h, label, fill, cex = 0.72) {
    rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = "#6B7C85", lwd = 1.2)
    text(x, y, label, cex = cex, col = "#1F2933")
  }
  arrow <- function(x0, y0, x1, y1) arrows(x0, y0, x1, y1, length = 0.08, lwd = 1.3, col = "#596870")
  plot.new()
  plot.window(xlim = c(0, 16), ylim = c(0, 16))
  box(5.2, 14.8, 5.8, 1.2, paste0("Records assessed for phenotyping\nn = ", n[["Records assessed for phenotyping"]]), "#EEF3F6")
  box(12.6, 13.4, 5.2, 1.1, paste0("Age below 18 years\nn = ", n[["Excluded: age below 18 years"]]), "#F4F5F6", 0.66)
  box(5.2, 12.0, 5.8, 1.2, paste0("Adults eligible for phenotyping\nn = ", n[["Adults eligible for phenotyping"]]), "#E6EFF4")
  box(12.6, 10.6, 5.2, 1.3, paste0("Missing data for at least one\nclustering variable\nn = ",
                                   n[["Excluded: missing clustering variables"]]), "#F4F5F6", 0.64)
  box(5.2, 9.2, 5.8, 1.2, paste0("Phenotype population\nn = ", n[["Phenotype population"]]), "#DCEAF3")
  box(12.6, 7.8, 5.2, 1.3, paste0("IABP initiated before V-A ECMO\nn = ",
                                  n[["Excluded: IABP initiated before V-A ECMO"]]), "#F4F5F6", 0.64)
  box(5.2, 6.4, 5.8, 1.2, paste0("Target trial emulation population\nn = ", n[["TTE population"]]), "#DCEAF3")
  box(2.7, 2.6, 4.6, 1.4, paste0("mAFP-first\nn = ", n[["mAFP-first"]]), "#D7ECF7")
  box(8.6, 2.6, 4.8, 1.4, paste0("V-A ECMO with LV unloading\nn = ", n[["V-A ECMO with LV unloading"]]), "#F8E1D8")
  arrow(5.2, 14.2, 5.2, 12.6)
  arrow(8.1, 14.5, 10.0, 13.95)
  arrow(5.2, 11.4, 5.2, 9.8)
  arrow(8.1, 11.7, 10.0, 11.25)
  arrow(5.2, 8.6, 5.2, 7.0)
  arrow(8.1, 8.9, 10.0, 8.45)
  arrow(4.6, 5.8, 3.1, 3.3)
  arrow(5.8, 5.8, 8.1, 3.3)
}

plot_propensity_overlap <- function(propensity) {
  data <- propensity$data
  score <- propensity$propensity
  par(mfrow = c(2, 2), mar = c(4.6, 4.4, 2.4, 0.8), oma = c(0, 0, 2.2, 0))
  for (population in populations) {
    keep <- in_population(data, population)
    a <- density(score[keep & data$treatment == 0], from = 0, to = 1)
    b <- density(score[keep & data$treatment == 1], from = 0, to = 1)
    plot(NA, xlim = c(0, 1), ylim = c(0, max(a$y, b$y) * 1.08), axes = FALSE, main = population,
         xlab = "Propensity for V-A ECMO with LV unloading", ylab = "Density")
    axis(1, at = seq(0, 1, 0.2))
    axis(2)
    polygon(c(a$x, rev(a$x)), c(a$y, rep(0, length(a$y))), col = adjustcolor(strategy_colors[1], 0.25), border = NA)
    polygon(c(b$x, rev(b$x)), c(b$y, rep(0, length(b$y))), col = adjustcolor(strategy_colors[2], 0.25), border = NA)
    lines(a$x, a$y, col = strategy_colors[1], lwd = 2.2)
    lines(b$x, b$y, col = strategy_colors[2], lwd = 2.2)
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("top", legend = strategy_labels, col = strategy_colors, lwd = 2.5, horiz = TRUE, bty = "n", cex = 0.8)
}

draw_favors <- function(reference, higher_is_better) {
  labels <- c("Favors ECMO-first", "Favors mAFP-first")
  if (higher_is_better) labels <- rev(labels)
  usr <- par("usr")
  y <- usr[3] + 0.05 * (usr[4] - usr[3])
  gap <- 0.015 * (usr[2] - usr[1])
  width <- max(strwidth(labels, cex = 0.72))
  arrows(reference - gap, y, reference - gap - width, y, length = 0.05, angle = 25, lwd = 0.9, col = "#6B7C85")
  arrows(reference + gap, y, reference + gap + width, y, length = 0.05, angle = 25, lwd = 0.9, col = "#6B7C85")
  text(reference - gap, y, labels[1], adj = c(1, -0.6), cex = 0.72, col = "#6B7C85")
  text(reference + gap, y, labels[2], adj = c(0, -0.6), cex = 0.72, col = "#6B7C85")
}

plot_effect_forest <- function(table, metric, xlab, header, p_value, p_values, scale = 100, digits = 1,
                               higher_is_better = FALSE) {
  table <- table[match(populations, table$population), ]
  y <- rev(seq_len(nrow(table)))
  estimate <- scale * table[[metric]]
  lower <- scale * table[[paste0(metric, "_lower")]]
  upper <- scale * table[[paste0(metric, "_upper")]]
  reference <- if (grepl("ratio", metric)) 1 else 0
  span <- diff(range(c(lower, upper, reference)))
  ticks <- pretty(c(lower, upper, reference - 0.3 * span, reference + 0.3 * span), n = 5)
  colors <- c(palette[["neutral"]], phenotype_colors)
  layout(matrix(1:2, nrow = 1), widths = c(3.4, 1.8))
  par(mar = c(5.5, 8.6, 2.8, 0.5))
  plot(NA, xlim = range(ticks), ylim = c(0.1, nrow(table) + 0.75), axes = FALSE, xlab = xlab, ylab = "")
  axis(1, at = ticks)
  axis(2, at = y, labels = table$population, tick = FALSE)
  abline(v = reference, col = "#9A9A9A", lty = 3)
  draw_favors(reference, higher_is_better)
  segments(lower, y, upper, y, lwd = 2, col = colors)
  points(estimate, y, pch = 16, cex = 1.15, col = colors)
  mtext(paste("Global phenotype interaction", format_p(p_value)), side = 3, line = 0.45, adj = 0, cex = 0.82)
  par(mar = c(5.5, 0.5, 2.8, 0.5))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0.1, nrow(table) + 0.75))
  labels <- sprintf(paste0("%.", digits, "f (%.", digits, "f to %.", digits, "f)"), estimate, lower, upper)
  text(0, y, labels, adj = 0, cex = 0.8)
  text(0, nrow(table) + 0.55, header, adj = 0, font = 2, cex = 0.8)
  p_column <- max(strwidth(labels, cex = 0.8), strwidth(header, cex = 0.8, font = 2)) + 0.08
  p_values <- p_values[table$population]
  text(p_column, y, ifelse(p_values < 0.001, "< 0.001", sprintf("%.3f", p_values)), adj = 0, cex = 0.8)
  text(p_column, nrow(table) + 0.55, "p-value", adj = 0, font = 2, cex = 0.8)
  layout(1)
}

plot_standardized_cif <- function(curves, table, p_value, p_values, analysis = "Primary") {
  par(mfrow = c(2, 2), mar = c(4.4, 4.4, 2.0, 0.8), oma = c(0.2, 0.2, 2.0, 0))
  for (population in populations) {
    plot(NA, xlim = c(0, settings$horizon), ylim = c(0, 1), axes = FALSE, main = population,
         xlab = "Days from T0", ylab = "Cumulative incidence of in-hospital death")
    axis(1, at = seq(0, 350, 50))
    axis(2, at = seq(0, 1, 0.2))
    for (treatment in 0:1) {
      line <- curves[curves$analysis == analysis & curves$population == population & curves$treatment == treatment, ]
      lines(line$day, line$death, col = strategy_colors[treatment + 1], lwd = 2.3)
    }
    if (population == "Overall") {
      legend("topleft", legend = strategy_labels, col = strategy_colors, lwd = 2.3, bty = "n", cex = 0.72)
    }
    row <- table[table$population == population, ]
    text(settings$horizon, 0.03,
         sprintf("RMTL difference, V-A ECMO minus mAFP-first\n%.1f days (95%% CI %.1f to %.1f), %s",
                 row$rmtl_difference, row$rmtl_difference_lower, row$rmtl_difference_upper,
                 format_p(p_values[[population]])),
         adj = c(1, 0), cex = 0.6, col = palette[["neutral"]])
  }
  mtext(paste("Global phenotype interaction for RMTL:", format_p(p_value)), side = 3, outer = TRUE,
        line = 0.25, cex = 0.82, font = 2)
}

plot_state_occupation <- function(states) {
  colors <- c("#457B9D", "#B08968", "#6A4C93", "#9E2A2B")
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.62))
  for (strategy in strategy_labels) {
    par(mar = c(4, 4, 2.4, 0.8))
    x <- states[states$strategy == strategy, ]
    days <- unique(x$day)
    plot(NA, xlim = c(0, settings$horizon), ylim = c(0, 1), axes = FALSE, main = strategy,
         xlab = "Days from T0", ylab = "State occupation probability")
    axis(1, at = seq(0, 350, 50))
    axis(2, at = seq(0, 1, 0.2))
    lower <- rep(0, length(days))
    for (j in seq_along(state_names)) {
      upper <- lower + x$probability[x$state == state_names[j]]
      polygon(c(days, rev(days)), c(lower, rev(upper)), col = colors[j], border = NA)
      lower <- upper
    }
  }
  par(mar = c(4, 0.2, 2.4, 0.2))
  plot.new()
  legend("left", legend = state_names, fill = colors, border = NA, bty = "n", cex = 0.82, xpd = NA)
  layout(1)
}

plot_scai_heatmap <- function(outcomes) {
  stages <- c("C", "D", "E")
  panels <- list(
    list(variable = "in_hospital_mortality", title = "In-hospital mortality by day 365",
         colors = colorRampPalette(c("#F7FBFF", "#2166AC"))(101)),
    list(variable = "recovery", title = "Recovery at one year",
         colors = colorRampPalette(c("#F7FCF5", "#238B45"))(101))
  )
  par(mfrow = c(1, 2), mar = c(4.2, 7.4, 2.6, 1.0))
  for (panel in panels) {
    plot.new()
    plot.window(xlim = c(0.5, 3.5), ylim = c(0.5, 3.5))
    for (i in seq_along(phenotype_labels)) {
      for (j in seq_along(stages)) {
        row <- outcomes[outcomes$phenotype == phenotype_labels[i] & outcomes$scai == stages[j], ]
        if (row$n == 0) {
          rect(j - 0.47, 3.53 - i, j + 0.47, 4.47 - i, col = "#EFEFEF", border = "white", lwd = 2)
          text(j, 4 - i, "n = 0", cex = 0.84, col = "#1F2933")
          next
        }
        value <- row[[panel$variable]]
        rect(j - 0.47, 3.53 - i, j + 0.47, 4.47 - i, col = panel$colors[round(100 * value) + 1],
             border = "white", lwd = 2)
        text(j, 4 - i, sprintf("%.1f%%\nn = %d", 100 * value, row$n), cex = 0.84,
             col = if (value >= 0.58) "white" else "#1F2933")
      }
    }
    axis(1, at = 1:3, labels = stages, tick = FALSE)
    axis(2, at = 3:1, labels = phenotype_labels, tick = FALSE, cex.axis = 0.82)
    mtext("SCAI stage", side = 1, line = 2.3, cex = 0.9)
    title(main = panel$title, cex.main = 1)
  }
}

plot_sensitivity_forest <- function(table, outcome, xlab) {
  labels <- c("AMICS", "2017 or later", "V-A ECMO with mAFP\nversus mAFP", "V-A ECMO with IABP\nversus mAFP",
              "Excluding delayed\nLV unloading")
  table <- table[table$outcome == outcome & table$population %in% phenotype_labels, ]
  rows <- rev(seq_along(sensitivity_analyses))
  offsets <- setNames(c(0.22, 0, -0.22), phenotype_labels)
  span <- 100 * diff(range(c(table$difference_lower, table$difference_upper, 0)))
  ticks <- pretty(c(100 * table$difference_lower, 100 * table$difference_upper, -0.3 * span, 0.3 * span), n = 6)
  layout(matrix(1:2, nrow = 1), widths = c(4.3, 1.25))
  par(mar = c(6.5, 11.5, 3.7, 0.5))
  plot(NA, xlim = range(ticks), ylim = c(0.2, length(rows) + 1.05), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = ticks)
  axis(2, at = rows, labels = labels, tick = FALSE, cex.axis = 0.72)
  mtext(xlab, side = 1, line = 3.0, cex = 0.9)
  mtext("V-A ECMO with LV unloading minus mAFP-first", side = 1, line = 4.4, cex = 0.78)
  abline(v = 0, col = "#9A9A9A", lty = 3)
  draw_favors(0, outcome == "Recovery")
  for (i in seq_along(sensitivity_analyses)) {
    for (p in phenotype_labels) {
      row <- table[table$analysis == sensitivity_analyses[i] & table$population == p, ]
      y <- rows[i] + offsets[[p]]
      segments(100 * row$difference_lower, y, 100 * row$difference_upper, y, lwd = 1.6, col = phenotype_colors[[p]])
      points(100 * row$difference, y, pch = 16, cex = 0.82, col = phenotype_colors[[p]])
    }
  }
  legend("top", legend = phenotype_labels, col = phenotype_colors, pch = 16, horiz = TRUE, bty = "n", cex = 0.76)
  par(mar = c(6.5, 0.5, 3.7, 0.5))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0.2, length(rows) + 1.05))
  text(0, length(rows) + 0.68, "Global HTE test", adj = 0, font = 2, cex = 0.74)
  for (i in seq_along(sensitivity_analyses)) {
    p <- table$hte_p_value_difference[table$analysis == sensitivity_analyses[i]][1]
    text(0, rows[i], format_p(p), adj = 0, cex = 0.68)
  }
  layout(1)
}

plot_calibration <- function(calibration, label) {
  observed <- calibration$observed
  predicted <- calibration$predicted
  par(mfrow = c(1, 2), mar = c(4.4, 4.4, 2.6, 0.8), oma = c(0, 0, 2.0, 0))
  decile <- cut(predicted, quantile(predicted, seq(0, 1, 0.1)), include.lowest = TRUE)
  plot(tapply(predicted, decile, mean), tapply(observed, decile, mean), xlim = c(0, 1), ylim = c(0, 1),
       pch = 16, cex = 1.1, col = strategy_colors[2], axes = FALSE,
       xlab = "Mean predicted risk", ylab = "Observed risk", main = "Calibration by decile")
  axis(1, at = seq(0, 1, 0.2))
  axis(2, at = seq(0, 1, 0.2))
  abline(0, 1, col = "#9A9A9A", lty = 3)
  thresholds <- seq(1, 0, length.out = 101)
  sensitivity <- sapply(thresholds, function(t) mean(predicted[observed == 1] >= t))
  specificity <- sapply(thresholds, function(t) mean(predicted[observed == 0] < t))
  plot(1 - specificity, sensitivity, type = "l", lwd = 2.2, col = strategy_colors[2], xlim = c(0, 1),
       ylim = c(0, 1), axes = FALSE, xlab = "1 - specificity", ylab = "Sensitivity", main = "ROC curve")
  axis(1, at = seq(0, 1, 0.2))
  axis(2, at = seq(0, 1, 0.2))
  abline(0, 1, col = "#9A9A9A", lty = 3)
  mtext(sprintf("C-statistic = %.3f", c_statistic(observed, predicted)), side = 3, line = 0.3, cex = 0.82)
  mtext(label, side = 3, outer = TRUE, line = 0.25, cex = 0.92, font = 2)
}

plot_phenotype_radar <- function(data) {
  variables <- c("egfr", "lactate", "alt", "platelets")
  standardized <- scale(data[variables])
  means <- t(sapply(phenotype_labels, function(p) colMeans(standardized[data$phenotype == p, ])))
  angles <- seq(0, 2 * pi, length.out = 5)[-5]
  extent <- ceiling(max(abs(means)) * 10) / 10
  radius <- means + extent
  circle <- seq(0, 2 * pi, length.out = 100)
  par(mar = c(2, 2, 2, 2))
  plot.new()
  plot.window(xlim = c(-2, 2) * extent * 1.3, ylim = c(-2, 2) * extent * 1.3, asp = 1)
  for (ring in pretty(c(0, 2 * extent), n = 4)) lines(ring * cos(circle), ring * sin(circle), col = "#D9E2E8")
  lines(extent * cos(circle), extent * sin(circle), col = "#9A9A9A", lty = 2)
  for (i in 1:4) {
    segments(0, 0, 2 * extent * cos(angles[i]), 2 * extent * sin(angles[i]), col = "#D9E2E8")
    text(2.24 * extent * cos(angles[i]), 2.24 * extent * sin(angles[i]), c("eGFR", "Lactate", "ALT", "Platelets")[i], cex = 0.82)
  }
  for (i in 1:3) {
    x <- radius[i, ] * cos(angles)
    y <- radius[i, ] * sin(angles)
    polygon(c(x, x[1]), c(y, y[1]), border = phenotype_colors[i], lwd = 2.2)
    points(x, y, col = phenotype_colors[i], pch = 16, cex = 0.9)
  }
  legend("bottomright", legend = phenotype_labels, col = phenotype_colors, lwd = 2.2, pch = 16, bty = "n", cex = 0.78)
  title("Phenotype metabolic profile, standardized means")
  mtext("Dashed ring marks the population average (z = 0)", side = 1, line = -1, cex = 0.68, col = "#6B7C85")
}

plot_phenotype_mortality <- function(outcomes, cumulative_incidence) {
  par(mfrow = c(1, 2), mar = c(5.4, 4.4, 2.6, 0.8))
  rows <- outcomes[outcomes$scai == "All" & outcomes$phenotype %in% phenotype_labels, ]
  heights <- 100 * rows$in_hospital_mortality[match(phenotype_labels, rows$phenotype)]
  bars <- barplot(heights, col = phenotype_colors, ylim = c(0, max(heights) * 1.25), names.arg = phenotype_labels,
                  cex.names = 0.72, ylab = "In-hospital mortality by day 365 (%)")
  text(bars, heights, sprintf("%.1f%%", heights), pos = 3, cex = 0.82)
  title("In-hospital mortality by phenotype", cex.main = 0.92)
  curves <- cumulative_incidence$curves
  plot(NA, xlim = c(0, settings$horizon), ylim = c(0, 1), axes = FALSE, xlab = "Days from T0",
       ylab = "Cumulative incidence of in-hospital death")
  axis(1, at = seq(0, 350, 50))
  axis(2, at = seq(0, 1, 0.2))
  for (p in phenotype_labels) {
    line <- curves[curves$series == paste(p, 1), ]
    lines(line$day, line$cumulative_incidence, type = "s", col = phenotype_colors[[p]], lwd = 2.3)
  }
  legend("topleft", legend = phenotype_labels, col = phenotype_colors, lwd = 2.3, bty = "n", cex = 0.7)
  text(settings$horizon, 1, paste("Gray test", format_p(cumulative_incidence$gray_test$median_p_value)),
       adj = c(1, 1), cex = 0.78)
  title("Aalen-Johansen estimate by phenotype", cex.main = 0.92)
}

plot_cluster_quality <- function(quality) {
  metrics <- c(within_ss = "Within-cluster sum of squares", silhouette = "Mean silhouette width",
               calinski_harabasz = "Calinski-Harabasz index", davies_bouldin = "Davies-Bouldin index",
               ambiguous_proportion = "Proportion of ambiguous clustering")
  par(mfrow = c(2, 3), mar = c(4.2, 4.6, 2.5, 1))
  for (m in names(metrics)) {
    plot(quality$k, quality[[m]], type = "b", pch = 16, col = palette[["neutral"]], axes = FALSE,
         xlab = "Number of clusters", ylab = "", main = metrics[[m]], cex.main = 0.9)
    axis(1, at = 2:7)
    axis(2)
    abline(v = 3, col = "#9A9A9A", lty = 3)
  }
}

plot_scai_violin <- function(data) {
  stage <- as.numeric(data$scai)
  densities <- lapply(phenotype_labels, function(p) {
    x <- stage[data$phenotype == p]
    density(x, bw = 0.15, from = min(x) - 0.45, to = max(x) + 0.45)
  })
  top <- max(sapply(densities, function(d) max(d$y)))
  par(mar = c(3.6, 4.6, 1.5, 0.8))
  plot(NA, xlim = c(0.5, 3.5), ylim = c(0.55, 3.45), axes = FALSE, xlab = "", ylab = "SCAI stage")
  axis(1, at = 1:3, labels = phenotype_labels, tick = FALSE)
  axis(2, at = 1:3, labels = c("C", "D", "E"))
  for (i in 1:3) {
    d <- densities[[i]]
    width <- 0.42 * d$y / top
    polygon(c(i - width, rev(i + width)), c(d$x, rev(d$x)),
            col = adjustcolor(phenotype_colors[i], 0.45), border = phenotype_colors[i], lwd = 1.2)
  }
}

plot_unadjusted_cif <- function(aalen_johansen) {
  par(mfrow = c(2, 2), mar = c(4.4, 4.4, 2.0, 0.8))
  for (population in populations) {
    plot(NA, xlim = c(0, settings$horizon), ylim = c(0, 1), axes = FALSE, main = population,
         xlab = "Days from T0", ylab = "Cumulative incidence of in-hospital death")
    axis(1, at = seq(0, 350, 50))
    axis(2, at = seq(0, 1, 0.2))
    for (treatment in 1:2) {
      series <- if (population == "Overall") paste(strategy_labels[treatment], 1) else
        paste(strategy_labels[treatment], "|", population, 1)
      line <- aalen_johansen[aalen_johansen$series == series, ]
      lines(line$day, line$cumulative_incidence, type = "s", col = strategy_colors[treatment], lwd = 2.3)
    }
    if (population == "Overall") {
      legend("topleft", legend = strategy_labels, col = strategy_colors, lwd = 2.3, bty = "n", cex = 0.72)
    }
  }
}

plot_rescue_cif <- function(rescue) {
  par(mar = c(4.6, 4.6, 2.6, 0.8))
  plot(NA, xlim = c(0, 60), ylim = c(0, 0.4), axes = FALSE, xlab = "Days from T0",
       ylab = "Cumulative incidence of rescue V-A ECMO", main = "Rescue V-A ECMO in the mAFP-first strategy",
       cex.main = 0.92)
  axis(1, at = seq(0, 60, 10))
  axis(2, at = seq(0, 0.4, 0.1))
  colors <- c(Overall = palette[["neutral"]], phenotype_colors)
  for (series in names(colors)) {
    line <- rescue[rescue$series == paste(series, 1), ]
    lines(line$day, line$cumulative_incidence, type = "s", col = colors[[series]], lwd = 2.3)
  }
  legend("topleft", legend = names(colors), col = colors, lwd = 2.3, bty = "n", cex = 0.75)
}

plot_tsne <- function(tsne) {
  par(mar = c(4.4, 4.4, 2.6, 0.8))
  plot(tsne$dimension_1, tsne$dimension_2, col = phenotype_colors[tsne$phenotype], pch = 16, cex = 0.6,
       xlab = "t-SNE 1", ylab = "t-SNE 2", main = "t-SNE of the clustering space")
  legend("topright", legend = phenotype_labels, col = phenotype_colors, pch = 16, bty = "n", cex = 0.78)
}

save_all_figures <- function(output) {
  with(output, {
    mortality <- effect_table(results, hte, "Primary", "In-hospital death", c("difference", "ratio"))
    recovery <- effect_table(results, hte, "Primary", "Recovery", c("difference", "ratio"))
    competing <- effect_table(results, hte, "Primary", "In-hospital death, competing risks", c("difference", "rmtl_difference"))
    sensitivity <- sensitivity_table(results, hte)
    mortality_p <- population_p_values(bootstrap, "Primary", "In-hospital death", "difference")
    recovery_p <- population_p_values(bootstrap, "Primary", "Recovery", "difference")
    rmtl_p <- population_p_values(bootstrap, "Primary", "In-hospital death, competing risks", "rmtl_difference")

    save_figure("figure_1_participant_flow", function() plot_participant_flow(flow), 9.0, 7.2)
    save_figure("figure_2_propensity_overlap", function() plot_propensity_overlap(propensity), 8.4, 7.2)
    save_figure("figure_3_in_hospital_mortality_risk_difference", function() {
      plot_effect_forest(mortality, "difference",
                         "In-hospital mortality risk difference, percentage points\nV-A ECMO with LV unloading minus mAFP-first",
                         "Percentage points (95% CI)", mortality$hte_p_value_difference[1], mortality_p)
    }, 9.2, 5.2)
    save_figure("figure_4_standardized_cumulative_incidence", function() {
      plot_standardized_cif(point$curves, competing, competing$hte_p_value_rmtl_difference[1], rmtl_p)
    }, 8.0, 7.2)
    save_figure("figure_5_recovery_risk_difference", function() {
      plot_effect_forest(recovery, "difference",
                         "Recovery risk difference, percentage points\nV-A ECMO with LV unloading minus mAFP-first",
                         "Percentage points (95% CI)", recovery$hte_p_value_difference[1], recovery_p,
                         higher_is_better = TRUE)
    }, 9.2, 5.2)
    save_figure("figure_6_state_occupation", function() plot_state_occupation(states), 11.0, 5.2)
    save_figure("figure_7_scai_phenotype_outcomes", function() plot_scai_heatmap(phenotype_results$outcomes), 9.4, 4.8)
    save_figure("figure_s1_mortality_sensitivity", function() {
      plot_sensitivity_forest(sensitivity, "In-hospital death", "In-hospital mortality risk difference, percentage points")
    }, 9.5, 7.6)
    save_figure("figure_s2_recovery_sensitivity", function() {
      plot_sensitivity_forest(sensitivity, "Recovery", "Recovery risk difference, percentage points")
    }, 9.5, 7.6)
    save_figure("figure_s3_mortality_calibration_discrimination", function() {
      plot_calibration(calibration$in_hospital_death, "In-hospital mortality by day 365")
    }, 9.0, 4.6)
    save_figure("figure_s4_recovery_calibration_discrimination", function() {
      plot_calibration(calibration$recovery_365, "Recovery at one year")
    }, 9.0, 4.6)
    save_figure("figure_s5_phenotype_metabolic_profile", function() plot_phenotype_radar(phenotype_data), 6.5, 6.0)
    save_figure("figure_s6_phenotype_in_hospital_mortality", function() {
      plot_phenotype_mortality(phenotype_results$outcomes, phenotype_results$cumulative_incidence)
    }, 9.0, 4.8)
    save_figure("figure_s7_cluster_quality", function() plot_cluster_quality(phenotypes$quality), 9.0, 6.0)
    save_figure("figure_s8_scai_phenotype_violin", function() plot_scai_violin(phenotype_data), 7.0, 5.4)
    save_figure("figure_s9_unadjusted_cumulative_incidence", function() plot_unadjusted_cif(point$aalen_johansen), 8.0, 7.2)
    save_figure("figure_s10_rescue_ecmo_cumulative_incidence", function() plot_rescue_cif(point$rescue), 6.5, 5.2)
    save_figure("phenotype_tsne", function() plot_tsne(phenotypes$tsne), 7.0, 6.0, folder = "private")
  })
}
