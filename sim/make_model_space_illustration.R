source(file.path("sim", "design.R"))
source(file.path("sim", "exact_bma.R"))
source(file.path("sim", "repfam_bma.R"))

output_dir <- file.path("sim", "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(920007L)
data <- simulate_block_regression(
  n_train = 45,
  n_test = 200,
  K = 5,
  m = 3,
  rho = 0.98,
  active_groups = c(1, 3),
  active_reps = c(1, 2),
  beta_values = c(1.0, 0.8),
  sigma = 1.5
)

exact <- fit_exact_bma(
  X = data$X_train,
  y = data$y_train,
  theta = 0.10,
  tau2 = 4,
  a0 = 1,
  b0 = 1
)
repfam <- fit_repfam_bma(
  exact_fit = exact,
  group_id = data$group_id,
  lambda = 1.4
)

support_summary <- posterior_summary(exact$posterior)
family_summary <- posterior_summary(repfam$family_posterior)
oracle_index <- family_index(repfam$families, data$active_groups)
oracle_retained <- repfam$alpha[oracle_index]

top_support <- order(exact$posterior, decreasing = TRUE)[seq_len(45)]
top_family <- order(repfam$family_posterior, decreasing = TRUE)[seq_len(12)]

family_label <- function(row) {
  active <- which(row == 1L)
  if (length(active) == 0L) {
    return("{}")
  }
  paste0("{", paste(active, collapse = ","), "}")
}

draw_cluster_panel <- function() {
  set.seed(42L)
  centers <- matrix(
    c(0.25, 0.67,
      0.47, 0.42,
      0.68, 0.66,
      0.77, 0.34),
    ncol = 2,
    byrow = TRUE
  )
  family_mass <- c(0.34, 0.25, 0.18, 0.08)
  family_mass <- family_mass / sum(family_mass)
  cols <- c("#2f6f9f", "#ca6f1e", "#5b8e3e", "#7b5aa6")
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "Model space coarse graining"
  )
  box(col = "gray75")
  for (g in seq_len(nrow(centers))) {
    theta <- seq(0, 2 * pi, length.out = 120)
    polygon(
      centers[g, 1] + 0.17 * cos(theta),
      centers[g, 2] + 0.13 * sin(theta),
      border = adjustcolor(cols[g], 0.45),
      col = adjustcolor(cols[g], 0.08),
      lwd = 1.2
    )
    n_points <- c(28, 24, 20, 14)[g]
    xy <- cbind(
      rnorm(n_points, centers[g, 1], 0.055),
      rnorm(n_points, centers[g, 2], 0.045)
    )
    point_size <- rexp(n_points, rate = 5)
    point_size <- 0.45 + 1.6 * point_size / max(point_size)
    points(
      xy[, 1],
      xy[, 2],
      pch = 19,
      cex = point_size,
      col = adjustcolor(cols[g], 0.78)
    )
    text(
      centers[g, 1],
      centers[g, 2] - 0.18,
      labels = paste0("region ", g),
      cex = 0.75,
      col = cols[g]
    )
  }
  text(
    0.5,
    0.05,
    "many support states become a few region states",
    cex = 0.78,
    col = "gray30"
  )
}

draw_support_panel <- function() {
  weights <- exact$posterior[top_support]
  barplot(
    weights,
    col = "#7f8fa6",
    border = NA,
    main = "Unrestricted support posterior",
    xlab = "support rank",
    ylab = "posterior probability",
    names.arg = rep("", length(weights))
  )
  mtext(
    sprintf("95%% count = %d", support_summary$n_mass),
    side = 3,
    line = -1.1,
    adj = 0.98,
    cex = 0.82,
    col = "gray25"
  )
}

draw_family_panel <- function() {
  weights <- repfam$family_posterior[top_family]
  labels <- vapply(
    top_family,
    function(i) family_label(repfam$families[i, ]),
    character(1)
  )
  barplot(
    weights,
    col = "#2f6f9f",
    border = NA,
    main = "Support region summary",
    xlab = "region state",
    ylab = "posterior probability",
    names.arg = labels,
    las = 2,
    cex.names = 0.65
  )
  mtext(
    sprintf("95%% count = %d, best retained = %.2f",
            family_summary$n_mass,
            oracle_retained),
    side = 3,
    line = -1.1,
    adj = 0.98,
    cex = 0.82,
    col = "gray25"
  )
}

draw_figure <- function(path, device) {
  device(path, width = 10.5, height = 3.7)
  op <- par(mfrow = c(1, 3), mar = c(4.3, 4.1, 2.5, 0.8), oma = c(0, 0, 0, 0))
  on.exit({
    par(op)
    dev.off()
  })
  draw_cluster_panel()
  draw_support_panel()
  draw_family_panel()
}

draw_figure(file.path(output_dir, "model_space_compression_illustration.svg"), svg)
draw_figure(file.path(output_dir, "model_space_compression_illustration.pdf"), pdf)

writeLines(
  c(
    sprintf("Support 95%% count: %d", support_summary$n_mass),
    sprintf("Region 95%% count: %d", family_summary$n_mass),
    sprintf("Oracle retained mass: %.4f", oracle_retained),
    sprintf("Top support probability: %.4f", max(exact$posterior)),
    sprintf("Top family probability: %.4f", max(repfam$family_posterior))
  ),
  file.path(output_dir, "model_space_compression_illustration_report.txt")
)
