source(file.path("sim", "design.R"))
source(file.path("sim", "exact_bma.R"))
source(file.path("sim", "src", "family_dictionary.R"))
source(file.path("sim", "src", "adaptive_support_kernel_compression.R"))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) sub("^--mode=", "", args[1]) else "smoke"
if (!mode %in% c("smoke", "medium", "full")) stop("mode must be smoke, medium, or full")

out_root <- file.path("sim", "output")
table_dir <- file.path(out_root, "tables")
fig_dir <- file.path(out_root, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260524)

rep_count <- switch(mode, smoke = 1L, medium = 3L, full = 30L)
rho_grid <- switch(mode, smoke = 0.95, medium = c(0.8, 0.95), full = c(0.5, 0.8, 0.95))
scenario_grid <- switch(mode,
  smoke = c("one_representative"),
  medium = c("one_representative", "multi_representative", "weak_signal"),
  full = c("one_representative", "multi_representative", "weak_signal", "noisy_grouping")
)

simulate_custom <- function(scenario, rho, seed, K = if (mode == "full") 5L else 4L, m = 3L) {
  set.seed(seed)
  n_train <- 70L
  n_test <- 220L
  p <- K * m
  group_id <- make_group_id(K, m)
  X_train_raw <- simulate_block_design(n_train, K, m, rho)
  X_test_raw <- simulate_block_design(n_test, K, m, rho)
  z <- standardize_with_training(X_train_raw, X_test_raw)
  beta <- rep(0, p)
  sigma <- 1.35
  idx3 <- if (K >= 5L) c(1, 7, 14) else c(1, 7, 11)
  if (scenario == "one_representative" || scenario == "noisy_grouping") {
    beta[idx3] <- c(1.4, 1.1, 0.9)
  } else if (scenario == "multi_representative") {
    beta[c(1, 2, idx3[2], idx3[3])] <- c(0.9, 0.75, 1.1, 0.9)
  } else if (scenario == "weak_signal") {
    beta[idx3] <- c(0.75, 0.65, 0.55)
    sigma <- 1.8
  }
  if (scenario == "noisy_grouping") {
    group_id <- sample(group_id)
  }
  y_train_raw <- as.numeric(z$X_train %*% beta + rnorm(n_train, sd = sigma))
  y_test_raw <- as.numeric(z$X_test %*% beta + rnorm(n_test, sd = sigma))
  center <- mean(y_train_raw)
  list(
    X_train = z$X_train,
    X_test = z$X_test,
    y_train = y_train_raw - center,
    y_test = y_test_raw - center,
    group_id = group_id,
    beta = beta,
    scenario = scenario,
    rho = rho,
    K = K,
    m = m
  )
}

mixture_weights_from_q <- function(fit, W, alpha, q) {
  h <- mixture_h(q, W, alpha)
  w <- fit$posterior * h
  w / sum(w)
}

predict_from_q <- function(fit, dat, W, alpha, q) {
  predict_mixture(
    fit,
    X_train = dat$X_train,
    y_train = dat$y_train,
    X_test = dat$X_test,
    y_test = dat$y_test,
    weights = mixture_weights_from_q(fit, W, alpha, q)
  )
}

evaluate_kernel_fit <- function(method, fit, dat, W, alpha, q, costs, reference_pred,
                                extra = list()) {
  d <- mixture_distortions(q, W, alpha, fit$posterior)
  pred <- predict_from_q(fit, dat, W, alpha, q)
  entropy <- -sum(ifelse(q > 0, q * log(q), 0))
  out <- data.frame(
    method = method,
    tv = d$tv,
    fkl = d$kl_base_to_compressed,
    rkl = d$kl_compressed_to_base,
    mean_h_error = abs(d$mean_h - 1),
    expected_code = sum(q * costs),
    q0 = q[1],
    q_effective_kernels = exp(entropy),
    q_kernel95 = which(cumsum(sort(q, decreasing = TRUE)) >= 0.95)[1],
    rmse = pred$rmse,
    logscore = pred$mean_log_score,
    rmse_gap = pred$rmse - reference_pred$rmse,
    logscore_gap = pred$mean_log_score - reference_pred$mean_log_score,
    stringsAsFactors = FALSE
  )
  if (length(extra)) out <- cbind(out, as.data.frame(extra, stringsAsFactors = FALSE))
  out
}

fit_fixed_hard <- function(dat, fit, capacity = 1L) {
  active_sets <- candidate_active_sets_from_supports(fit$supports, dat$group_id, fit$posterior, max_sets = 32L)
  dict <- make_representative_dictionary(dat$group_id, active_sets, capacity = capacity, include_safety = TRUE)
  W <- family_membership_matrix(fit$supports, dict)
  alpha <- estimate_family_alpha(W, fit$posterior, alpha_floor = 1e-8)$alpha_truncated
  costs <- dict$families$cost
  qfit <- optimize_family_mixture(W, alpha, costs, fit$posterior, beta = 0.02, tau = 1e-3, distortion = "fkl", q0_min = 1e-3)
  list(W = W, alpha = alpha, costs = costs, q = qfit$q, n_kernels = ncol(W), kkt_residual = qfit$kkt_residual, converged = qfit$converged)
}

fit_topm <- function(fit, top_m = 32L) {
  dict <- make_support_atom_dictionary(fit$supports, fit$posterior, top_m = top_m, include_safety = TRUE)
  W <- family_membership_matrix(fit$supports, dict)
  alpha <- estimate_family_alpha(W, fit$posterior, alpha_floor = 1e-8)$alpha_truncated
  costs <- dict$families$cost
  qfit <- optimize_family_mixture(W, alpha, costs, fit$posterior, beta = 0.02, tau = 1e-3, distortion = "fkl", q0_min = 1e-3)
  list(W = W, alpha = alpha, costs = costs, q = qfit$q, n_kernels = ncol(W), kkt_residual = qfit$kkt_residual, converged = qfit$converged)
}

run_exact_identity <- function() {
  dat <- simulate_custom("one_representative", 0.9, 1001)
  fit <- fit_exact_bma(dat$X_train, dat$y_train, theta = 0.08, tau2 = 4, a0 = 1, b0 = 1)
  kernels <- c(
    list(safety_kernel(cost = 30)),
    generate_active_set_kernel_pool(fit$supports, dat$group_id, fit$posterior, max_sets = 8L, capacities = 1L, rho_grid = c(1), include_hard = TRUE),
    generate_posterior_cluster_kernels(fit$supports, fit$posterior, dat$group_id, top_medoids = 4L, rho_grid = c(1)),
    generate_interval_kernel_pool(ncol(fit$supports), lengths = c(3, 4, 6), rho_grid = c(1), capacities = c(1L, 2L))
  )
  dict <- new_support_kernel_dictionary(dedupe_kernels(kernels))
  W <- support_kernel_weight_matrix(fit$supports, dict, dat$group_id)
  alpha <- estimate_kernel_alpha(W, fit$posterior, alpha_floor = 1e-8)$alpha_truncated
  costs <- kernel_costs(dict)
  qfit <- optimize_family_mixture(W, alpha, costs, fit$posterior, beta = 0.02, tau = 1e-3, distortion = "fkl", q0_min = 1e-3)
  err <- verify_mixture_identities_exact(fit$posterior, W, alpha, qfit$q)
  err$mode <- mode
  err$kernels <- ncol(W)
  err$supports <- nrow(W)
  err$q0 <- qfit$q[1]
  err$tv <- qfit$distortions$tv
  err$fkl <- qfit$distortions$kl_base_to_compressed
  err$rkl <- qfit$distortions$kl_compressed_to_base
  write.csv(err, file.path(table_dir, "table_soft_kernel_identity.csv"), row.names = FALSE)
  err
}

run_topm_separation <- function() {
  rows <- list()
  id <- 0L
  for (m in c(3, 5, 8, 12)) {
    for (s in c(2, 3, 4, 5)) {
      M_star <- m^s
      M_grid <- unique(pmax(1, round(exp(seq(0, log(M_star), length.out = 40)))))
      id <- id + 1L
      rows[[id]] <- topm_atom_tv_uniform(m, s, M_grid)
    }
  }
  out <- do.call(rbind, rows)
  out$mode <- mode
  write.csv(out, file.path(table_dir, "table_topm_atom_separation.csv"), row.names = FALSE)
  pdf(file.path(fig_dir, "fig_topm_atom_separation.pdf"), width = 6.4, height = 4.2)
  keep <- subset(out, active_groups %in% c(2, 4, 5) & group_size %in% c(3, 8))
  plot(
    NA,
    xlim = range(keep$atoms),
    ylim = c(0, 1),
    log = "x",
    xlab = "Stored support atoms M",
    ylab = "TV distortion",
    main = "Top-M atom truncation under exchangeable redundancy"
  )
  cols <- c("3" = "steelblue4", "8" = "firebrick3")
  for (key in unique(interaction(keep$group_size, keep$active_groups))) {
    z <- keep[interaction(keep$group_size, keep$active_groups) == key, ]
    lines(z$atoms, z$topM_tv, col = cols[as.character(z$group_size[1])], lwd = 1.5, lty = z$active_groups[1] - 1)
  }
  abline(h = 0, lty = 2)
  legend("topright", legend = c("m=3", "m=8"), col = cols, lwd = 1.5, bty = "n")
  grid()
  dev.off()
  out
}

run_exact_adaptive <- function() {
  rows <- list()
  details <- list()
  id <- 0L
  did <- 0L
  for (scenario in scenario_grid) {
    for (rho in rho_grid) {
      for (rep_id in seq_len(rep_count)) {
        seed <- 20260524 + rep_id + round(1000 * rho) + match(scenario, scenario_grid) * 10000
        dat <- simulate_custom(scenario, rho, seed)
        fit <- fit_exact_bma(dat$X_train, dat$y_train, theta = 0.08, tau2 = 4, a0 = 1, b0 = 1)
        ref <- predict_mixture(fit, dat$X_train, dat$y_train, dat$X_test, dat$y_test, weights = fit$posterior)
        fixed <- fit_fixed_hard(dat, fit, capacity = if (scenario == "multi_representative") 2L else 1L)
        topm <- fit_topm(fit, top_m = if (mode == "smoke") 12L else 32L)
        pool <- make_default_kernel_pool(
          supports = fit$supports,
          weights = fit$posterior,
          group_id = dat$group_id,
          X = dat$X_train,
          include_intervals = TRUE,
          include_graph = TRUE,
          include_clusters = TRUE,
          include_active_sets = TRUE,
          mode = mode
        )
        adapt <- learn_adaptive_kernel_dictionary(
          supports = fit$supports,
          weights = fit$posterior,
          group_id = dat$group_id,
          X = dat$X_train,
          candidate_kernels = pool,
          beta = 0.02,
          tau = 1e-3,
          max_iter = if (mode == "smoke") 5L else 12L,
          eta = 1e-4,
          mode = mode
        )
        id <- id + 1L
        rows[[id]] <- rbind(
          evaluate_kernel_fit("fixed_hard_dictionary", fit, dat, fixed$W, fixed$alpha, fixed$q, fixed$costs, ref, list(n_kernels = fixed$n_kernels, kkt_residual = fixed$kkt_residual, active_kernels_001 = NA_real_)),
          evaluate_kernel_fit("adaptive_support_kernel", fit, dat, adapt$W, adapt$alpha, adapt$q, kernel_costs(adapt$dictionary), ref, list(n_kernels = adapt$summary$n_kernels, kkt_residual = adapt$summary$kkt_residual, active_kernels_001 = adapt$summary$active_kernels_001)),
          evaluate_kernel_fit("topM_support_atoms", fit, dat, topm$W, topm$alpha, topm$q, topm$costs, ref, list(n_kernels = topm$n_kernels, kkt_residual = topm$kkt_residual, active_kernels_001 = NA_real_))
        )
        rows[[id]]$scenario <- scenario
        rows[[id]]$rho <- rho
        rows[[id]]$replication_id <- rep_id
        rows[[id]]$mode <- mode
        if (nrow(adapt$trace)) {
          did <- did + 1L
          details[[did]] <- cbind(scenario = scenario, rho = rho, replication_id = rep_id, adapt$trace)
        }
      }
    }
  }
  detail <- do.call(rbind, rows)
  write.csv(detail, file.path(table_dir, "table_adaptive_kernel_exact_detail.csv"), row.names = FALSE)
  if (length(details)) {
    write.csv(do.call(rbind, details), file.path(table_dir, "table_adaptive_kernel_column_trace.csv"), row.names = FALSE)
  }
  aggregate_summary(detail, c("scenario", "rho", "method"), file.path(table_dir, "table_adaptive_kernel_exact.csv"))
}

aggregate_summary <- function(detail, group_cols, outfile) {
  metrics <- c("tv", "fkl", "rkl", "q0", "expected_code", "q_effective_kernels", "q_kernel95", "rmse_gap", "logscore_gap", "n_kernels", "kkt_residual", "active_kernels_001")
  metrics <- intersect(metrics, names(detail))
  se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
  parts <- split(detail, interaction(detail[group_cols], drop = TRUE))
  out <- do.call(rbind, lapply(parts, function(d) {
    key <- d[1, group_cols, drop = FALSE]
    vals <- unlist(lapply(metrics, function(m) c(mean = mean(d[[m]], na.rm = TRUE), se = se(d[[m]]))))
    names(vals) <- paste0(rep(metrics, each = 2), rep(c("_mean", "_se"), length(metrics)))
    data.frame(key, as.list(vals), n_replications = length(unique(d$replication_id)), mode = d$mode[1], stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  write.csv(out, outfile, row.names = FALSE)
  out
}

parse_file_meta <- function(path) {
  nm <- basename(path)
  m <- regexec("^(.*)_rho([0-9]p[0-9]+)_rep([0-9]+)_unrestricted_bma_bNA_lNA\\.rds$", nm)
  z <- regmatches(nm, m)[[1]]
  data.frame(
    scenario = z[2],
    rho = as.numeric(gsub("p", ".", z[3], fixed = TRUE)),
    replication_id = as.integer(z[4]),
    stringsAsFactors = FALSE
  )
}

run_full_draw_adaptive <- function() {
  files <- sort(Sys.glob(file.path("sim", "output", "direct_full", "*_unrestricted_bma_bNA_lNA.rds")))
  if (!length(files)) {
    return(data.frame(status = "missing_direct_full_unrestricted_rds"))
  }
  if (mode %in% c("smoke", "medium")) {
    meta_all <- do.call(rbind, lapply(files, parse_file_meta))
    meta_all$file <- files
    per_cell <- if (mode == "smoke") 1L else 2L
    files <- unlist(lapply(split(meta_all, interaction(meta_all$scenario, meta_all$rho, drop = TRUE)), function(z) head(z$file, per_cell)))
  }
  rows <- list()
  for (i in seq_along(files)) {
    fit <- readRDS(files[[i]])
    meta <- parse_file_meta(files[[i]])
    weights <- normalize_weights(fit$posterior)
    pool <- make_default_kernel_pool(
      supports = fit$supports,
      weights = weights,
      group_id = fit$group_id,
      X = NULL,
      include_intervals = FALSE,
      include_graph = FALSE,
      include_clusters = TRUE,
      include_active_sets = TRUE,
      mode = if (mode == "full") "medium" else mode
    )
    adapt <- learn_adaptive_kernel_dictionary(
      supports = fit$supports,
      weights = weights,
      group_id = fit$group_id,
      candidate_kernels = pool,
      beta = 0.02,
      max_iter = if (mode == "full") 8L else 6L,
      eta = 1e-4,
      mode = "smoke"
    )
    fixed_active <- candidate_active_sets_from_supports(fit$supports, fit$group_id, weights, max_sets = 64L)
    fixed_dict <- make_representative_dictionary(fit$group_id, fixed_active, capacity = 1L, include_safety = TRUE)
    W_fixed <- family_membership_matrix(fit$supports, fixed_dict)
    a_fixed <- estimate_family_alpha(W_fixed, weights, alpha_floor = 1e-8)$alpha_truncated
    c_fixed <- fixed_dict$families$cost
    q_fixed_fit <- optimize_family_mixture(W_fixed, a_fixed, c_fixed, weights, beta = 0.02, tau = 1e-3, distortion = "fkl", max_iter = 35L, q0_min = 1e-3)
    q_fixed <- q_fixed_fit$q
    d_adapt <- adapt$summary
    d_fixed <- data.frame(method = "fixed_hard_dictionary", t(unlist(mixture_distortions(q_fixed, W_fixed, a_fixed, weights)[c("tv", "kl_base_to_compressed", "kl_compressed_to_base")])))
    names(d_fixed)[2:4] <- c("tv", "fkl", "rkl")
    d_fixed$q0 <- q_fixed[1]
    d_fixed$expected_code <- sum(q_fixed * c_fixed)
    d_fixed$q_effective_kernels <- exp(-sum(ifelse(q_fixed > 0, q_fixed * log(q_fixed), 0)))
    d_fixed$q_kernel95 <- which(cumsum(sort(q_fixed, decreasing = TRUE)) >= 0.95)[1]
    d_fixed$n_kernels <- ncol(W_fixed)
    d_fixed$kkt_residual <- q_fixed_fit$kkt_residual
    common <- c("method", "tv", "fkl", "rkl", "q0", "expected_code", "q_effective_kernels", "q_kernel95", "n_kernels", "kkt_residual")
    z <- rbind(d_fixed[, common], d_adapt[, common])
    z <- cbind(experiment = "saved_p100_unrestricted_draws", meta, z, mode = mode)
    rows[[i]] <- z
  }
  detail <- do.call(rbind, rows)
  write.csv(detail, file.path(table_dir, "table_adaptive_kernel_full_draws_detail.csv"), row.names = FALSE)
  aggregate_summary(detail, c("scenario", "rho", "method"), file.path(table_dir, "table_adaptive_kernel_full_draws.csv"))
}

load_real_x <- function(dataset) {
  if (dataset == "tecator") {
    if (!requireNamespace("caret", quietly = TRUE)) return(NULL)
    data("tecator", package = "caret")
    as.matrix(absorp)
  } else if (dataset == "gasoline") {
    if (!requireNamespace("pls", quietly = TRUE)) return(NULL)
    data("gasoline", package = "pls")
    as.matrix(gasoline$NIR)
  } else {
    NULL
  }
}

run_semisynthetic_realx <- function() {
  X_raw <- load_real_x("tecator")
  if (is.null(X_raw) || mode == "smoke") {
    out <- data.frame(status = if (is.null(X_raw)) "tecator_unavailable" else "skipped_in_smoke")
    write.csv(out, file.path(table_dir, "table_adaptive_kernel_semisynthetic_realx.csv"), row.names = FALSE)
    return(out)
  }
  n <- nrow(X_raw)
  p0 <- if (mode == "full") 15L else 12L
  cols <- unique(round(seq(max(1, floor(ncol(X_raw) * 0.25)), min(ncol(X_raw), ceiling(ncol(X_raw) * 0.75)), length.out = p0)))
  X0 <- X_raw[, cols, drop = FALSE]
  rows <- list()
  for (rep_id in seq_len(rep_count)) {
    set.seed(5000 + rep_id)
    test_id <- sample(seq_len(n), max(10L, floor(0.25 * n)))
    train_id <- setdiff(seq_len(n), test_id)
    z <- standardize_with_training(X0[train_id, , drop = FALSE], X0[test_id, , drop = FALSE])
    beta <- rep(0, p0)
    band <- if (p0 >= 15L) 6:9 else 5:7
    beta[band] <- seq(0.8, 1.0, length.out = length(band))
    ytr <- as.numeric(z$X_train %*% beta + rnorm(nrow(z$X_train), sd = 1.25))
    yte <- as.numeric(z$X_test %*% beta + rnorm(nrow(z$X_test), sd = 1.25))
    center <- mean(ytr)
    dat <- list(X_train = z$X_train, X_test = z$X_test, y_train = ytr - center, y_test = yte - center, group_id = make_group_id(as.integer(p0 / 3L), 3L))
    fit <- fit_exact_bma(dat$X_train, dat$y_train, theta = 0.10, tau2 = 4, a0 = 1, b0 = 1)
    ref <- predict_mixture(fit, dat$X_train, dat$y_train, dat$X_test, dat$y_test, weights = fit$posterior)
    fixed <- fit_fixed_hard(dat, fit, capacity = 2L)
    pool <- make_default_kernel_pool(fit$supports, fit$posterior, dat$group_id, dat$X_train, mode = mode)
    adapt <- learn_adaptive_kernel_dictionary(fit$supports, fit$posterior, dat$group_id, dat$X_train, pool, beta = 0.02, max_iter = 10L, mode = mode)
    topm <- fit_topm(fit, top_m = 32L)
    rows[[rep_id]] <- rbind(
      evaluate_kernel_fit("fixed_hard_dictionary", fit, dat, fixed$W, fixed$alpha, fixed$q, fixed$costs, ref, list(n_kernels = fixed$n_kernels, kkt_residual = fixed$kkt_residual, active_kernels_001 = NA_real_)),
      evaluate_kernel_fit("adaptive_support_kernel", fit, dat, adapt$W, adapt$alpha, adapt$q, kernel_costs(adapt$dictionary), ref, list(n_kernels = adapt$summary$n_kernels, kkt_residual = adapt$summary$kkt_residual, active_kernels_001 = adapt$summary$active_kernels_001)),
      evaluate_kernel_fit("topM_support_atoms", fit, dat, topm$W, topm$alpha, topm$q, topm$costs, ref, list(n_kernels = topm$n_kernels, kkt_residual = topm$kkt_residual, active_kernels_001 = NA_real_))
    )
    rows[[rep_id]]$dataset <- "tecator_semisynthetic"
    rows[[rep_id]]$replication_id <- rep_id
    rows[[rep_id]]$mode <- mode
  }
  detail <- do.call(rbind, rows)
  write.csv(detail, file.path(table_dir, "table_adaptive_kernel_semisynthetic_realx_detail.csv"), row.names = FALSE)
  aggregate_summary(detail, c("dataset", "method"), file.path(table_dir, "table_adaptive_kernel_semisynthetic_realx.csv"))
}

make_figures <- function() {
  exact <- read.csv(file.path(table_dir, "table_adaptive_kernel_exact_detail.csv"))
  method_labels <- c(
    adaptive_support_kernel = "Adaptive kernels",
    fixed_hard_dictionary = "Fixed hard dictionary",
    topM_support_atoms = "Top-M atoms"
  )
  method_cols <- c(
    adaptive_support_kernel = "#1F77B4",
    fixed_hard_dictionary = "#D95F02",
    topM_support_atoms = "#2CA02C"
  )
  pdf(file.path(fig_dir, "fig_adaptive_kernel_rate_distortion.pdf"), width = 6.7, height = 4.6)
  op <- par(mar = c(4.5, 4.8, 2.0, 1.0), mgp = c(2.7, 0.8, 0))
  plot(
    exact$expected_code, exact$fkl,
    pch = 19,
    col = method_cols[exact$method],
    xlab = "Reporting cost",
    ylab = "Forward KL",
    main = "Fixed versus adaptive support-kernel compression",
    cex = 0.95,
    cex.main = 0.95,
    cex.lab = 0.95
  )
  grid(col = "grey88", lty = "dotted")
  present <- names(method_labels)[names(method_labels) %in% unique(exact$method)]
  legend(
    "topleft",
    inset = c(0.035, 0.035),
    legend = unname(method_labels[present]),
    col = method_cols[present],
    pch = 19,
    bty = "o",
    bg = grDevices::adjustcolor("white", alpha.f = 0.9),
    box.col = "grey85",
    cex = 0.82
  )
  par(op)
  dev.off()
  if (file.exists(file.path(table_dir, "table_adaptive_kernel_semisynthetic_realx_detail.csv"))) {
    realx <- read.csv(file.path(table_dir, "table_adaptive_kernel_semisynthetic_realx_detail.csv"))
    if ("method" %in% names(realx)) {
      pdf(file.path(fig_dir, "fig_adaptive_kernel_semisynthetic_realx.pdf"), width = 6.7, height = 4.5)
      op <- par(mar = c(4.5, 4.8, 2.0, 1.0), mgp = c(2.7, 0.8, 0))
      plot(
        realx$expected_code, realx$fkl,
        pch = 19,
        col = method_cols[realx$method],
        xlab = "Reporting cost",
        ylab = "Forward KL",
        main = "Semi-synthetic Tecator rate-distortion",
        cex = 1.05,
        cex.main = 0.95,
        cex.lab = 0.95
      )
      grid(col = "grey88", lty = "dotted")
      present <- names(method_labels)[names(method_labels) %in% unique(realx$method)]
      legend(
        "topleft",
        inset = c(0.035, 0.035),
        legend = unname(method_labels[present]),
        col = method_cols[present],
        pch = 19,
        bty = "o",
        bg = grDevices::adjustcolor("white", alpha.f = 0.9),
        box.col = "grey85",
        cex = 0.82
      )
      par(op)
      dev.off()
    }
  }
}

invisible(run_exact_identity())
invisible(run_topm_separation())
invisible(run_exact_adaptive())
invisible(run_full_draw_adaptive())
invisible(run_semisynthetic_realx())
invisible(make_figures())

cat("Adaptive support-kernel experiments complete in", mode, "mode\n")
