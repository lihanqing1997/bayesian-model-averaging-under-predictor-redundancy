source(file.path("sim", "design.R"))
source(file.path("sim", "exact_bma.R"))
source(file.path("sim", "comparison_methods.R"))
source(file.path("sim", "src", "family_dictionary.R"))
source(file.path("sim", "src", "adaptive_support_kernel_compression.R"))

benchmark_mode_config <- function(mode = c("smoke", "pilot", "medium", "full")) {
  mode <- match.arg(mode)
  list(
    mode = mode,
    rep_count = switch(mode, smoke = 1L, pilot = 2L, medium = 3L, full = 5L),
    rho_grid = switch(mode, smoke = 0.9, pilot = 0.9, medium = c(0.7, 0.9), full = c(0.7, 0.9)),
    scenario_grid = switch(
      mode,
      smoke = c("one_representative"),
      pilot = c("one_representative", "weak_signal"),
      medium = c("one_representative", "multi_representative", "weak_signal", "noisy_grouping"),
      full = c(
        "one_representative",
        "multi_representative",
        "weak_signal",
        "noisy_grouping",
        "misspecified_grouping",
        "ordered_interval",
        "graph_community"
      )
    ),
    beta_grid = switch(mode, smoke = 0.02, pilot = 0.02, medium = c(0.005, 0.02, 0.08), full = c(0.005, 0.02, 0.08)),
    topm_grid = switch(mode, smoke = 16L, pilot = c(16L, 32L), medium = c(8L, 16L, 32L), full = c(8L, 16L, 32L)),
    credible_grid = switch(mode, smoke = 0.90, pilot = c(0.90, 0.95), medium = c(0.80, 0.90, 0.95), full = c(0.80, 0.90, 0.95)),
    cluster_grid = switch(mode, smoke = 6L, pilot = 8L, medium = c(4L, 8L, 12L), full = c(4L, 8L, 12L)),
    dilution_power_grid = switch(mode, smoke = 0.5, pilot = 0.5, medium = c(0.25, 0.5, 1), full = c(0.25, 0.5, 1)),
    dpp_power_grid = switch(mode, smoke = 1, pilot = 1, medium = c(0.5, 1), full = c(0.5, 1))
  )
}

simulate_exact_benchmark_data <- function(scenario,
                                          rho,
                                          seed,
                                          mode = "smoke") {
  set.seed(seed)
  K <- 4L
  m <- 3L
  n_train <- if (identical(mode, "full")) 80L else 70L
  n_test <- if (identical(mode, "full")) 260L else 220L
  p <- K * m
  group_id <- make_group_id(K, m)
  X_train_raw <- simulate_block_design(n_train, K, m, rho)
  X_test_raw <- simulate_block_design(n_test, K, m, rho)

  if (identical(scenario, "ordered_interval")) {
    p <- 12L
    K <- p
    m <- 1L
    group_id <- seq_len(p)
    z_train <- as.numeric(stats::filter(rnorm(n_train + p), filter = rep(1, 4) / 4, sides = 1))
    z_test <- as.numeric(stats::filter(rnorm(n_test + p), filter = rep(1, 4) / 4, sides = 1))
    z_train[!is.finite(z_train)] <- rnorm(sum(!is.finite(z_train)))
    z_test[!is.finite(z_test)] <- rnorm(sum(!is.finite(z_test)))
    X_train_raw <- matrix(0, n_train, p)
    X_test_raw <- matrix(0, n_test, p)
    for (j in seq_len(p)) {
      X_train_raw[, j] <- sqrt(rho) * z_train[j:(j + n_train - 1L)] + sqrt(1 - rho) * rnorm(n_train)
      X_test_raw[, j] <- sqrt(rho) * z_test[j:(j + n_test - 1L)] + sqrt(1 - rho) * rnorm(n_test)
    }
  }

  z <- standardize_with_training(X_train_raw, X_test_raw)
  beta <- rep(0, ncol(z$X_train))
  sigma <- 1.35

  if (scenario %in% c("one_representative", "noisy_grouping", "misspecified_grouping")) {
    idx <- c(1, 7, min(11, length(beta)))
    beta[idx] <- c(1.4, 1.1, 0.9)
  } else if (identical(scenario, "multi_representative")) {
    idx <- c(1, 2, 7, min(11, length(beta)))
    beta[idx] <- c(0.9, 0.75, 1.1, 0.9)
  } else if (identical(scenario, "weak_signal")) {
    idx <- c(1, 7, min(11, length(beta)))
    beta[idx] <- c(0.75, 0.65, 0.55)
    sigma <- 1.8
  } else if (identical(scenario, "ordered_interval")) {
    band <- 5:7
    beta[band] <- c(0.8, 1.1, 0.8)
  } else if (identical(scenario, "graph_community")) {
    idx <- c(1, 2, 6, 7, min(11, length(beta)))
    beta[idx] <- c(0.65, 0.55, 0.85, 0.65, 0.7)
  }

  if (identical(scenario, "noisy_grouping")) {
    swap <- sample(seq_along(group_id), size = max(1L, floor(0.25 * length(group_id))))
    group_id[swap] <- sample(group_id[swap])
  }
  if (identical(scenario, "misspecified_grouping")) {
    group_id <- sample(group_id)
  }

  y_train_raw <- as.numeric(z$X_train %*% beta + rnorm(nrow(z$X_train), sd = sigma))
  y_test_raw <- as.numeric(z$X_test %*% beta + rnorm(nrow(z$X_test), sd = sigma))
  center <- mean(y_train_raw)
  active_groups <- sort(unique(group_id[beta != 0]))

  list(
    X_train = z$X_train,
    X_test = z$X_test,
    y_train = y_train_raw - center,
    y_test = y_test_raw - center,
    beta = beta,
    true_support = as.integer(beta != 0),
    group_id = group_id,
    active_groups = active_groups,
    K = max(group_id),
    m = m,
    rho = rho,
    scenario = scenario
  )
}

posterior_distribution_metrics <- function(reference_weights,
                                           approx_weights,
                                           h_floor = 1e-12) {
  p0 <- normalize_weights(reference_weights)
  p1 <- normalize_weights(approx_weights)
  h <- p1 / p0
  data.frame(
    tv = 0.5 * sum(abs(p1 - p0)),
    fkl = sum(p0 * (-log(pmax(h, h_floor)))),
    rkl = sum(ifelse(p1 > 0, p1 * log(pmax(h, h_floor)), 0)),
    mean_h_error = abs(sum(p0 * h) - 1),
    stringsAsFactors = FALSE
  )
}

mixture_weights_from_kernel_fit <- function(reference_fit, W, alpha, q) {
  h <- mixture_h(q, W, alpha)
  out <- reference_fit$posterior * h
  normalize_weights(out)
}

evaluate_weight_summary <- function(method,
                                    reference_fit,
                                    approx_weights,
                                    dat,
                                    reference_pred,
                                    expected_code = NA_real_,
                                    q0 = NA_real_,
                                    active_kernels = NA_real_,
                                    effective_kernels = NA_real_,
                                    n_kernels = NA_real_,
                                    stored_atoms = NA_real_,
                                    objective = NA_real_,
                                    kkt_residual = NA_real_,
                                    optimizer_status = NA_character_,
                                    optimizer_iterations = NA_real_,
                                    objective_change = NA_real_,
                                    raw_candidates = NA_real_,
                                    key_deduped_candidates = NA_real_,
                                    near_deduped_candidates = NA_real_,
                                    objective_consistent = NA,
                                    runtime_sec = NA_real_,
                                    beta = NA_real_,
                                    tuning = NA_character_,
                                    method_family = NA_character_) {
  d <- posterior_distribution_metrics(reference_fit$posterior, approx_weights)
  pred <- predict_mixture(
    reference_fit,
    X_train = dat$X_train,
    y_train = dat$y_train,
    X_test = dat$X_test,
    y_test = dat$y_test,
    weights = approx_weights
  )
  data.frame(
    method = method,
    method_family = method_family,
    tuning = tuning,
    beta = beta,
    tv = d$tv,
    fkl = d$fkl,
    rkl = d$rkl,
    mean_h_error = d$mean_h_error,
    q0 = q0,
    expected_code = expected_code,
    active_kernels_001 = active_kernels,
    q_effective_kernels = effective_kernels,
    n_kernels = n_kernels,
    stored_atoms = stored_atoms,
    objective = objective,
    rmse = pred$rmse,
    logscore = pred$mean_log_score,
    rmse_gap = pred$rmse - reference_pred$rmse,
    logscore_gap = pred$mean_log_score - reference_pred$mean_log_score,
    kkt_residual = kkt_residual,
    optimizer_status = optimizer_status,
    optimizer_iterations = optimizer_iterations,
    objective_change = objective_change,
    raw_candidates = raw_candidates,
    key_deduped_candidates = key_deduped_candidates,
    near_deduped_candidates = near_deduped_candidates,
    objective_consistent = objective_consistent,
    runtime_sec = runtime_sec,
    status = "ok",
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

evaluate_predictive_summary <- function(method,
                                        dat,
                                        reference_pred,
                                        pred,
                                        selected_count,
                                        runtime_sec = NA_real_,
                                        tuning = NA_character_,
                                        method_family = "predictive selection") {
  data.frame(
    method = method,
    method_family = method_family,
    tuning = tuning,
    beta = NA_real_,
    tv = NA_real_,
    fkl = NA_real_,
    rkl = NA_real_,
    mean_h_error = NA_real_,
    q0 = NA_real_,
    expected_code = selected_count,
    active_kernels_001 = NA_real_,
    q_effective_kernels = NA_real_,
    n_kernels = NA_real_,
    stored_atoms = selected_count,
    objective = NA_real_,
    rmse = pred$rmse,
    logscore = pred$mean_log_score,
    rmse_gap = pred$rmse - reference_pred$rmse,
    logscore_gap = pred$mean_log_score - reference_pred$mean_log_score,
    kkt_residual = NA_real_,
    optimizer_status = NA_character_,
    optimizer_iterations = NA_real_,
    objective_change = NA_real_,
    raw_candidates = NA_real_,
    key_deduped_candidates = NA_real_,
    near_deduped_candidates = NA_real_,
    objective_consistent = NA,
    runtime_sec = runtime_sec,
    status = "ok",
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

failure_summary <- function(method,
                            method_family = NA_character_,
                            tuning = NA_character_,
                            beta = NA_real_,
                            error) {
  data.frame(
    method = method,
    method_family = method_family,
    tuning = tuning,
    beta = beta,
    tv = NA_real_,
    fkl = NA_real_,
    rkl = NA_real_,
    mean_h_error = NA_real_,
    q0 = NA_real_,
    expected_code = NA_real_,
    active_kernels_001 = NA_real_,
    q_effective_kernels = NA_real_,
    n_kernels = NA_real_,
    stored_atoms = NA_real_,
    objective = NA_real_,
    rmse = NA_real_,
    logscore = NA_real_,
    rmse_gap = NA_real_,
    logscore_gap = NA_real_,
    kkt_residual = NA_real_,
    optimizer_status = NA_character_,
    optimizer_iterations = NA_real_,
    objective_change = NA_real_,
    raw_candidates = NA_real_,
    key_deduped_candidates = NA_real_,
    near_deduped_candidates = NA_real_,
    objective_consistent = NA,
    runtime_sec = NA_real_,
    status = "failed",
    error = conditionMessage(error),
    stringsAsFactors = FALSE
  )
}

cell_failure_result <- function(scenario, rho, replication_id, config, error) {
  out <- failure_summary(
    method = "Cell failure",
    method_family = "reference posterior or simulation",
    tuning = NA_character_,
    beta = NA_real_,
    error = error
  )
  out$scenario <- scenario
  out$rho <- rho
  out$replication_id <- replication_id
  out$mode <- config$mode
  out$reference_support_entropy <- NA_real_
  out$reference_support_n95 <- NA_real_
  out
}

evaluate_kernel_summary <- function(method,
                                    reference_fit,
                                    dat,
                                    W,
                                    alpha,
                                    q,
                                    costs,
                                    reference_pred,
                                    beta = NA_real_,
                                    tuning = NA_character_,
                                    method_family = NA_character_,
                                    objective = NA_real_,
                                    kkt_residual = NA_real_,
                                    optimizer_status = NA_character_,
                                    optimizer_iterations = NA_real_,
                                    objective_change = NA_real_,
                                    raw_candidates = NA_real_,
                                    key_deduped_candidates = NA_real_,
                                    near_deduped_candidates = NA_real_,
                                    objective_consistent = NA) {
  q <- normalize_weights(q)
  approx_weights <- mixture_weights_from_kernel_fit(reference_fit, W, alpha, q)
  entropy <- -sum(ifelse(q > 0, q * log(q), 0))
  evaluate_weight_summary(
    method = method,
    reference_fit = reference_fit,
    approx_weights = approx_weights,
    dat = dat,
    reference_pred = reference_pred,
    expected_code = sum(q * costs),
    q0 = q[1],
    active_kernels = sum(q > 0.001),
    effective_kernels = exp(entropy),
    n_kernels = ncol(W),
    stored_atoms = NA_real_,
    objective = objective,
    kkt_residual = kkt_residual,
    optimizer_status = optimizer_status,
    optimizer_iterations = optimizer_iterations,
    objective_change = objective_change,
    raw_candidates = raw_candidates,
    key_deduped_candidates = key_deduped_candidates,
    near_deduped_candidates = near_deduped_candidates,
    objective_consistent = objective_consistent,
    beta = beta,
    tuning = tuning,
    method_family = method_family
  )
}

fit_fixed_hard_benchmark <- function(dat,
                                     reference_fit,
                                     beta = 0.02,
                                     q0_max = 1,
                                     capacity = 1L,
                                     max_sets = 48L) {
  active_sets <- candidate_active_sets_from_supports(
    supports = reference_fit$supports,
    group_id = dat$group_id,
    weights = reference_fit$posterior,
    max_sets = max_sets
  )
  dict <- make_representative_dictionary(dat$group_id, active_sets, capacity = capacity, include_safety = TRUE)
  W <- family_membership_matrix(reference_fit$supports, dict)
  alpha <- estimate_family_alpha(W, reference_fit$posterior, alpha_floor = 1e-8)$alpha_truncated
  costs <- dict$families$cost
  qfit <- optimize_family_mixture(
    W,
    alpha,
    costs,
    reference_fit$posterior,
    beta = beta,
    tau = 1e-3,
    distortion = "fkl",
    q0_min = 1e-3,
    q0_max = q0_max,
    max_iter = 120L
  )
  list(W = W, alpha = alpha, costs = costs, q = qfit$q, fit = qfit, n_kernels = ncol(W))
}

fit_topm_benchmark <- function(reference_fit,
                               beta = 0.02,
                               q0_max = 1,
                               top_m = 32L) {
  dict <- make_support_atom_dictionary(reference_fit$supports, reference_fit$posterior, top_m = top_m, include_safety = TRUE)
  W <- family_membership_matrix(reference_fit$supports, dict)
  alpha <- estimate_family_alpha(W, reference_fit$posterior, alpha_floor = 1e-8)$alpha_truncated
  costs <- dict$families$cost
  qfit <- optimize_family_mixture(
    W,
    alpha,
    costs,
    reference_fit$posterior,
    beta = beta,
    tau = 1e-3,
    distortion = "fkl",
    q0_min = 1e-3,
    q0_max = q0_max,
    max_iter = 120L
  )
  list(W = W, alpha = alpha, costs = costs, q = qfit$q, fit = qfit, n_kernels = ncol(W), stored_atoms = top_m)
}

fit_credible_support_benchmark <- function(reference_fit,
                                           beta = 0.02,
                                           q0_max = 1,
                                           coverage = 0.95) {
  ord <- order(reference_fit$posterior, decreasing = TRUE)
  cum_mass <- cumsum(reference_fit$posterior[ord])
  top_m <- which(cum_mass >= coverage)[1]
  if (!is.finite(top_m)) {
    top_m <- length(ord)
  }
  fit <- fit_topm_benchmark(reference_fit, beta = beta, q0_max = q0_max, top_m = top_m)
  fit$coverage <- coverage
  fit$stored_atoms <- top_m
  fit
}

fit_posterior_clustering_benchmark <- function(dat,
                                               reference_fit,
                                               beta = 0.02,
                                               q0_max = 1,
                                               n_clusters = 8L,
                                               n_draws = 1500L,
                                               seed = 1L) {
  set.seed(seed)
  weights <- normalize_weights(reference_fit$posterior)
  draw_id <- sample(seq_along(weights), size = n_draws, replace = TRUE, prob = weights)
  Xs <- reference_fit$supports[draw_id, , drop = FALSE]
  k <- min(as.integer(n_clusters), nrow(unique(Xs)))
  if (k <= 1L) {
    centers <- matrix(round(colMeans(Xs) > 0.5), nrow = 1)
  } else {
    km <- stats::kmeans(Xs, centers = k, iter.max = 50)
    centers <- (km$centers >= 0.5) * 1L
  }
  centers <- unique(as.matrix(centers))
  kernels <- list(safety_kernel(cost = max(10, 4 * log(max(ncol(reference_fit$supports), 2)))))
  for (i in seq_len(nrow(centers))) {
    kernels[[length(kernels) + 1L]] <- soft_hamming_kernel(centers[i, ], rho = 1, name = paste0("cluster_atom_", i))
    if (!is.null(dat$group_id)) {
      cnt <- support_group_counts(centers[i, , drop = FALSE], dat$group_id)
      kernels[[length(kernels) + 1L]] <- soft_group_hamming_kernel(which(cnt[1, ] > 0), dat$group_id, rho = 1, name = paste0("cluster_group_", i))
    }
  }
  dict <- new_support_kernel_dictionary(dedupe_kernels(kernels))
  W <- support_kernel_weight_matrix(reference_fit$supports, dict, group_id = dat$group_id)
  alpha <- estimate_kernel_alpha(W, reference_fit$posterior, alpha_floor = 1e-8)$alpha_truncated
  costs <- kernel_costs(dict)
  qfit <- optimize_family_mixture(
    W,
    alpha,
    costs,
    reference_fit$posterior,
    beta = beta,
    tau = 1e-3,
    distortion = "fkl",
    q0_min = 1e-3,
    q0_max = q0_max,
    max_iter = 120L
  )
  list(W = W, alpha = alpha, costs = costs, q = qfit$q, fit = qfit, n_kernels = ncol(W), centers = centers)
}

fit_dilution_benchmark <- function(dat,
                                   reference_fit,
                                   power = 0.5) {
  log_det <- support_logdet_correlation(reference_fit$supports, dat$X_train)
  log_prior <- log_prior_bernoulli(reference_fit$supports, reference_fit$theta) + power * log_det
  fit_bma_from_log_prior(
    supports = reference_fit$supports,
    log_marginal = reference_fit$log_marginal,
    log_prior = log_prior,
    theta = reference_fit$theta,
    tau2 = reference_fit$tau2,
    a0 = reference_fit$a0,
    b0 = reference_fit$b0,
    n = reference_fit$n
  )
}

fit_dpp_benchmark <- function(dat,
                              reference_fit,
                              power = 1) {
  log_det <- support_logdet_correlation(reference_fit$supports, dat$X_train)
  log_prior <- log_prior_bernoulli(reference_fit$supports, reference_fit$theta) + power * log_det
  fit_bma_from_log_prior(
    supports = reference_fit$supports,
    log_marginal = reference_fit$log_marginal,
    log_prior = log_prior,
    theta = reference_fit$theta,
    tau2 = reference_fit$tau2,
    a0 = reference_fit$a0,
    b0 = reference_fit$b0,
    n = reference_fit$n
  )
}

fit_glmnet_predictive_benchmark <- function(dat,
                                            alpha,
                                            method,
                                            seed) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet is required for predictive-selection baselines")
  }
  set.seed(seed)
  elapsed <- system.time({
    fit <- glmnet::cv.glmnet(
      x = dat$X_train,
      y = dat$y_train,
      family = "gaussian",
      alpha = alpha,
      nfolds = 5,
      intercept = TRUE,
      standardize = FALSE
    )
    pred_train <- as.numeric(stats::predict(fit, newx = dat$X_train, s = "lambda.min"))
    pred_test <- as.numeric(stats::predict(fit, newx = dat$X_test, s = "lambda.min"))
  })[["elapsed"]]
  resid <- dat$y_train - pred_train
  sigma <- sqrt(mean(resid^2))
  sigma <- max(sigma, 1e-6)
  coef_vec <- as.matrix(stats::coef(fit, s = "lambda.min"))[-1, 1]
  selected_count <- sum(abs(coef_vec) > 1e-8)
  pred <- list(
    rmse = sqrt(mean((dat$y_test - pred_test)^2)),
    mean_log_score = mean(stats::dnorm(dat$y_test, mean = pred_test, sd = sigma, log = TRUE))
  )
  list(
    method = method,
    pred = pred,
    selected_count = selected_count,
    runtime_sec = elapsed,
    tuning = paste0("alpha=", alpha, ", lambda.min")
  )
}

fit_group_lasso_predictive_benchmark <- function(dat, seed) {
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("grpreg is required for group-lasso predictive baseline")
  }
  if (is.null(dat$group_id)) {
    stop("group_id is required for group lasso")
  }
  set.seed(seed)
  elapsed <- system.time({
    fit <- grpreg::cv.grpreg(
      X = dat$X_train,
      y = dat$y_train,
      group = dat$group_id,
      family = "gaussian",
      penalty = "grLasso",
      nfolds = 5
    )
    pred_train <- as.numeric(stats::predict(fit, X = dat$X_train, type = "response"))
    pred_test <- as.numeric(stats::predict(fit, X = dat$X_test, type = "response"))
  })[["elapsed"]]
  resid <- dat$y_train - pred_train
  sigma <- sqrt(mean(resid^2))
  sigma <- max(sigma, 1e-6)
  coef_vec <- as.numeric(stats::coef(fit))[-1]
  selected_vars <- which(abs(coef_vec) > 1e-8)
  selected_count <- length(selected_vars)
  selected_groups <- if (selected_count) length(unique(dat$group_id[selected_vars])) else 0L
  pred <- list(
    rmse = sqrt(mean((dat$y_test - pred_test)^2)),
    mean_log_score = mean(stats::dnorm(dat$y_test, mean = pred_test, sd = sigma, log = TRUE))
  )
  list(
    method = "Group lasso",
    pred = pred,
    selected_count = selected_groups,
    runtime_sec = elapsed,
    tuning = "group lasso, lambda.min"
  )
}

run_one_exact_benchmark_cell <- function(scenario,
                                         rho,
                                         replication_id,
                                         config) {
  seed <- 20260527 + replication_id + round(1000 * rho) + 10000 * match(scenario, config$scenario_grid)
  dat <- simulate_exact_benchmark_data(scenario, rho, seed, mode = config$mode)
  reference_fit <- fit_exact_bma(dat$X_train, dat$y_train, theta = 0.08, tau2 = 4, a0 = 1, b0 = 1)
  reference_pred <- predict_mixture(
    reference_fit,
    dat$X_train,
    dat$y_train,
    dat$X_test,
    dat$y_test,
    weights = reference_fit$posterior
  )
  rows <- list()
  add_row <- function(x) {
    rows[[length(rows) + 1L]] <<- x
  }
  add_failure <- function(method, method_family, tuning, beta, error) {
    add_row(failure_summary(method, method_family, tuning, beta, error))
  }

  for (beta in config$beta_grid) {
    ask_mode <- if (identical(config$mode, "full")) "medium" else config$mode
    tuning_beta <- paste0("beta=", signif(beta, 3))
    ask <- tryCatch(
      fit_askpc_pooled_pruned(
        supports = reference_fit$supports,
        weights = reference_fit$posterior,
        group_id = dat$group_id,
        X = dat$X_train,
        beta = beta,
        tau = 1e-3,
        q0_min = 1e-3,
        prune_mass = 0.99,
        mode = ask_mode,
        max_iter = if (identical(config$mode, "smoke")) 80L else if (identical(config$mode, "full")) 120L else NULL
      ),
      error = function(e) e
    )
    if (inherits(ask, "error")) {
      add_failure("ASK-PC pooled union", "support-kernel compression", tuning_beta, beta, ask)
      add_failure("ASK-PC pooled-pruned 99%", "support-kernel compression", tuning_beta, beta, ask)
    } else {
      add_row(evaluate_kernel_summary(
        method = "ASK-PC pooled union",
        reference_fit = reference_fit,
        dat = dat,
        W = ask$pooled$W,
        alpha = ask$pooled$alpha,
        q = ask$pooled$fit$q,
        costs = ask$pooled$costs,
        reference_pred = reference_pred,
        beta = beta,
        tuning = tuning_beta,
        method_family = "support-kernel compression",
        objective = ask$pooled$fit$objective,
        kkt_residual = ask$pooled$fit$kkt_residual,
        optimizer_status = ask$pooled$fit$status,
        optimizer_iterations = ask$pooled$fit$iterations,
        objective_change = ask$pooled$fit$objective_change,
        raw_candidates = ask$pool_diagnostics$raw_candidates,
        key_deduped_candidates = ask$pool_diagnostics$key_deduped_candidates,
        near_deduped_candidates = ask$pool_diagnostics$near_deduped_candidates,
        objective_consistent = ask$pool_diagnostics$objective_consistent
      ))
      add_row(evaluate_kernel_summary(
        method = "ASK-PC pooled-pruned 99%",
        reference_fit = reference_fit,
        dat = dat,
        W = ask$pruned$W,
        alpha = ask$pruned$alpha,
        q = ask$pruned$fit$q,
        costs = ask$pruned$costs,
        reference_pred = reference_pred,
        beta = beta,
        tuning = tuning_beta,
        method_family = "support-kernel compression",
        objective = ask$pruned$fit$objective,
        kkt_residual = ask$pruned$fit$kkt_residual,
        optimizer_status = ask$pruned$fit$status,
        optimizer_iterations = ask$pruned$fit$iterations,
        objective_change = ask$pruned$fit$objective_change,
        raw_candidates = ask$pool_diagnostics$raw_candidates,
        key_deduped_candidates = ask$pool_diagnostics$key_deduped_candidates,
        near_deduped_candidates = ask$pool_diagnostics$near_deduped_candidates,
        objective_consistent = ask$pool_diagnostics$objective_consistent
      ))
    }

    capacity <- if (identical(scenario, "multi_representative")) 2L else 1L
    fixed <- tryCatch(fit_fixed_hard_benchmark(dat, reference_fit, beta = beta, capacity = capacity), error = function(e) e)
    if (inherits(fixed, "error")) {
      add_failure("Fixed hard dictionary", "fixed region dictionary", tuning_beta, beta, fixed)
    } else {
      add_row(evaluate_kernel_summary(
        method = "Fixed hard dictionary",
        reference_fit = reference_fit,
        dat = dat,
        W = fixed$W,
        alpha = fixed$alpha,
        q = fixed$q,
        costs = fixed$costs,
        reference_pred = reference_pred,
        beta = beta,
        tuning = tuning_beta,
        method_family = "fixed region dictionary",
        objective = fixed$fit$objective,
        kkt_residual = fixed$fit$kkt_residual,
        optimizer_status = fixed$fit$status,
        optimizer_iterations = fixed$fit$iterations,
        objective_change = fixed$fit$objective_change
      ))
    }

    for (top_m in config$topm_grid) {
      tuning_topm <- paste0("M=", top_m, ", beta=", signif(beta, 3))
      topm <- tryCatch(fit_topm_benchmark(reference_fit, beta = beta, top_m = top_m), error = function(e) e)
      if (inherits(topm, "error")) {
        add_failure("Top-M support atoms", "atom truncation", tuning_topm, beta, topm)
      } else {
        z <- evaluate_kernel_summary(
          method = "Top-M support atoms",
          reference_fit = reference_fit,
          dat = dat,
          W = topm$W,
          alpha = topm$alpha,
          q = topm$q,
          costs = topm$costs,
          reference_pred = reference_pred,
          beta = beta,
          tuning = tuning_topm,
          method_family = "atom truncation",
          objective = topm$fit$objective,
          kkt_residual = topm$fit$kkt_residual,
          optimizer_status = topm$fit$status,
          optimizer_iterations = topm$fit$iterations,
          objective_change = topm$fit$objective_change
        )
        z$stored_atoms <- top_m
        add_row(z)
      }
    }

    for (coverage in config$credible_grid) {
      tuning_credible <- paste0("coverage=", signif(coverage, 3), ", beta=", signif(beta, 3))
      credible <- tryCatch(
        fit_credible_support_benchmark(reference_fit, beta = beta, coverage = coverage),
        error = function(e) e
      )
      if (inherits(credible, "error")) {
        add_failure("Credible support set", "credible support summary", tuning_credible, beta, credible)
      } else {
        z <- evaluate_kernel_summary(
          method = "Credible support set",
          reference_fit = reference_fit,
          dat = dat,
          W = credible$W,
          alpha = credible$alpha,
          q = credible$q,
          costs = credible$costs,
          reference_pred = reference_pred,
          beta = beta,
          tuning = tuning_credible,
          method_family = "credible support summary",
          objective = credible$fit$objective,
          kkt_residual = credible$fit$kkt_residual,
          optimizer_status = credible$fit$status,
          optimizer_iterations = credible$fit$iterations,
          objective_change = credible$fit$objective_change
        )
        z$stored_atoms <- credible$stored_atoms
        add_row(z)
      }
    }

    for (k in config$cluster_grid) {
      tuning_cluster <- paste0("clusters=", k, ", beta=", signif(beta, 3))
      cl <- tryCatch(
        fit_posterior_clustering_benchmark(
          dat,
          reference_fit,
          beta = beta,
          n_clusters = k,
          seed = seed + k
        ),
        error = function(e) e
      )
      if (inherits(cl, "error")) {
        add_failure("Posterior clustering", "posterior clustering", tuning_cluster, beta, cl)
      } else {
        add_row(evaluate_kernel_summary(
          method = "Posterior clustering",
          reference_fit = reference_fit,
          dat = dat,
          W = cl$W,
          alpha = cl$alpha,
          q = cl$q,
          costs = cl$costs,
          reference_pred = reference_pred,
          beta = beta,
          tuning = tuning_cluster,
          method_family = "posterior clustering",
          objective = cl$fit$objective,
          kkt_residual = cl$fit$kkt_residual,
          optimizer_status = cl$fit$status,
          optimizer_iterations = cl$fit$iterations,
          objective_change = cl$fit$objective_change
        ))
      }
    }
  }

  for (power in config$dilution_power_grid) {
    tuning_dilution <- paste0("power=", signif(power, 3))
    dil <- tryCatch(fit_dilution_benchmark(dat, reference_fit, power = power), error = function(e) e)
    if (inherits(dil, "error")) {
      add_failure("Dilution-prior BMA", "redundancy-aware prior", tuning_dilution, NA_real_, dil)
    } else {
      support_sum <- posterior_summary(dil$posterior)
      add_row(evaluate_weight_summary(
        method = "Dilution-prior BMA",
        reference_fit = reference_fit,
        approx_weights = dil$posterior,
        dat = dat,
        reference_pred = reference_pred,
        expected_code = support_sum$n_mass,
        n_kernels = NA_real_,
        stored_atoms = support_sum$n_mass,
        objective = NA_real_,
        tuning = tuning_dilution,
        method_family = "redundancy-aware prior"
      ))
    }
  }

  for (power in config$dpp_power_grid) {
    tuning_dpp <- paste0("power=", signif(power, 3))
    dpp <- tryCatch(fit_dpp_benchmark(dat, reference_fit, power = power), error = function(e) e)
    if (inherits(dpp, "error")) {
      add_failure("DPP-prior BMA", "redundancy-aware prior", tuning_dpp, NA_real_, dpp)
    } else {
      support_sum <- posterior_summary(dpp$posterior)
      add_row(evaluate_weight_summary(
        method = "DPP-prior BMA",
        reference_fit = reference_fit,
        approx_weights = dpp$posterior,
        dat = dat,
        reference_pred = reference_pred,
        expected_code = support_sum$n_mass,
        n_kernels = NA_real_,
        stored_atoms = support_sum$n_mass,
        objective = NA_real_,
        tuning = tuning_dpp,
        method_family = "redundancy-aware prior"
      ))
    }
  }

  predictive_specs <- list(
    list(method = "Lasso", alpha = 1),
    list(method = "Elastic net", alpha = 0.5)
  )
  for (spec in predictive_specs) {
    pred_fit <- tryCatch(
      fit_glmnet_predictive_benchmark(dat, alpha = spec$alpha, method = spec$method, seed = seed + round(100 * spec$alpha)),
      error = function(e) e
    )
    if (inherits(pred_fit, "error")) {
      add_failure(spec$method, "predictive selection", paste0("alpha=", spec$alpha), NA_real_, pred_fit)
    } else {
      add_row(evaluate_predictive_summary(
        method = spec$method,
        dat = dat,
        reference_pred = reference_pred,
        pred = pred_fit$pred,
        selected_count = pred_fit$selected_count,
        runtime_sec = pred_fit$runtime_sec,
        tuning = pred_fit$tuning
      ))
    }
  }
  group_pred <- tryCatch(fit_group_lasso_predictive_benchmark(dat, seed = seed + 333L), error = function(e) e)
  if (inherits(group_pred, "error")) {
    add_failure("Group lasso", "predictive selection", "group lasso", NA_real_, group_pred)
  } else {
    add_row(evaluate_predictive_summary(
      method = "Group lasso",
      dat = dat,
      reference_pred = reference_pred,
      pred = group_pred$pred,
      selected_count = group_pred$selected_count,
      runtime_sec = group_pred$runtime_sec,
      tuning = group_pred$tuning
    ))
  }

  out <- do.call(rbind, rows)
  out$scenario <- scenario
  out$rho <- rho
  out$replication_id <- replication_id
  out$mode <- config$mode
  out$reference_support_entropy <- posterior_summary(reference_fit$posterior)$entropy
  out$reference_support_n95 <- posterior_summary(reference_fit$posterior)$n_mass
  out
}

summarize_benchmark_detail <- function(detail) {
  metrics <- c(
    "tv", "fkl", "rkl", "q0", "expected_code", "active_kernels_001",
    "q_effective_kernels", "n_kernels", "stored_atoms", "objective",
    "rmse_gap", "logscore_gap", "kkt_residual", "optimizer_iterations",
    "objective_change", "raw_candidates", "key_deduped_candidates",
    "near_deduped_candidates", "runtime_sec", "reference_support_n95"
  )
  metrics <- intersect(metrics, names(detail))
  se <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  }
  groups <- c("scenario", "rho", "method", "method_family", "tuning", "mode")
  parts <- split(detail, interaction(detail[groups], drop = TRUE), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(d) {
    key <- d[1, groups, drop = FALSE]
    vals <- unlist(lapply(metrics, function(m) c(mean = mean(d[[m]], na.rm = TRUE), se = se(d[[m]]))))
    names(vals) <- paste0(rep(metrics, each = 2), rep(c("_mean", "_se"), length(metrics)))
    data.frame(key, as.list(vals), n_replications = length(unique(d$replication_id)), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

benchmark_cell_grid <- function(config) {
  expand.grid(
    scenario = config$scenario_grid,
    rho = config$rho_grid,
    replication_id = seq_len(config$rep_count),
    stringsAsFactors = FALSE
  )
}

append_benchmark_checkpoint <- function(path, rows) {
  write.table(
    rows,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(path),
    append = file.exists(path),
    quote = TRUE
  )
}

run_exact_competitor_benchmark <- function(mode = c("smoke", "pilot", "medium", "full"),
                                           table_dir = file.path("sim", "output", "tables"),
                                           checkpoint = FALSE,
                                           resume = TRUE,
                                           workers = 1L,
                                           chunk_size = 8L) {
  config <- benchmark_mode_config(mode)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_path <- file.path(
    table_dir,
    sprintf("table_support_kernel_competitor_benchmark_detail_%s_checkpoint.csv", config$mode)
  )
  cell_key <- function(scenario, rho, replication_id) {
    paste(scenario, sprintf("%.6f", rho), replication_id, sep = "|")
  }
  done <- character(0)
  if (checkpoint && !resume && file.exists(checkpoint_path)) {
    backup_path <- sub("\\.csv$", paste0("_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"), checkpoint_path)
    file.rename(checkpoint_path, backup_path)
    message(sprintf("starting fresh checkpoint; moved old checkpoint to %s", backup_path))
  }
  if (checkpoint && resume && file.exists(checkpoint_path)) {
    previous <- read.csv(checkpoint_path, stringsAsFactors = FALSE)
    if ("status" %in% names(previous)) {
      previous_ok <- previous[previous$status == "ok" | is.na(previous$status), , drop = FALSE]
    } else {
      previous_ok <- previous
    }
    done <- unique(cell_key(previous_ok$scenario, previous_ok$rho, previous_ok$replication_id))
    message(sprintf("resuming from %s with %d completed cells", checkpoint_path, length(done)))
  }
  cells <- benchmark_cell_grid(config)
  cells$key <- cell_key(cells$scenario, cells$rho, cells$replication_id)
  remaining <- cells[!(cells$key %in% done), , drop = FALSE]
  rows <- list()
  skipped <- nrow(cells) - nrow(remaining)
  workers <- max(1L, as.integer(workers))
  chunk_size <- max(1L, as.integer(chunk_size))
  if (nrow(remaining) > 0L && workers > 1L) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("parallel package is required for workers > 1")
    }
    workdir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, "workdir", envir = environment())
    parallel::clusterEvalQ(cl, {
      setwd(workdir)
      source(file.path("sim", "src", "support_kernel_benchmark.R"))
      NULL
    })
    for (start in seq(1L, nrow(remaining), by = chunk_size)) {
      end <- min(nrow(remaining), start + chunk_size - 1L)
      chunk <- remaining[start:end, , drop = FALSE]
      message(sprintf(
        "benchmark chunk cells %d-%d of %d using %d workers",
        start, end, nrow(remaining), workers
      ))
      chunk_rows <- parallel::parLapplyLB(
        cl,
        seq_len(nrow(chunk)),
        function(i, chunk, config) {
          z <- chunk[i, , drop = FALSE]
          tryCatch(
            run_one_exact_benchmark_cell(z$scenario[[1]], z$rho[[1]], z$replication_id[[1]], config),
            error = function(e) cell_failure_result(z$scenario[[1]], z$rho[[1]], z$replication_id[[1]], config, e)
          )
        },
        chunk = chunk,
        config = config
      )
      for (cell in chunk_rows) {
        rows[[length(rows) + 1L]] <- cell
        if (checkpoint) {
          append_benchmark_checkpoint(checkpoint_path, cell)
        }
      }
    }
  } else {
    for (i in seq_len(nrow(remaining))) {
      z <- remaining[i, , drop = FALSE]
      message(sprintf(
        "benchmark cell %d of %d: %s rho=%.2f rep=%d",
        i, nrow(remaining), z$scenario[[1]], z$rho[[1]], z$replication_id[[1]]
      ))
      cell <- tryCatch(
        run_one_exact_benchmark_cell(z$scenario[[1]], z$rho[[1]], z$replication_id[[1]], config),
        error = function(e) cell_failure_result(z$scenario[[1]], z$rho[[1]], z$replication_id[[1]], config, e)
      )
      rows[[length(rows) + 1L]] <- cell
      if (checkpoint) {
        append_benchmark_checkpoint(checkpoint_path, cell)
      }
    }
  }
  if (checkpoint && file.exists(checkpoint_path)) {
    detail <- read.csv(checkpoint_path, stringsAsFactors = FALSE)
  } else {
    detail <- do.call(rbind, rows)
  }
  if ("status" %in% names(detail)) {
    keys <- cell_key(detail$scenario, detail$rho, detail$replication_id)
    ok_keys <- unique(keys[detail$status == "ok" | is.na(detail$status)])
    drop_replaced_failures <- detail$status == "failed" &
      detail$method == "Cell failure" &
      keys %in% ok_keys
    if (any(drop_replaced_failures, na.rm = TRUE)) {
      detail <- detail[!drop_replaced_failures, , drop = FALSE]
    }
  }
  summary <- summarize_benchmark_detail(detail)
  write.csv(detail, file.path(table_dir, "table_support_kernel_competitor_benchmark_detail.csv"), row.names = FALSE)
  write.csv(summary, file.path(table_dir, "table_support_kernel_competitor_benchmark_summary.csv"), row.names = FALSE)
  if (checkpoint) {
    message(sprintf("checkpoint path: %s", checkpoint_path))
    message(sprintf("completed cells in this run: %d, skipped: %d", length(rows), skipped))
  }
  list(detail = detail, summary = summary)
}
