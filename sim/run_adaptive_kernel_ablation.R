source(file.path("sim", "design.R"))
source(file.path("sim", "exact_bma.R"))
source(file.path("sim", "src", "family_dictionary.R"))
source(file.path("sim", "src", "adaptive_support_kernel_compression.R"))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) sub("^--mode=", "", args[1]) else "medium"
if (!mode %in% c("smoke", "medium", "full")) stop("mode must be smoke, medium, or full")

out_root <- file.path("sim", "output")
table_dir <- file.path(out_root, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260524)

beta_ablation <- 0.02
tau_ablation <- 1e-3
q0_min_ablation <- 1e-3

rep_count <- switch(mode, smoke = 1L, medium = 3L, full = 12L)
scenario_grid <- switch(mode,
  smoke = "one_representative",
  medium = c("one_representative", "weak_signal"),
  full = c("one_representative", "multi_representative", "weak_signal", "noisy_grouping")
)
rho_grid <- switch(mode, smoke = 0.95, medium = 0.95, full = c(0.8, 0.95))

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
  code_penalty <- beta_ablation * sum(q * costs)
  entropy_penalty <- tau_ablation * sum(ifelse(q > 0, q * log(q), 0))
  out <- data.frame(
    method = method,
    tv = d$tv,
    fkl = d$kl_base_to_compressed,
    rkl = d$kl_compressed_to_base,
    mean_h_error = abs(d$mean_h - 1),
    expected_code = sum(q * costs),
    code_penalty = code_penalty,
    entropy_penalty = entropy_penalty,
    objective_check = d$kl_base_to_compressed + code_penalty + entropy_penalty,
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
  qfit <- optimize_family_mixture(W, alpha, costs, fit$posterior, beta = beta_ablation, tau = tau_ablation, distortion = "fkl", q0_min = q0_min_ablation)
  list(W = W, alpha = alpha, costs = costs, q = qfit$q, n_kernels = ncol(W), objective = qfit$objective, kkt_residual = qfit$kkt_residual, q_iterations = qfit$iterations)
}

fit_topm <- function(fit, top_m = 32L) {
  dict <- make_support_atom_dictionary(fit$supports, fit$posterior, top_m = top_m, include_safety = TRUE)
  W <- family_membership_matrix(fit$supports, dict)
  alpha <- estimate_family_alpha(W, fit$posterior, alpha_floor = 1e-8)$alpha_truncated
  costs <- dict$families$cost
  qfit <- optimize_family_mixture(W, alpha, costs, fit$posterior, beta = beta_ablation, tau = tau_ablation, distortion = "fkl", q0_min = q0_min_ablation)
  list(W = W, alpha = alpha, costs = costs, q = qfit$q, n_kernels = ncol(W), objective = qfit$objective, kkt_residual = qfit$kkt_residual, q_iterations = qfit$iterations)
}

variant_specs <- list(
  full_adaptive = list(label = "Full adaptive", active_sets = TRUE, intervals = TRUE, graph = TRUE, clusters = TRUE, residual = TRUE),
  no_residual_cover = list(label = "No residual", active_sets = TRUE, intervals = TRUE, graph = TRUE, clusters = TRUE, residual = FALSE),
  no_posterior_cluster = list(label = "No cluster", active_sets = TRUE, intervals = TRUE, graph = TRUE, clusters = FALSE, residual = TRUE),
  posterior_cluster_only = list(label = "Cluster only", active_sets = FALSE, intervals = FALSE, graph = FALSE, clusters = TRUE, residual = FALSE),
  residual_cover_only = list(label = "Residual only", active_sets = FALSE, intervals = FALSE, graph = FALSE, clusters = FALSE, residual = TRUE),
  intervals_only = list(label = "Intervals only", active_sets = FALSE, intervals = TRUE, graph = FALSE, clusters = FALSE, residual = FALSE),
  graph_only = list(label = "Graph only", active_sets = FALSE, intervals = FALSE, graph = TRUE, clusters = FALSE, residual = FALSE)
)

fit_pooled_union <- function(fit, dat) {
  base_pool <- make_default_kernel_pool(
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
  residual_pool <- generate_residual_cover_kernels(
    supports = fit$supports,
    residual_weights = fit$posterior,
    group_id = dat$group_id,
    top_medoids = if (mode == "full") 24L else 8L,
    rho_grid = c(0.5, 1, 2),
    group_level = !is.null(dat$group_id)
  )
  kernels <- dedupe_kernels(c(
    list(safety_kernel(cost = max(10, 4 * log(max(ncol(fit$supports), 2))))),
    base_pool,
    residual_pool
  ))
  dict <- new_support_kernel_dictionary(kernels)
  elapsed <- system.time({
    W <- support_kernel_weight_matrix(fit$supports, dict, group_id = dat$group_id)
    alpha <- estimate_kernel_alpha(W, fit$posterior, alpha_floor = 1e-8)$alpha_truncated
    costs <- kernel_costs(dict)
    qfit <- optimize_family_mixture(
      membership = W,
      alpha = alpha,
      costs = costs,
      base_weights = fit$posterior,
      beta = beta_ablation,
      tau = tau_ablation,
      distortion = "fkl",
      max_iter = if (mode == "full") 200L else 140L,
      tol = 1e-7,
      q0_min = q0_min_ablation,
      safety_index = 1L
    )
  })[["elapsed"]]
  list(
    dict = dict,
    W = W,
    alpha = alpha,
    costs = costs,
    q = qfit$q,
    n_kernels = ncol(W),
    objective = qfit$objective,
    kkt_residual = qfit$kkt_residual,
    q_iterations = qfit$iterations,
    elapsed = elapsed
  )
}

fit_pooled_pruned <- function(pooled, fit, dat, mass = 0.99) {
  q <- as.numeric(pooled$q)
  keep <- rep(FALSE, length(q))
  keep[1L] <- TRUE
  for (j in order(q, decreasing = TRUE)) {
    keep[j] <- TRUE
    if (sum(q[keep]) >= mass) break
  }
  dict <- new_support_kernel_dictionary(pooled$dict$kernels[keep])
  q_init <- q[keep]
  q_init <- q_init / sum(q_init)
  elapsed <- system.time({
    W <- support_kernel_weight_matrix(fit$supports, dict, group_id = dat$group_id)
    alpha <- estimate_kernel_alpha(W, fit$posterior, alpha_floor = 1e-8)$alpha_truncated
    costs <- kernel_costs(dict)
    qfit <- optimize_family_mixture(
      membership = W,
      alpha = alpha,
      costs = costs,
      base_weights = fit$posterior,
      beta = beta_ablation,
      tau = tau_ablation,
      distortion = "fkl",
      max_iter = if (mode == "full") 200L else 140L,
      tol = 1e-7,
      q_init = q_init,
      q0_min = q0_min_ablation,
      safety_index = 1L
    )
  })[["elapsed"]]
  list(
    dict = dict,
    W = W,
    alpha = alpha,
    costs = costs,
    q = qfit$q,
    n_kernels = ncol(W),
    retained_q_mass = sum(q[keep]),
    objective = qfit$objective,
    kkt_residual = qfit$kkt_residual,
    q_iterations = qfit$iterations,
    elapsed = elapsed
  )
}

refit_pooled_union_from_pruned <- function(pooled, pruned, fit) {
  pooled_keys <- vapply(pooled$dict$kernels, function(k) k$key, character(1))
  pruned_keys <- vapply(pruned$dict$kernels, function(k) k$key, character(1))
  idx <- match(pruned_keys, pooled_keys)
  q_init <- rep(0, length(pooled_keys))
  q_init[idx] <- pruned$q
  q_init <- q_init / sum(q_init)
  elapsed <- system.time({
    qfit <- optimize_family_mixture(
      membership = pooled$W,
      alpha = pooled$alpha,
      costs = pooled$costs,
      base_weights = fit$posterior,
      beta = beta_ablation,
      tau = tau_ablation,
      distortion = "fkl",
      max_iter = if (mode == "full") 40L else 5L,
      tol = 1e-8,
      q_init = q_init,
      q0_min = q0_min_ablation,
      safety_index = 1L
    )
  })[["elapsed"]]
  pooled$q <- qfit$q
  pooled$objective <- qfit$objective
  pooled$kkt_residual <- qfit$kkt_residual
  pooled$q_iterations <- qfit$iterations
  pooled$elapsed <- pooled$elapsed + elapsed
  pooled$warm_started_from_pruned <- TRUE
  pooled
}

fit_adaptive_variant <- function(variant, spec, fit, dat) {
  pool <- make_default_kernel_pool(
    supports = fit$supports,
    weights = fit$posterior,
    group_id = dat$group_id,
    X = dat$X_train,
    include_intervals = spec$intervals,
    include_graph = spec$graph,
    include_clusters = spec$clusters,
    include_active_sets = spec$active_sets,
    mode = mode
  )
  elapsed <- system.time({
    adapt <- learn_adaptive_kernel_dictionary(
      supports = fit$supports,
      weights = fit$posterior,
      group_id = dat$group_id,
      X = dat$X_train,
      candidate_kernels = pool,
      beta = beta_ablation,
      tau = tau_ablation,
      max_iter = if (mode == "smoke") 5L else if (mode == "medium") 12L else 16L,
      eta = 1e-4,
      residual_cover = spec$residual,
      prune_threshold = 1e-4,
      q0_min = q0_min_ablation,
      mode = mode
    )
  })[["elapsed"]]
  list(
    fit = adapt,
    elapsed = elapsed,
    column_iterations = nrow(adapt$selected),
    trace_rows = nrow(adapt$trace),
    selected_types = if (nrow(adapt$selected)) paste(unique(adapt$selected$type), collapse = ";") else ""
  )
}

se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))

aggregate_summary <- function(detail, group_cols, outfile) {
  metrics <- c(
    "tv", "fkl", "rkl", "q0", "expected_code", "code_penalty", "entropy_penalty",
    "objective_check", "q_effective_kernels",
    "q_kernel95", "rmse_gap", "logscore_gap", "n_kernels",
    "active_kernels_001", "objective", "kkt_residual", "runtime_sec", "column_iterations",
    "retained_q_mass",
    "current_objective", "current_fkl", "current_tv", "current_q0", "current_code",
    "current_active_kernels_001", "best_score", "best_gain", "best_cost"
  )
  metrics <- intersect(metrics, names(detail))
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

rows <- list()
trace_rows <- list()
rid <- 0L
tid <- 0L
for (scenario in scenario_grid) {
  for (rho in rho_grid) {
    for (rep_id in seq_len(rep_count)) {
      seed <- 760000 + rep_id + round(1000 * rho) + match(scenario, scenario_grid) * 10000
      dat <- simulate_custom(scenario, rho, seed)
      fit <- fit_exact_bma(dat$X_train, dat$y_train, theta = 0.08, tau2 = 4, a0 = 1, b0 = 1)
      ref <- predict_mixture(fit, dat$X_train, dat$y_train, dat$X_test, dat$y_test, weights = fit$posterior)

      elapsed <- system.time({ fixed <- fit_fixed_hard(dat, fit, capacity = if (scenario == "multi_representative") 2L else 1L) })[["elapsed"]]
      rid <- rid + 1L
      rows[[rid]] <- evaluate_kernel_fit(
        "fixed_hard_dictionary", fit, dat, fixed$W, fixed$alpha, fixed$q, fixed$costs, ref,
        list(
          label = "Fixed hard",
          n_kernels = fixed$n_kernels,
          active_kernels_001 = NA_real_,
          objective = fixed$objective,
          kkt_residual = fixed$kkt_residual,
          runtime_sec = elapsed,
          column_iterations = 0L,
          q_iterations = fixed$q_iterations,
          selected_types = ""
        )
      )

      elapsed <- system.time({ topm <- fit_topm(fit, top_m = if (mode == "smoke") 12L else 32L) })[["elapsed"]]
      rid <- rid + 1L
      rows[[rid]] <- evaluate_kernel_fit(
        "topM_support_atoms", fit, dat, topm$W, topm$alpha, topm$q, topm$costs, ref,
        list(
          label = "Top-$M$ atoms",
          n_kernels = topm$n_kernels,
          active_kernels_001 = NA_real_,
          objective = topm$objective,
          kkt_residual = topm$kkt_residual,
          runtime_sec = elapsed,
          column_iterations = 0L,
          q_iterations = topm$q_iterations,
          selected_types = ""
        )
      )

      pooled <- fit_pooled_union(fit, dat)
      pooled_pruned <- fit_pooled_pruned(pooled, fit, dat, mass = 0.99)
      pooled <- refit_pooled_union_from_pruned(pooled, pooled_pruned, fit)
      rid <- rid + 1L
      rows[[rid]] <- evaluate_kernel_fit(
        "pooled_candidate_union", fit, dat, pooled$W, pooled$alpha, pooled$q, pooled$costs, ref,
        list(
          label = "Pooled union",
          n_kernels = pooled$n_kernels,
          active_kernels_001 = sum(pooled$q > 0.001),
          objective = pooled$objective,
          kkt_residual = pooled$kkt_residual,
          runtime_sec = pooled$elapsed,
          column_iterations = 0L,
          q_iterations = pooled$q_iterations,
          selected_types = "pooled_union_warmstarted"
        )
      )

      rid <- rid + 1L
      rows[[rid]] <- evaluate_kernel_fit(
        "pooled_pruned_99", fit, dat, pooled_pruned$W, pooled_pruned$alpha, pooled_pruned$q, pooled_pruned$costs, ref,
        list(
          label = "Pooled-pruned 99%",
          n_kernels = pooled_pruned$n_kernels,
          active_kernels_001 = sum(pooled_pruned$q > 0.001),
          objective = pooled_pruned$objective,
          kkt_residual = pooled_pruned$kkt_residual,
          runtime_sec = pooled_pruned$elapsed,
          column_iterations = 0L,
          q_iterations = pooled_pruned$q_iterations,
          retained_q_mass = pooled_pruned$retained_q_mass,
          selected_types = "pooled_pruned_99"
        )
      )

      for (variant in names(variant_specs)) {
        spec <- variant_specs[[variant]]
        va <- fit_adaptive_variant(variant, spec, fit, dat)
        rid <- rid + 1L
        rows[[rid]] <- evaluate_kernel_fit(
          variant, fit, dat, va$fit$W, va$fit$alpha, va$fit$q, kernel_costs(va$fit$dictionary), ref,
          list(
            label = spec$label,
            n_kernels = va$fit$summary$n_kernels,
            active_kernels_001 = va$fit$summary$active_kernels_001,
            objective = va$fit$fit$objective,
            kkt_residual = va$fit$summary$kkt_residual,
            runtime_sec = va$elapsed,
            column_iterations = va$column_iterations,
            q_iterations = va$fit$fit$iterations,
            selected_types = va$selected_types
          )
        )
        if (nrow(va$fit$trace)) {
          tid <- tid + 1L
          trace_rows[[tid]] <- cbind(
            scenario = scenario,
            rho = rho,
            replication_id = rep_id,
            mode = mode,
            method = variant,
            label = spec$label,
            va$fit$trace
          )
        }
      }

      idx <- seq.int(length(rows) - length(variant_specs) - 3L, length(rows))
      for (j in idx) {
        rows[[j]]$scenario <- scenario
        rows[[j]]$rho <- rho
        rows[[j]]$replication_id <- rep_id
        rows[[j]]$mode <- mode
      }
    }
  }
}

all_cols <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(d) {
  missing <- setdiff(all_cols, names(d))
  for (m in missing) d[[m]] <- NA
  d[, all_cols, drop = FALSE]
})
detail <- do.call(rbind, rows)
write.csv(detail, file.path(table_dir, "table_adaptive_kernel_ablation_detail.csv"), row.names = FALSE)
summary <- aggregate_summary(detail, c("scenario", "rho", "method", "label"), file.path(table_dir, "table_adaptive_kernel_ablation.csv"))
write.csv(summary, file.path(table_dir, "table_adaptive_kernel_ablation_summary.csv"), row.names = FALSE)
write.csv(
  subset(summary, method %in% c("pooled_candidate_union", "pooled_pruned_99")),
  file.path(table_dir, "table_adaptive_kernel_ablation_pooled_pruned.csv"),
  row.names = FALSE
)

if (length(trace_rows)) {
  trace <- do.call(rbind, trace_rows)
  write.csv(trace, file.path(table_dir, "table_adaptive_kernel_ablation_convergence.csv"), row.names = FALSE)
  conv <- aggregate_summary(trace, c("scenario", "rho", "method", "label", "iteration"), file.path(table_dir, "table_adaptive_kernel_ablation_convergence_summary.csv"))
}

cat("Adaptive support-kernel ablation complete in", mode, "mode\n")
