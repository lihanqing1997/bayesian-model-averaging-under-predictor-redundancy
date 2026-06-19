source(file.path("sim", "src", "family_mixture_compression.R"))
source(file.path("sim", "src", "support_kernels.R"))

estimate_kernel_alpha <- function(W, base_weights, alpha_floor = 1e-8) {
  estimate_family_alpha(W, base_weights, alpha_floor = alpha_floor)
}

support_kernel_distortions <- function(q, W, alpha, base_weights, h_floor = 1e-12) {
  mixture_distortions(q, W, alpha, base_weights = base_weights, h_floor = h_floor)
}

score_candidate_kernels <- function(W_pool, alpha_pool, h_current, weights, costs,
                                    beta = 0, current_cost = 0,
                                    alpha_floor = 1e-8, h_floor = 1e-10) {
  weights <- normalize_weights(weights)
  alpha_pool <- pmax(alpha_pool, alpha_floor)
  h_current <- pmax(h_current, h_floor)
  G <- sweep(W_pool, 2, alpha_pool, "/")
  gain <- as.numeric(crossprod(weights / h_current, G)) - 1
  score <- gain + beta * (current_cost - costs)
  data.frame(
    candidate_index = seq_along(score),
    gain = gain,
    cost = costs,
    score = score,
    alpha = alpha_pool,
    stringsAsFactors = FALSE
  )
}

best_candidate_from_pool <- function(W_pool, alpha_pool, h_current, weights, costs,
                                     beta = 0, current_cost = 0,
                                     alpha_floor = 1e-8, h_floor = 1e-10) {
  scores <- score_candidate_kernels(W_pool, alpha_pool, h_current, weights, costs, beta, current_cost, alpha_floor, h_floor)
  scores[which.max(scores$score), , drop = FALSE]
}

make_default_kernel_pool <- function(supports, weights, group_id = NULL, X = NULL,
                                     include_intervals = TRUE,
                                     include_graph = TRUE,
                                     include_clusters = TRUE,
                                     include_active_sets = TRUE,
                                     include_metric_balls = TRUE,
                                     mode = "medium") {
  p <- ncol(supports)
  kernels <- list()
  if (include_active_sets && !is.null(group_id)) {
    kernels <- c(
      kernels,
      generate_active_set_kernel_pool(
        supports = supports,
        group_id = group_id,
        weights = weights,
        max_sets = if (mode == "full") 96L else 32L,
        capacities = c(1L, 2L),
        rho_grid = c(1, 2),
        include_hard = TRUE
      )
    )
  }
  if (include_clusters) {
    kernels <- c(
      kernels,
      generate_posterior_cluster_kernels(
        supports = supports,
        weights = weights,
        group_id = group_id,
        top_medoids = if (mode == "full") 48L else 16L,
        rho_grid = c(0.5, 1, 2),
        group_level = !is.null(group_id)
      )
    )
  }
  if (include_metric_balls) {
    kernels <- c(
      kernels,
      generate_metric_ball_kernel_pool(
        supports = supports,
        weights = weights,
        group_id = group_id,
        top_centers = if (mode == "full") 48L else 16L,
        radii = 0:3,
        rho_grid = c(0.5, 1, 2, 4),
        include_hard = TRUE,
        include_soft = TRUE,
        group_level = !is.null(group_id)
      )
    )
  }
  if (include_intervals) {
    kernels <- c(
      kernels,
      generate_interval_kernel_pool(
        p = p,
        lengths = unique(pmin(c(2, 3, 4, 6, 8, 12, 16), p)),
        rho_grid = c(0.5, 1.5),
        capacities = c(1L, 2L, Inf)
      )
    )
  }
  if (include_graph && !is.null(X)) {
    kernels <- c(kernels, generate_graph_kernel_pool(X, thresholds = c(0.5, 0.7, 0.9), rho_grid = c(0.5, 1.5)))
  }
  dedupe_kernels(kernels)
}

learn_adaptive_kernel_dictionary <- function(supports,
                                             weights,
                                             group_id = NULL,
                                             X = NULL,
                                             candidate_kernels = NULL,
                                             seed_kernels = NULL,
                                             beta = 0.02,
                                             tau = 1e-3,
                                             alpha_floor = 1e-8,
                                             max_iter = 12L,
                                             eta = 1e-4,
                                             q0_min = 1e-3,
                                             q0_max = 1,
                                             prune_threshold = 1e-4,
                                             residual_cover = TRUE,
                                             mode = "medium",
                                             verbose = FALSE) {
  weights <- normalize_weights(weights)
  p <- ncol(supports)
  if (is.null(seed_kernels)) {
    seed_kernels <- list(safety_kernel(cost = max(10, 4 * log(max(p, 2)))))
  } else if (!any(vapply(seed_kernels, function(k) identical(k$type, "safety"), logical(1)))) {
    seed_kernels <- c(list(safety_kernel(cost = max(10, 4 * log(max(p, 2))))), seed_kernels)
  }
  if (is.null(candidate_kernels)) {
    candidate_kernels <- make_default_kernel_pool(supports, weights, group_id, X, mode = mode)
  }
  current <- dedupe_kernels(seed_kernels)
  pool <- dedupe_kernels(candidate_kernels)
  trace <- list()
  selected <- data.frame()

  fit_current <- function(kernels, q_init = NULL) {
    dict <- new_support_kernel_dictionary(kernels)
    W <- support_kernel_weight_matrix(supports, dict, group_id = group_id)
    alpha <- estimate_kernel_alpha(W, weights, alpha_floor = alpha_floor)$alpha_truncated
    costs <- kernel_costs(dict)
    if (!is.null(q_init) && length(q_init) != length(kernels)) {
      q_init <- NULL
    }
    fit <- optimize_family_mixture(
      membership = W,
      alpha = alpha,
      costs = costs,
      base_weights = weights,
      beta = beta,
      tau = tau,
      distortion = "fkl",
      max_iter = if (mode == "full") 120L else 80L,
      tol = 1e-6,
      q_init = q_init,
      q0_min = q0_min,
      q0_max = q0_max,
      safety_index = 1L
    )
    list(dict = dict, W = W, alpha = alpha, costs = costs, fit = fit)
  }

  prune_current <- function(kernels, state) {
    if (prune_threshold <= 0 || length(kernels) <= 1L) {
      return(list(kernels = kernels, state = state, pruned = 0L))
    }
    keep <- vapply(kernels, function(k) identical(k$type, "safety"), logical(1)) |
      state$fit$q > prune_threshold
    if (all(keep) || !any(keep)) {
      return(list(kernels = kernels, state = state, pruned = 0L))
    }
    kernels2 <- kernels[keep]
    q_init <- state$fit$q[keep]
    q_init <- q_init / sum(q_init)
    list(kernels = kernels2, state = fit_current(kernels2, q_init = q_init), pruned = sum(!keep))
  }

  state <- fit_current(current)
  p0 <- prune_current(current, state)
  current <- p0$kernels
  state <- p0$state
  for (iter in seq_len(max_iter)) {
    existing_keys <- vapply(current, function(k) k$key, character(1))
    pool_iter <- pool
    if (residual_cover) {
      h_for_residual <- pmax(mixture_h(state$fit$q, state$W, state$alpha), max(q0_min, 1e-10))
      residual_weights <- weights / h_for_residual
      pool_iter <- dedupe_kernels(c(
        pool_iter,
        generate_residual_cover_kernels(
          supports = supports,
          residual_weights = residual_weights,
          group_id = group_id,
          top_medoids = if (mode == "full") 24L else 8L,
          rho_grid = c(0.5, 1, 2),
          group_level = !is.null(group_id)
        )
      ))
    }
    pool_remaining <- pool_iter[!vapply(pool_iter, function(k) k$key %in% existing_keys, logical(1))]
    if (!length(pool_remaining)) break
    pool_dict <- new_support_kernel_dictionary(pool_remaining)
    W_pool <- support_kernel_weight_matrix(supports, pool_dict, group_id = group_id)
    alpha_pool <- estimate_kernel_alpha(W_pool, weights, alpha_floor = alpha_floor)$alpha_truncated
    h_current <- mixture_h(state$fit$q, state$W, state$alpha)
    scores <- score_candidate_kernels(
      W_pool = W_pool,
      alpha_pool = alpha_pool,
      h_current = h_current,
      weights = weights,
      costs = kernel_costs(pool_dict),
      beta = beta,
      current_cost = sum(state$fit$q * state$costs),
      alpha_floor = alpha_floor
    )
    best <- scores[which.max(scores$score), , drop = FALSE]
    best_kernel <- pool_remaining[[best$candidate_index]]
    trace[[length(trace) + 1L]] <- data.frame(
      iteration = iter,
      current_kernels = length(current),
      best_kernel = best_kernel$name,
      best_type = best_kernel$type,
      best_gain = best$gain,
      best_cost = best$cost,
      best_score = best$score,
      current_objective = state$fit$objective,
      current_fkl = state$fit$distortions$kl_base_to_compressed,
      current_tv = state$fit$distortions$tv,
      current_q0 = state$fit$q[1],
      current_code = state$fit$expected_cost,
      current_active_kernels_001 = sum(state$fit$q > 0.001),
      kkt_residual = state$fit$kkt_residual,
      stringsAsFactors = FALSE
    )
    if (verbose) {
      message(sprintf("iter %d score %.4f kernel %s", iter, best$score, best_kernel$name))
    }
    if (!is.finite(best$score) || best$score <= eta) break
    old_q <- state$fit$q
    q_seed <- min(0.05, max(0.01, 1 / (10 * (length(old_q) + 1))))
    q_init <- c((1 - q_seed) * old_q, q_seed)
    current <- c(current, list(best_kernel))
    selected <- rbind(
      selected,
      data.frame(
        iteration = iter,
        kernel = best_kernel$name,
        type = best_kernel$type,
        score = best$score,
        gain = best$gain,
        cost = best$cost,
        pruned_after_refit = NA_integer_,
        stringsAsFactors = FALSE
      )
    )
    state <- fit_current(current, q_init = q_init)
    pruned <- prune_current(current, state)
    current <- pruned$kernels
    state <- pruned$state
    if (nrow(selected)) {
      selected$pruned_after_refit[nrow(selected)] <- pruned$pruned
    }
  }

  d <- state$fit$distortions
  q <- state$fit$q
  entropy <- -sum(ifelse(q > 0, q * log(q), 0))
  summary <- data.frame(
    method = "adaptive_support_kernel",
    tv = d$tv,
    fkl = d$kl_base_to_compressed,
    rkl = d$kl_compressed_to_base,
    mean_h_error = abs(d$mean_h - 1),
    expected_code = sum(q * state$costs),
    q0 = q[1],
    q_entropy = entropy,
    q_effective_kernels = exp(entropy),
    q_kernel95 = which(cumsum(sort(q, decreasing = TRUE)) >= 0.95)[1],
    active_kernels_001 = sum(q > 0.001),
    n_kernels = length(state$dict$kernels),
    beta = beta,
    tau = tau,
    q0_min = q0_min,
    q0_max = q0_max,
    prune_threshold = prune_threshold,
    kkt_residual = state$fit$kkt_residual,
    optimizer_status = state$fit$status,
    objective_change = state$fit$objective_change,
    converged = state$fit$converged,
    stringsAsFactors = FALSE
  )
  list(
    dictionary = state$dict,
    W = state$W,
    alpha = state$alpha,
    q = q,
    fit = state$fit,
    summary = summary,
    trace = if (length(trace)) do.call(rbind, trace) else data.frame(),
    selected = selected
  )
}

topm_atom_tv_uniform <- function(group_size, active_groups, M_grid) {
  M_star <- group_size^active_groups
  data.frame(
    group_size = group_size,
    active_groups = active_groups,
    M_star = M_star,
    atoms = M_grid,
    topM_tv = pmax(0, 1 - pmin(M_grid, M_star) / M_star),
    family_tv = 0,
    family_regions = 1L,
    stringsAsFactors = FALSE
  )
}

prune_kernel_indices_by_q_mass <- function(q, safety_index = 1L, mass = 0.99) {
  q <- normalize_weights(q)
  mass <- min(max(as.numeric(mass), 0), 1)
  keep <- rep(FALSE, length(q))
  keep[safety_index] <- TRUE
  nonsafety <- setdiff(seq_along(q), safety_index)
  target <- mass * sum(q[nonsafety])
  carried <- 0
  for (j in nonsafety[order(q[nonsafety], decreasing = TRUE)]) {
    keep[j] <- TRUE
    carried <- carried + q[j]
    if (carried >= target) break
  }
  which(keep)
}

summarize_kernel_fit <- function(label, fit, costs, safety_index = 1L) {
  q <- fit$q
  entropy <- -sum(ifelse(q > 0, q * log(q), 0))
  data.frame(
    label = label,
    tv = fit$distortions$tv,
    fkl = fit$distortions$kl_base_to_compressed,
    rkl = fit$distortions$kl_compressed_to_base,
    mean_h_error = abs(fit$distortions$mean_h - 1),
    expected_code = sum(q * costs),
    q0 = q[safety_index],
    q_effective_kernels = exp(entropy),
    active_kernels_001 = sum(q > 0.001),
    objective = fit$objective,
    kkt_residual = fit$kkt_residual,
    iterations = fit$iterations,
    converged = fit$converged,
    optimizer_status = fit$status,
    objective_change = fit$objective_change,
    stringsAsFactors = FALSE
  )
}

kernel_pool_count_table <- function(kernels, stage) {
  if (!length(kernels)) {
    return(data.frame(stage = stage, origin = character(0), type = character(0), n = integer(0)))
  }
  x <- data.frame(
    stage = stage,
    origin = vapply(kernels, kernel_origin, character(1)),
    type = vapply(kernels, function(k) k$type, character(1)),
    stringsAsFactors = FALSE
  )
  out <- aggregate(rep(1L, nrow(x)), x, sum)
  names(out)[ncol(out)] <- "n"
  out[order(out$stage, out$origin, out$type), , drop = FALSE]
}

map_q_to_kernel_pool <- function(source_kernels, source_q, target_kernels) {
  source_keys <- vapply(source_kernels, function(k) k$key, character(1))
  target_keys <- vapply(target_kernels, function(k) k$key, character(1))
  q <- rep(0, length(target_kernels))
  idx <- match(source_keys, target_keys)
  ok <- !is.na(idx)
  q[idx[ok]] <- source_q[ok]
  if (sum(q) <= 0) {
    q <- rep(1 / length(q), length(q))
  } else {
    q <- q / sum(q)
  }
  q
}

fit_askpc_pooled_pruned <- function(supports,
                                    weights,
                                    group_id = NULL,
                                    X = NULL,
                                    beta = 0.02,
                                    tau = 1e-3,
                                    alpha_floor = 1e-8,
                                    q0_min = 1e-3,
                                    q0_max = 1,
                                    prune_mass = 0.99,
                                    mode = "medium",
                                    include_intervals = TRUE,
                                    include_graph = TRUE,
                                    include_clusters = TRUE,
                                    include_active_sets = TRUE,
                                    include_metric_balls = TRUE,
                                    include_residual = TRUE,
                                    safety_cost = NULL,
                                    max_iter = NULL,
                                    tol = 1e-7,
                                    kkt_tol = 1e-4,
                                    near_duplicate_tol = 1e-10,
                                    refit_rounds = 2L,
                                    objective_tol = 1e-7,
                                    polish_pooled = FALSE,
                                    fit_mode = c("exact_weighted", "sample_split")) {
  fit_mode <- match.arg(fit_mode)
  supports <- as.matrix(supports)
  weights <- normalize_weights(weights)
  p <- ncol(supports)
  if (is.null(safety_cost)) {
    safety_cost <- max(10, 4 * log(max(p, 2)))
  }
  base_pool <- make_default_kernel_pool(
    supports = supports,
    weights = weights,
    group_id = group_id,
    X = X,
    include_intervals = include_intervals,
    include_graph = include_graph,
    include_clusters = include_clusters,
    include_active_sets = include_active_sets,
    include_metric_balls = include_metric_balls,
    mode = mode
  )
  residual_pool <- list()
  if (include_residual) {
    residual_pool <- generate_residual_cover_kernels(
      supports = supports,
      residual_weights = weights,
      group_id = group_id,
      top_medoids = if (mode == "full") 24L else 8L,
      rho_grid = c(0.5, 1, 2),
      group_level = !is.null(group_id)
    )
  }
  raw_kernels <- c(list(safety_kernel(cost = safety_cost)), base_pool, residual_pool)
  deduped_kernels <- dedupe_kernels(raw_kernels)
  near_dict <- new_support_kernel_dictionary(deduped_kernels)
  W_near <- support_kernel_weight_matrix(supports, near_dict, group_id = group_id)
  keep_near <- near_duplicate_kernel_keep(
    W = W_near,
    base_weights = weights,
    alpha_floor = alpha_floor,
    tol = near_duplicate_tol,
    always_keep = 1L
  )
  kernels <- deduped_kernels[keep_near]
  pool_counts <- rbind(
    kernel_pool_count_table(raw_kernels, "raw"),
    kernel_pool_count_table(deduped_kernels, "deduped"),
    kernel_pool_count_table(kernels, "near_deduped")
  )
  if (is.null(max_iter)) {
    max_iter <- if (mode == "full") 240L else 160L
  }

  fit_dictionary <- function(kernels, q_init = NULL, iter = max_iter, polish = polish_pooled) {
    dict <- new_support_kernel_dictionary(kernels)
    W <- support_kernel_weight_matrix(supports, dict, group_id = group_id)
    alpha <- estimate_kernel_alpha(W, weights, alpha_floor = alpha_floor)$alpha_truncated
    costs <- kernel_costs(dict)
    fit <- optimize_family_mixture(
      membership = W,
      alpha = alpha,
      costs = costs,
      base_weights = weights,
      beta = beta,
      tau = tau,
      distortion = "fkl",
      max_iter = iter,
      tol = tol,
      kkt_tol = kkt_tol,
      polish = polish,
      polish_maxit = if (mode == "full") 120L else 80L,
      q_init = q_init,
      q0_min = q0_min,
      q0_max = q0_max,
      safety_index = 1L
    )
    list(dict = dict, W = W, alpha = alpha, costs = costs, fit = fit)
  }

  evaluate_dictionary <- function(kernels, q, status = "feasible_evaluation") {
    dict <- new_support_kernel_dictionary(kernels)
    W <- support_kernel_weight_matrix(supports, dict, group_id = group_id)
    alpha <- estimate_kernel_alpha(W, weights, alpha_floor = alpha_floor)$alpha_truncated
    costs <- kernel_costs(dict)
    fit <- evaluate_family_mixture_fit(
      membership = W,
      alpha = alpha,
      costs = costs,
      q = q,
      base_weights = weights,
      beta = beta,
      tau = tau,
      distortion = "fkl",
      q0_min = q0_min,
      q0_max = q0_max,
      safety_index = 1L,
      status = status
    )
    list(dict = dict, W = W, alpha = alpha, costs = costs, fit = fit)
  }

  pooled <- fit_dictionary(kernels)
  keep_idx <- integer(0)
  pruned <- NULL
  refit_rounds <- max(1L, as.integer(refit_rounds))
  for (round_id in seq_len(refit_rounds)) {
    keep_idx <- prune_kernel_indices_by_q_mass(pooled$fit$q, safety_index = 1L, mass = prune_mass)
    pruned_kernels <- pooled$dict$kernels[keep_idx]
    q_init_pruned <- pooled$fit$q[keep_idx]
    q_init_pruned <- q_init_pruned / sum(q_init_pruned)
    pruned <- fit_dictionary(pruned_kernels, q_init = q_init_pruned)
    q_init_pooled <- map_q_to_kernel_pool(pruned$dict$kernels, pruned$fit$q, pooled$dict$kernels)
    pooled_new <- fit_dictionary(kernels, q_init = q_init_pooled, iter = max_iter)
    if (pooled_new$fit$objective <= pooled$fit$objective + objective_tol) {
      pooled <- pooled_new
    }
  }
  keep_idx <- prune_kernel_indices_by_q_mass(pooled$fit$q, safety_index = 1L, mass = prune_mass)
  pruned_kernels <- pooled$dict$kernels[keep_idx]
  q_init_pruned <- pooled$fit$q[keep_idx]
  q_init_pruned <- q_init_pruned / sum(q_init_pruned)
  pruned <- fit_dictionary(pruned_kernels, q_init = q_init_pruned)
  q_init_pooled <- map_q_to_kernel_pool(pruned$dict$kernels, pruned$fit$q, pooled$dict$kernels)
  pooled_warm <- fit_dictionary(kernels, q_init = q_init_pooled, iter = max_iter * 2L)
  if (pooled_warm$fit$objective <= pooled$fit$objective + objective_tol) {
    pooled <- pooled_warm
    keep_idx <- prune_kernel_indices_by_q_mass(pooled$fit$q, safety_index = 1L, mass = prune_mass)
    pruned_kernels <- pooled$dict$kernels[keep_idx]
    q_init_pruned <- pooled$fit$q[keep_idx]
    q_init_pruned <- q_init_pruned / sum(q_init_pruned)
    pruned <- fit_dictionary(pruned_kernels, q_init = q_init_pruned)
  }
  q_embed <- map_q_to_kernel_pool(pruned$dict$kernels, pruned$fit$q, pooled$dict$kernels)
  pooled_solution_source <- "pooled_optimizer"
  pooled_embedded <- evaluate_dictionary(
    kernels = pooled$dict$kernels,
    q = q_embed,
    status = "embedded_pruned_feasible_solution"
  )
  if (pooled_embedded$fit$objective <= pooled$fit$objective + objective_tol) {
    pooled <- pooled_embedded
    pooled_solution_source <- "embedded_pruned_fallback"
  }
  retained_q_mass <- sum(pooled$fit$q[keep_idx])
  objective_gap <- pooled$fit$objective - pruned$fit$objective
  pooled_tol <- if (is.finite(pooled$fit$tol)) pooled$fit$tol else objective_tol
  objective_consistent <- objective_gap <= max(objective_tol, 10 * pooled_tol * (1 + abs(pooled$fit$objective)))

  pool_diagnostics <- data.frame(
    raw_candidates = length(raw_kernels),
    key_deduped_candidates = length(deduped_kernels),
    near_deduped_candidates = length(kernels),
    key_duplicates_removed = length(raw_kernels) - length(deduped_kernels),
    near_duplicates_removed = length(deduped_kernels) - length(kernels),
    prune_mass = prune_mass,
    retained_q_mass = retained_q_mass,
    pooled_objective = pooled$fit$objective,
    pruned_objective = pruned$fit$objective,
    objective_gap_pooled_minus_pruned = objective_gap,
    objective_consistent = objective_consistent,
    pooled_solution_source = pooled_solution_source,
    fit_mode = fit_mode,
    stringsAsFactors = FALSE
  )

  list(
    pooled = pooled,
    pruned = pruned,
    prune_mass = prune_mass,
    retained_q_mass = retained_q_mass,
    pool_counts = pool_counts,
    pool_diagnostics = pool_diagnostics,
    pooled_top_kernels = kernel_summary_table(pooled$dict, pooled$fit$q, pooled$alpha, top_n = 25L),
    pruned_top_kernels = kernel_summary_table(pruned$dict, pruned$fit$q, pruned$alpha, top_n = 25L),
    summary = rbind(
      cbind(
        method = "pooled_candidate_union",
        summarize_kernel_fit("Pooled union", pooled$fit, pooled$costs),
        pool_diagnostics[, c("raw_candidates", "key_deduped_candidates", "near_deduped_candidates", "objective_consistent", "pooled_solution_source", "fit_mode"), drop = FALSE]
      ),
      cbind(
        method = "pooled_pruned_99",
        summarize_kernel_fit("Pooled-pruned 99%", pruned$fit, pruned$costs),
        pool_diagnostics[, c("raw_candidates", "key_deduped_candidates", "near_deduped_candidates", "objective_consistent", "pooled_solution_source", "fit_mode"), drop = FALSE]
      )
    )
  )
}
