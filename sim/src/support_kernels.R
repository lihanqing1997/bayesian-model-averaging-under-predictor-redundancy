if (!exists("support_group_counts")) {
  support_group_counts <- function(supports, group_id) {
    K <- max(group_id)
    counts <- matrix(0, nrow = nrow(supports), ncol = K)
    for (k in seq_len(K)) {
      counts[, k] <- rowSums(supports[, group_id == k, drop = FALSE])
    }
    counts
  }
}

kernel_key <- function(type, payload) {
  paste(type, paste(names(payload), unlist(payload), sep = "=", collapse = ","), sep = ":")
}

new_support_kernel_dictionary <- function(kernels) {
  if (!length(kernels)) {
    stop("kernel dictionary cannot be empty")
  }
  keys <- vapply(kernels, function(k) k$key, character(1))
  kernels <- kernels[!duplicated(keys)]
  structure(list(kernels = kernels), class = "support_kernel_dictionary")
}

safety_kernel <- function(cost = 100) {
  list(
    type = "safety",
    name = "safety_unrestricted",
    cost = cost,
    key = "safety",
    params = list()
  )
}

hard_active_set_kernel <- function(active_groups, capacity = 1L, group_id, cost = NULL, name = NULL) {
  active_groups <- sort(unique(as.integer(active_groups)))
  K <- max(group_id)
  cap <- rep(0L, K)
  cap[active_groups] <- as.integer(capacity)
  if (is.null(cost)) {
    cost <- length(active_groups) * log(max(K, 2L)) + length(active_groups) * log(max(capacity + 1L, 2L))
  }
  if (is.null(name)) {
    name <- paste0("hard_active_", if (length(active_groups)) paste(active_groups, collapse = "-") else "empty")
  }
  list(
    type = "hard_active",
    name = name,
    cost = cost,
    key = kernel_key("hard_active", list(S = paste(active_groups, collapse = "-"), b = capacity)),
    params = list(active_groups = active_groups, capacity = cap)
  )
}

soft_capacity_kernel <- function(active_groups, capacity = 1L, group_id, rho = 2, cost = NULL, name = NULL) {
  active_groups <- sort(unique(as.integer(active_groups)))
  K <- max(group_id)
  cap <- rep(0L, K)
  cap[active_groups] <- as.integer(capacity)
  if (is.null(cost)) {
    cost <- length(active_groups) * log(max(K, 2L)) + log1p(rho) + length(active_groups) * log(max(capacity + 1L, 2L))
  }
  if (is.null(name)) {
    name <- paste0("soft_capacity_", if (length(active_groups)) paste(active_groups, collapse = "-") else "empty")
  }
  list(
    type = "soft_capacity",
    name = name,
    cost = cost,
    key = kernel_key("soft_capacity", list(S = paste(active_groups, collapse = "-"), b = capacity, rho = rho)),
    params = list(active_groups = active_groups, capacity = cap, rho = rho)
  )
}

soft_hamming_kernel <- function(center_support, rho = 1, cost = NULL, name = NULL) {
  center_support <- as.integer(center_support > 0)
  if (is.null(cost)) {
    cost <- log(length(center_support) + 1) + sum(center_support) * log(2)
  }
  if (is.null(name)) {
    active <- which(center_support > 0)
    name <- paste0("soft_hamming_", paste(head(active, 12), collapse = "-"))
  }
  list(
    type = "soft_hamming",
    name = name,
    cost = cost + log1p(rho),
    key = kernel_key("soft_hamming", list(center = paste(which(center_support > 0), collapse = "-"), rho = rho)),
    params = list(center_support = center_support, rho = rho)
  )
}

soft_group_hamming_kernel <- function(center_groups, group_id, rho = 1, cost = NULL, name = NULL) {
  center_groups <- sort(unique(as.integer(center_groups)))
  K <- max(group_id)
  center_vec <- integer(K)
  center_vec[center_groups] <- 1L
  if (is.null(cost)) {
    cost <- length(center_groups) * log(max(K, 2L)) + log1p(rho)
  }
  if (is.null(name)) {
    name <- paste0("soft_group_", if (length(center_groups)) paste(center_groups, collapse = "-") else "empty")
  }
  list(
    type = "soft_group_hamming",
    name = name,
    cost = cost,
    key = kernel_key("soft_group_hamming", list(center = paste(center_groups, collapse = "-"), rho = rho)),
    params = list(center_groups = center_groups, center_vec = center_vec, rho = rho)
  )
}

soft_interval_kernel <- function(interval, p, capacity = Inf, rho = 1, cost = NULL, name = NULL) {
  interval <- sort(unique(as.integer(interval)))
  interval <- interval[interval >= 1L & interval <= p]
  if (is.null(cost)) {
    cost <- log(max(p, 2L)) + log(max(length(interval), 2L)) + log1p(rho) + log1p(ifelse(is.finite(capacity), capacity, p))
  }
  if (is.null(name)) {
    name <- paste0("soft_interval_", min(interval), "_", max(interval))
  }
  list(
    type = "soft_interval",
    name = name,
    cost = cost,
    key = kernel_key("soft_interval", list(a = min(interval), z = max(interval), u = capacity, rho = rho)),
    params = list(interval = interval, p = p, capacity = capacity, rho = rho)
  )
}

soft_graph_neighborhood_kernel <- function(nodes, p, rho = 1, cost = NULL, name = NULL) {
  nodes <- sort(unique(as.integer(nodes)))
  nodes <- nodes[nodes >= 1L & nodes <= p]
  if (is.null(cost)) {
    cost <- log(max(p, 2L)) + length(nodes) * log(2) + log1p(rho)
  }
  if (is.null(name)) {
    name <- paste0("soft_graph_", paste(head(nodes, 8), collapse = "-"))
  }
  list(
    type = "soft_graph",
    name = name,
    cost = cost,
    key = kernel_key("soft_graph", list(nodes = paste(nodes, collapse = "-"), rho = rho)),
    params = list(nodes = nodes, p = p, rho = rho)
  )
}

kernel_weights_one <- function(supports, kernel, group_id = NULL) {
  supports <- as.matrix(supports)
  n <- nrow(supports)
  type <- kernel$type
  if (identical(type, "safety")) {
    return(rep(1, n))
  }
  if (type %in% c("hard_active", "soft_capacity", "soft_group_hamming")) {
    if (is.null(group_id)) {
      stop("group_id required for group-based kernels")
    }
    counts <- support_group_counts(supports, group_id)
    active <- counts > 0
  }
  if (identical(type, "hard_active")) {
    cap <- kernel$params$capacity
    ok <- matrix(TRUE, nrow = n, ncol = 1)
    for (k in seq_len(ncol(counts))) {
      if (cap[k] == 0L) {
        ok <- ok & counts[, k, drop = FALSE] == 0
      } else {
        ok <- ok & counts[, k, drop = FALSE] <= cap[k]
      }
    }
    return(as.numeric(ok[, 1]))
  }
  if (identical(type, "soft_capacity")) {
    cap <- kernel$params$capacity
    rho <- kernel$params$rho
    inactive_penalty <- rowSums(counts[, cap == 0L, drop = FALSE])
    over_penalty <- rowSums(pmax(sweep(counts[, cap > 0L, drop = FALSE], 2, cap[cap > 0L], "-"), 0))
    return(exp(-rho * (inactive_penalty + over_penalty)))
  }
  if (identical(type, "soft_hamming")) {
    center <- kernel$params$center_support
    rho <- kernel$params$rho
    d <- rowSums(abs(sweep(supports, 2, center, "-")))
    return(exp(-rho * d))
  }
  if (identical(type, "soft_group_hamming")) {
    center <- kernel$params$center_vec
    rho <- kernel$params$rho
    d <- rowSums(abs(sweep(active * 1, 2, center, "-")))
    return(exp(-rho * d))
  }
  if (identical(type, "soft_interval")) {
    interval <- kernel$params$interval
    rho <- kernel$params$rho
    capacity <- kernel$params$capacity
    outside <- rowSums(supports[, -interval, drop = FALSE])
    inside_count <- rowSums(supports[, interval, drop = FALSE])
    over <- if (is.finite(capacity)) pmax(inside_count - capacity, 0) else 0
    return(exp(-rho * (outside + over)))
  }
  if (identical(type, "soft_graph")) {
    nodes <- kernel$params$nodes
    rho <- kernel$params$rho
    outside <- rowSums(supports[, -nodes, drop = FALSE])
    return(exp(-rho * outside))
  }
  stop("unknown support kernel type: ", type)
}

support_kernel_weight_matrix <- function(supports, dictionary, group_id = NULL) {
  W <- sapply(dictionary$kernels, function(k) kernel_weights_one(supports, k, group_id = group_id))
  W <- as.matrix(W)
  storage.mode(W) <- "double"
  if (any(W < -1e-12 | W > 1 + 1e-12 | !is.finite(W))) {
    stop("kernel weights must be finite and in [0,1]")
  }
  pmin(pmax(W, 0), 1)
}

kernel_costs <- function(dictionary) {
  vapply(dictionary$kernels, kernel_code_length, numeric(1))
}

kernel_names <- function(dictionary) {
  vapply(dictionary$kernels, function(k) k$name, character(1))
}

generate_active_set_kernel_pool <- function(supports, group_id, weights, max_sets = 32L,
                                            capacities = c(1L, 2L), rho_grid = c(1, 2),
                                            include_hard = TRUE) {
  active_sets <- candidate_active_sets_from_supports(
    supports = supports,
    group_id = group_id,
    weights = weights,
    max_sets = max_sets,
    include_empty = TRUE
  )
  kernels <- list()
  for (S in active_sets) {
    for (b in capacities) {
      if (include_hard) {
        kernels[[length(kernels) + 1L]] <- hard_active_set_kernel(S, b, group_id)
      }
      for (rho in rho_grid) {
        kernels[[length(kernels) + 1L]] <- soft_capacity_kernel(S, b, group_id, rho = rho)
      }
    }
  }
  kernels
}

generate_posterior_cluster_kernels <- function(supports, weights, group_id = NULL,
                                               top_medoids = 16L, rho_grid = c(0.5, 1, 2),
                                               group_level = TRUE) {
  ord <- order(weights, decreasing = TRUE)
  centers <- supports[head(ord, min(top_medoids, length(ord))), , drop = FALSE]
  kernels <- list()
  for (i in seq_len(nrow(centers))) {
    for (rho in rho_grid) {
      kernels[[length(kernels) + 1L]] <- soft_hamming_kernel(centers[i, ], rho = rho)
      if (group_level && !is.null(group_id)) {
        cnt <- support_group_counts(centers[i, , drop = FALSE], group_id)
        kernels[[length(kernels) + 1L]] <- soft_group_hamming_kernel(which(cnt[1, ] > 0), group_id, rho = rho)
      }
    }
  }
  kernels
}

generate_residual_cover_kernels <- function(supports, residual_weights, group_id = NULL,
                                            top_medoids = 12L, rho_grid = c(0.5, 1, 2),
                                            group_level = TRUE) {
  residual_weights <- as.numeric(residual_weights)
  residual_weights[!is.finite(residual_weights) | residual_weights < 0] <- 0
  if (sum(residual_weights) <= 0) {
    residual_weights <- rep(1, nrow(supports))
  }
  kernels <- generate_posterior_cluster_kernels(
    supports = supports,
    weights = residual_weights,
    group_id = group_id,
    top_medoids = top_medoids,
    rho_grid = rho_grid,
    group_level = group_level
  )
  lapply(kernels, function(k) {
    k$name <- paste0("residual_", k$name)
    k$key <- paste0("residual:", k$key)
    k$origin <- "residual_cover"
    k$cost <- k$cost + log(2)
    k
  })
}

generate_interval_kernel_pool <- function(p, lengths = c(2, 3, 4, 6, 8, 12),
                                          rho_grid = c(0.5, 1, 2),
                                          capacities = c(1L, 2L, Inf)) {
  kernels <- list()
  for (len in lengths[lengths <= p]) {
    starts <- unique(c(seq(1, p, by = max(1L, floor(len / 2))), p - len + 1L))
    starts <- starts[starts >= 1L & starts + len - 1L <= p]
    for (s in starts) {
      interval <- s:(s + len - 1L)
      for (cap in capacities) {
        for (rho in rho_grid) {
          kernels[[length(kernels) + 1L]] <- soft_interval_kernel(interval, p = p, capacity = cap, rho = rho)
        }
      }
    }
  }
  kernels
}

generate_graph_kernel_pool <- function(X, thresholds = c(0.5, 0.7, 0.9), rho_grid = c(0.5, 1, 2)) {
  p <- ncol(X)
  C <- abs(stats::cor(X))
  C[!is.finite(C)] <- 0
  diag(C) <- 0
  kernels <- list()
  for (thr in thresholds) {
    seen <- rep(FALSE, p)
    adj <- C >= thr
    for (j in seq_len(p)) {
      if (seen[j]) next
      queue <- j
      comp <- integer(0)
      seen[j] <- TRUE
      while (length(queue)) {
        v <- queue[1]
        queue <- queue[-1]
        comp <- c(comp, v)
        nb <- which(adj[v, ] & !seen)
        if (length(nb)) {
          seen[nb] <- TRUE
          queue <- c(queue, nb)
        }
      }
      if (length(comp) >= 2L) {
        for (rho in rho_grid) {
          kernels[[length(kernels) + 1L]] <- soft_graph_neighborhood_kernel(comp, p, rho = rho)
        }
      }
    }
  }
  kernels
}

dedupe_kernels <- function(kernels) {
  if (!length(kernels)) return(kernels)
  keys <- vapply(kernels, function(k) k$key, character(1))
  kernels[!duplicated(keys)]
}

kernel_origin <- function(kernel) {
  if (!is.null(kernel$origin)) {
    return(as.character(kernel$origin))
  }
  switch(
    kernel$type,
    safety = "safety",
    hard_active = "response_independent_hard",
    soft_capacity = "response_independent_capacity",
    soft_interval = "ordered_interval",
    soft_graph = "graph_community",
    soft_hamming = "posterior_cluster",
    soft_group_hamming = "posterior_cluster",
    "unknown"
  )
}

kernel_code_length <- function(kernel) {
  as.numeric(kernel$cost)
}

kernel_param_summary <- function(kernel) {
  p <- kernel$params
  if (identical(kernel$type, "safety")) return("unrestricted")
  if (kernel$type %in% c("hard_active", "soft_capacity")) {
    active <- which(p$capacity > 0)
    return(paste0("groups=", paste(active, collapse = "-"), ", cap=", paste(unique(p$capacity[active]), collapse = "-")))
  }
  if (identical(kernel$type, "soft_group_hamming")) {
    return(paste0("groups=", paste(p$center_groups, collapse = "-"), ", rho=", signif(p$rho, 3)))
  }
  if (identical(kernel$type, "soft_hamming")) {
    return(paste0("support_size=", sum(p$center_support > 0), ", rho=", signif(p$rho, 3)))
  }
  if (identical(kernel$type, "soft_interval")) {
    return(paste0("interval=", min(p$interval), "-", max(p$interval), ", cap=", p$capacity, ", rho=", signif(p$rho, 3)))
  }
  if (identical(kernel$type, "soft_graph")) {
    return(paste0("nodes=", paste(head(p$nodes, 10), collapse = "-"), ", rho=", signif(p$rho, 3)))
  }
  ""
}

kernel_summary_table <- function(dictionary, q = NULL, alpha = NULL, top_n = Inf) {
  kernels <- dictionary$kernels
  n <- length(kernels)
  if (is.null(q)) q <- rep(NA_real_, n)
  if (is.null(alpha)) alpha <- rep(NA_real_, n)
  ord <- if (all(is.na(q))) seq_len(n) else order(q, decreasing = TRUE)
  ord <- head(ord, min(length(ord), top_n))
  data.frame(
    rank = seq_along(ord),
    kernel_index = ord,
    name = vapply(kernels[ord], function(k) k$name, character(1)),
    type = vapply(kernels[ord], function(k) k$type, character(1)),
    origin = vapply(kernels[ord], kernel_origin, character(1)),
    q = as.numeric(q[ord]),
    alpha = as.numeric(alpha[ord]),
    code = vapply(kernels[ord], kernel_code_length, numeric(1)),
    params = vapply(kernels[ord], kernel_param_summary, character(1)),
    key = vapply(kernels[ord], function(k) k$key, character(1)),
    stringsAsFactors = FALSE
  )
}

near_duplicate_kernel_keep <- function(W,
                                       base_weights,
                                       alpha_floor = 1e-8,
                                       tol = 1e-8,
                                       always_keep = integer(0)) {
  W <- as.matrix(W)
  weights <- normalize_weights(base_weights)
  G <- normalized_kernel_columns(W, weights, alpha_floor)$G
  keep <- rep(FALSE, ncol(G))
  always_keep <- sort(unique(as.integer(always_keep)))
  always_keep <- always_keep[always_keep >= 1L & always_keep <= ncol(G)]
  keep[always_keep] <- TRUE
  kept <- which(keep)
  for (j in seq_len(ncol(G))) {
    if (keep[j]) next
    if (!length(kept)) {
      keep[j] <- TRUE
      kept <- j
      next
    }
    d <- vapply(kept, function(k) sum(weights * abs(G[, j] - G[, k])), numeric(1))
    if (min(d) > tol) {
      keep[j] <- TRUE
      kept <- c(kept, j)
    }
  }
  keep
}

normalized_kernel_columns <- function(W, base_weights, alpha_floor = 1e-8) {
  W <- as.matrix(W)
  weights <- normalize_weights(base_weights)
  alpha <- as.numeric(crossprod(weights, W))
  alpha_floored <- pmax(alpha, alpha_floor)
  G <- sweep(W, 2, alpha_floored, "/")
  G[!is.finite(G)] <- 0
  list(G = G, alpha = alpha, alpha_floored = alpha_floored, alpha_floor = alpha_floor)
}

normalized_kernel_l1_distance <- function(W_oracle,
                                          W_pool,
                                          base_weights,
                                          alpha_floor = 1e-8) {
  W_oracle <- as.matrix(W_oracle)
  W_pool <- as.matrix(W_pool)
  if (nrow(W_oracle) != nrow(W_pool)) {
    stop("oracle and pool matrices must have the same number of rows")
  }
  weights <- normalize_weights(base_weights)
  G_oracle <- normalized_kernel_columns(W_oracle, weights, alpha_floor)$G
  G_pool <- normalized_kernel_columns(W_pool, weights, alpha_floor)$G
  D <- matrix(0, nrow = ncol(G_oracle), ncol = ncol(G_pool))
  for (j in seq_len(ncol(G_oracle))) {
    D[j, ] <- as.numeric(crossprod(weights, abs(sweep(G_pool, 1, G_oracle[, j], "-"))))
  }
  D
}

kernel_pool_cover_diagnostics <- function(W_oracle,
                                          W_pool,
                                          base_weights,
                                          alpha_floor = 1e-8,
                                          costs_oracle = NULL,
                                          costs_pool = NULL) {
  D <- normalized_kernel_l1_distance(W_oracle, W_pool, base_weights, alpha_floor)
  nearest <- max.col(-D, ties.method = "first")
  nearest_distance <- D[cbind(seq_len(nrow(D)), nearest)]
  out <- data.frame(
    oracle_column = seq_along(nearest),
    nearest_pool_column = nearest,
    normalized_l1 = nearest_distance,
    stringsAsFactors = FALSE
  )
  if (!is.null(costs_oracle) && !is.null(costs_pool)) {
    out$oracle_cost <- as.numeric(costs_oracle)
    out$nearest_pool_cost <- as.numeric(costs_pool)[nearest]
    out$code_excess <- out$nearest_pool_cost - out$oracle_cost
  }
  attr(out, "xi_max") <- max(nearest_distance)
  attr(out, "xi_mean") <- mean(nearest_distance)
  attr(out, "xi_median") <- stats::median(nearest_distance)
  if ("code_excess" %in% names(out)) {
    attr(out, "xi_code_max") <- max(out$code_excess, na.rm = TRUE)
  }
  out
}
