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

family_active_key <- function(active_groups) {
  if (length(active_groups) == 0L) {
    return("empty")
  }
  paste(sort(as.integer(active_groups)), collapse = "-")
}

parse_family_active_key <- function(key) {
  if (identical(key, "empty") || is.na(key) || !nzchar(key)) {
    return(integer(0))
  }
  as.integer(strsplit(key, "-", fixed = TRUE)[[1]])
}

candidate_active_sets_from_supports <- function(supports,
                                                group_id,
                                                weights = NULL,
                                                max_sets = 64L,
                                                include_empty = TRUE) {
  counts <- support_group_counts(supports, group_id)
  active_keys <- apply(counts > 0, 1, function(z) family_active_key(which(z)))
  if (is.null(weights)) {
    weights <- rep(1 / nrow(supports), nrow(supports))
  }
  score <- rowsum(weights, active_keys, reorder = FALSE)[, 1]
  ord <- order(score, decreasing = TRUE)
  keys <- names(score)[ord]
  if (include_empty && !("empty" %in% keys)) {
    keys <- c(keys, "empty")
  }
  keys <- head(unique(keys), max_sets)
  lapply(keys, parse_family_active_key)
}

make_representative_dictionary <- function(group_id,
                                           active_sets,
                                           capacity = 1L,
                                           include_safety = TRUE,
                                           cost_scale = NULL,
                                           name_prefix = "rep") {
  K <- max(group_id)
  if (is.null(cost_scale)) {
    cost_scale <- log(max(K, 2L))
  }
  active_sets <- lapply(active_sets, function(x) sort(unique(as.integer(x))))
  keys <- vapply(active_sets, family_active_key, character(1))
  keep <- !duplicated(keys)
  active_sets <- active_sets[keep]
  keys <- keys[keep]

  rows <- data.frame(
    id = seq_along(active_sets),
    name = paste0(name_prefix, "_", keys),
    type = "representative",
    cost = vapply(active_sets, length, integer(1)) * cost_scale,
    safety = FALSE,
    stringsAsFactors = FALSE
  )
  capacities <- replicate(length(active_sets), as.integer(rep(capacity, K)), simplify = FALSE)

  if (include_safety) {
    rows <- rbind(
      data.frame(
        id = 0L,
        name = "safety_unrestricted",
        type = "safety",
        cost = max(cost_scale * K * 4, 1),
        safety = TRUE,
        stringsAsFactors = FALSE
      ),
      rows
    )
    active_sets <- c(list(seq_len(K)), active_sets)
    capacities <- c(list(as.integer(rep(max(tabulate(group_id), na.rm = TRUE), K))), capacities)
    rows$id <- seq_len(nrow(rows))
  }

  structure(
    list(
      families = rows,
      active_sets = active_sets,
      capacities = capacities,
      group_id = as.integer(group_id)
    ),
    class = "family_dictionary"
  )
}

make_multiresolution_interval_dictionary <- function(p,
                                                     lengths = c(2, 4, 8, 16, 32, 64),
                                                     capacity = 1L,
                                                     max_unions = 1L,
                                                     include_safety = TRUE) {
  intervals <- list()
  for (len in lengths[lengths <= p]) {
    starts <- unique(c(seq(1, p, by = max(1L, floor(len / 2))), p - len + 1L))
    starts <- starts[starts >= 1L & starts + len - 1L <= p]
    for (s in starts) {
      intervals[[length(intervals) + 1L]] <- s:(s + len - 1L)
    }
  }
  active_sets <- lapply(intervals, identity)
  if (max_unions >= 2L && length(intervals) >= 2L) {
    pairs <- utils::combn(seq_along(intervals), 2)
    pair_sets <- lapply(seq_len(ncol(pairs)), function(j) {
      sort(unique(c(intervals[[pairs[1, j]]], intervals[[pairs[2, j]]])))
    })
    active_sets <- c(active_sets, pair_sets)
  }
  group_id <- seq_len(p)
  dict <- make_representative_dictionary(
    group_id = group_id,
    active_sets = active_sets,
    capacity = capacity,
    include_safety = include_safety,
    cost_scale = log(max(p, 2L)),
    name_prefix = "interval"
  )
  dict$families$type[!dict$families$safety] <- "interval"
  dict
}

make_correlation_graph_dictionary <- function(X,
                                              thresholds = c(0.5, 0.6, 0.7, 0.8, 0.9),
                                              capacity = 1L,
                                              include_safety = TRUE,
                                              max_components = 64L) {
  p <- ncol(X)
  group_id <- seq_len(p)
  active_sets <- list()
  names_out <- character(0)
  for (thr in thresholds) {
    component_id <- learn_correlation_groups(X, threshold = thr)
    for (k in seq_len(max(component_id))) {
      members <- which(component_id == k)
      if (!length(members)) next
      key <- paste(members, collapse = "-")
      if (!key %in% names_out) {
        active_sets[[length(active_sets) + 1L]] <- members
        names_out <- c(names_out, key)
      }
    }
  }
  if (!length(active_sets)) {
    stop("no correlation graph dictionaries were created")
  }
  if (length(active_sets) > max_components) {
    active_sets <- active_sets[seq_len(max_components)]
  }
  dict <- make_representative_dictionary(
    group_id = group_id,
    active_sets = active_sets,
    capacity = capacity,
    include_safety = include_safety,
    cost_scale = log(max(p, 2L)),
    name_prefix = "graph"
  )
  dict$families$type[!dict$families$safety] <- "graph_component"
  dict
}

make_adaptive_representative_dictionary <- function(X,
                                                    group_id,
                                                    active_sets,
                                                    u_max = 4L,
                                                    include_safety = TRUE,
                                                    name_prefix = "adaptive") {
  u <- adaptive_capacity_from_groups(X, group_id, u_max = u_max)
  dict <- make_representative_dictionary(
    group_id = group_id,
    active_sets = active_sets,
    capacity = 1L,
    include_safety = include_safety,
    name_prefix = name_prefix
  )
  for (i in seq_along(dict$capacities)) {
    if (!isTRUE(dict$families$safety[i])) {
      dict$capacities[[i]] <- as.integer(u)
    }
  }
  dict$adaptive_capacity <- u
  dict$families$type[!dict$families$safety] <- "adaptive"
  dict
}

make_support_atom_dictionary <- function(supports,
                                         posterior,
                                         top_m = 32L,
                                         include_safety = TRUE,
                                         name_prefix = "support") {
  posterior <- posterior / sum(posterior)
  ord <- order(posterior, decreasing = TRUE)
  keep <- head(ord, min(top_m, length(ord)))
  rows <- data.frame(
    id = seq_along(keep),
    name = paste0(name_prefix, "_", seq_along(keep)),
    type = "support_atom",
    cost = log(max(ncol(supports), 2L)) * pmax(1, rowSums(supports[keep, , drop = FALSE])),
    safety = FALSE,
    stringsAsFactors = FALSE
  )
  atom_supports <- supports[keep, , drop = FALSE]
  if (include_safety) {
    rows <- rbind(
      data.frame(
        id = 0L,
        name = "safety_unrestricted",
        type = "safety",
        cost = max(log(max(ncol(supports), 2L)) * ncol(supports) * 2, 1),
        safety = TRUE,
        stringsAsFactors = FALSE
      ),
      rows
    )
    atom_supports <- rbind(rep(NA_integer_, ncol(supports)), atom_supports)
    rows$id <- seq_len(nrow(rows))
  }
  structure(
    list(
      families = rows,
      atom_supports = atom_supports,
      group_id = seq_len(ncol(supports)),
      active_sets = vector("list", nrow(rows)),
      capacities = vector("list", nrow(rows))
    ),
    class = "family_dictionary"
  )
}

learn_correlation_groups <- function(X, threshold = 0.8) {
  p <- ncol(X)
  C <- abs(stats::cor(X))
  C[is.na(C)] <- 0
  adj <- C >= threshold
  diag(adj) <- TRUE
  group_id <- integer(p)
  current <- 0L
  for (j in seq_len(p)) {
    if (group_id[j] != 0L) next
    current <- current + 1L
    queue <- j
    group_id[j] <- current
    while (length(queue)) {
      node <- queue[1]
      queue <- queue[-1]
      nb <- which(adj[node, ] & group_id == 0L)
      if (length(nb)) {
        group_id[nb] <- current
        queue <- c(queue, nb)
      }
    }
  }
  group_id
}

effective_rank <- function(X_group) {
  if (ncol(X_group) == 0L) {
    return(0)
  }
  vals <- eigen(crossprod(X_group), symmetric = TRUE, only.values = TRUE)$values
  vals <- pmax(vals, 0)
  if (sum(vals) <= 0) {
    return(0)
  }
  sum(vals)^2 / sum(vals^2)
}

adaptive_capacity_from_groups <- function(X, group_id, u_max = 4L) {
  K <- max(group_id)
  u <- integer(K)
  for (k in seq_len(K)) {
    er <- effective_rank(X[, group_id == k, drop = FALSE])
    u[k] <- min(u_max, max(1L, ceiling(er)))
  }
  u
}

family_membership_matrix <- function(supports, dictionary) {
  if (!inherits(dictionary, "family_dictionary")) {
    stop("dictionary must have class family_dictionary")
  }
  if (!is.null(dictionary$atom_supports)) {
    n <- nrow(supports)
    M <- matrix(FALSE, nrow = n, ncol = nrow(dictionary$families))
    colnames(M) <- dictionary$families$name
    for (m in seq_len(ncol(M))) {
      if (isTRUE(dictionary$families$safety[m])) {
        M[, m] <- TRUE
      } else {
        atom <- dictionary$atom_supports[m, ]
        M[, m] <- rowSums(abs(sweep(supports, 2, atom, "-"))) == 0
      }
    }
    return(M)
  }
  group_id <- dictionary$group_id
  counts <- support_group_counts(supports, group_id)
  active <- counts > 0
  n <- nrow(supports)
  M <- matrix(FALSE, nrow = n, ncol = nrow(dictionary$families))
  colnames(M) <- dictionary$families$name

  for (m in seq_len(ncol(M))) {
    if (isTRUE(dictionary$families$safety[m])) {
      M[, m] <- TRUE
      next
    }
    S <- dictionary$active_sets[[m]]
    cap <- dictionary$capacities[[m]]
    outside <- if (length(S) == max(group_id)) {
      rep(FALSE, n)
    } else {
      rowSums(active[, setdiff(seq_len(max(group_id)), S), drop = FALSE]) > 0
    }
    capacity_ok <- apply(sweep(counts, 2, cap, "<="), 1, all)
    M[, m] <- !outside & capacity_ok
  }
  M
}
