source(file.path("sim", "src", "family_dictionary.R"))

family_index <- function(families, active_groups) {
  active_groups <- sort(unique(as.integer(active_groups)))
  pattern <- integer(ncol(families))
  if (length(active_groups)) {
    pattern[active_groups] <- 1L
  }
  hit <- which(rowSums(abs(sweep(families, 2, pattern, "-"))) == 0)
  if (!length(hit)) {
    stop("Requested family is not present in the family table")
  }
  hit[1]
}

fit_repfam_bma <- function(exact_fit, group_id, lambda = 1) {
  supports <- as.matrix(exact_fit$supports)
  weights <- exact_fit$posterior / sum(exact_fit$posterior)
  K <- max(group_id)
  families <- enumerate_supports(K)
  colnames(families) <- paste0("G", seq_len(K))

  counts <- support_group_counts(supports, group_id)
  active <- counts > 0
  membership <- matrix(FALSE, nrow = nrow(supports), ncol = nrow(families))
  for (j in seq_len(nrow(families))) {
    membership[, j] <- rowSums(abs(sweep(active * 1L, 2, families[j, ], "-"))) == 0
  }

  alpha <- as.numeric(crossprod(weights, membership))
  size <- rowSums(families)
  log_weight <- ifelse(alpha > 0, log(alpha), -Inf) - lambda * size * log(max(K, 2L))
  family_posterior <- exp(log_weight - log_sum_exp(log_weight))

  list(
    families = families,
    alpha = alpha,
    family_posterior = family_posterior,
    membership = membership,
    lambda = lambda,
    group_id = as.integer(group_id)
  )
}
