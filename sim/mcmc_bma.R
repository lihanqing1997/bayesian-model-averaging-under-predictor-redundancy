gamma_key <- function(gamma) {
  active <- which(gamma == 1L)
  if (length(active) == 0L) {
    return("0")
  }
  paste(active, collapse = ",")
}

sample_one <- function(x) {
  if (length(x) == 0L) {
    stop("Cannot sample from an empty vector")
  }
  x[sample.int(length(x), 1L)]
}

screen_start_support <- function(X, y, start_size) {
  p <- ncol(X)
  start_size <- min(start_size, p)
  scores <- abs(as.numeric(crossprod(X, y)))
  gamma <- integer(p)
  if (start_size > 0L) {
    gamma[order(scores, decreasing = TRUE)[seq_len(start_size)]] <- 1L
  }
  gamma
}

random_start_support <- function(X, y, start_size, chain_id) {
  p <- ncol(X)
  scores <- abs(as.numeric(crossprod(X, y)))
  pool_size <- min(p, max(3 * start_size, start_size + 5))
  pool <- order(scores, decreasing = TRUE)[seq_len(pool_size)]
  gamma <- integer(p)

  if (chain_id == 1L) {
    return(screen_start_support(X, y, start_size))
  }

  if (start_size > 0L) {
    gamma[pool[sample.int(length(pool), min(start_size, length(pool)))] ] <- 1L
  }
  gamma
}

logdet_corr_gamma <- function(gamma, corr, jitter = 1e-8) {
  active <- which(gamma == 1L)
  q <- length(active)
  if (q <= 1L) {
    return(0)
  }

  Rg <- corr[active, active, drop = FALSE] + diag(jitter, q)
  out <- tryCatch({
    chol_Rg <- chol(Rg)
    2 * sum(log(diag(chol_Rg)))
  }, error = function(e) -Inf)

  out
}

make_group_aware_proposal <- function(gamma,
                                      group_id,
                                      max_size,
                                      proposal_probs) {
  proposal_probs <- proposal_probs / sum(proposal_probs)
  type <- sample(names(proposal_probs), 1L, prob = proposal_probs)
  proposal <- gamma
  p <- length(gamma)
  active <- which(gamma == 1L)
  inactive <- which(gamma == 0L)
  counts <- tabulate(group_id[active], nbins = max(group_id))
  valid <- TRUE

  if (type == "flip_variable") {
    j <- sample.int(p, 1L)
    proposal[j] <- 1L - proposal[j]
    valid <- sum(proposal) <= max_size
  } else if (type == "swap_variable") {
    valid <- length(active) > 0L && length(inactive) > 0L
    if (valid) {
      drop_j <- sample_one(active)
      add_j <- sample_one(inactive)
      proposal[drop_j] <- 0L
      proposal[add_j] <- 1L
    }
  } else if (type == "swap_within_group") {
    eligible <- which(counts > 0L & counts < tabulate(group_id, nbins = max(group_id)))
    valid <- length(eligible) > 0L
    if (valid) {
      g <- sample_one(eligible)
      group_columns <- which(group_id == g)
      drop_j <- sample_one(intersect(group_columns, active))
      add_j <- sample_one(intersect(group_columns, inactive))
      proposal[drop_j] <- 0L
      proposal[add_j] <- 1L
    }
  } else if (type == "swap_singleton_group") {
    singleton_groups <- which(counts == 1L)
    inactive_groups <- which(counts == 0L)
    valid <- length(singleton_groups) > 0L && length(inactive_groups) > 0L
    if (valid) {
      drop_group <- sample_one(singleton_groups)
      add_group <- sample_one(inactive_groups)
      drop_j <- intersect(which(group_id == drop_group), active)
      add_j <- sample_one(which(group_id == add_group))
      proposal[drop_j] <- 0L
      proposal[add_j] <- 1L
    }
  } else {
    stop(sprintf("Unknown proposal type: %s", type))
  }

  list(gamma = proposal, type = type, valid = valid)
}

combine_sampled_supports <- function(samples) {
  if (nrow(samples) == 0L) {
    stop("No MCMC samples were saved")
  }

  sample_keys <- apply(samples, 1, gamma_key)
  tab <- table(sample_keys)
  unique_keys <- names(tab)
  p <- ncol(samples)
  unique_supports <- matrix(0L, nrow = length(unique_keys), ncol = p)
  for (i in seq_along(unique_keys)) {
    key <- unique_keys[i]
    if (key != "0") {
      unique_supports[i, as.integer(strsplit(key, ",", fixed = TRUE)[[1]])] <- 1L
    }
  }

  weights <- as.numeric(tab) / sum(tab)
  order_id <- order(weights, decreasing = TRUE)
  list(
    supports = unique_supports[order_id, , drop = FALSE],
    posterior = weights[order_id]
  )
}

group_pip_from_samples <- function(samples, group_id) {
  counts <- support_group_counts(samples, group_id)
  colMeans(counts > 0)
}

fit_mcmc_bma <- function(X,
                         y,
                         group_id = seq_len(ncol(X)),
                         theta = 0.04,
                         tau2 = 4,
                         a0 = 1,
                         b0 = 1,
                         prior_type = c("unrestricted", "dilution", "dpp"),
                         n_iter = 3000,
                         burn = 1000,
                         thin = 5,
                         max_size = 12,
                         start_size = 6,
                         n_chains = 1,
                         proposal_probs = c(
                           flip_variable = 0.25,
                           swap_variable = 0.25,
                           swap_within_group = 0.35,
                           swap_singleton_group = 0.15
                         ),
                         seed = NULL) {
  prior_type <- match.arg(prior_type)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  p <- ncol(X)
  corr <- stats::cor(X)
  corr[!is.finite(corr)] <- 0
  diag(corr) <- 1

  cache <- new.env(parent = emptyenv())
  log_post_gamma <- function(gamma) {
    if (sum(gamma) > max_size) {
      return(-Inf)
    }

    key <- gamma_key(gamma)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }

    size <- sum(gamma)
    log_prior <- size * log(theta) + (p - size) * log1p(-theta)
    extra <- switch(
      prior_type,
      unrestricted = 0,
      dilution = 0.5 * logdet_corr_gamma(gamma, corr),
      dpp = logdet_corr_gamma(gamma, corr)
    )
    value <- log_marginal_support(
      gamma = gamma,
      X = X,
      y = y,
      tau2 = tau2,
      a0 = a0,
      b0 = b0
    ) + log_prior + extra

    assign(key, value, envir = cache)
    value
  }

  n_saved <- floor((n_iter - burn) / thin)
  all_samples <- matrix(0L, nrow = 0L, ncol = p)
  samples_by_chain <- vector("list", n_chains)
  trace_rows <- list()
  chain_rows <- list()
  proposal_rows <- list()
  group_pips <- matrix(NA_real_, nrow = n_chains, ncol = max(group_id))

  for (chain in seq_len(n_chains)) {
    if (!is.null(seed)) {
      set.seed(seed + 997L * chain)
    }

    gamma <- random_start_support(X, y, start_size, chain)
    if (sum(gamma) > max_size) {
      keep <- which(gamma == 1L)[seq_len(max_size)]
      gamma <- integer(p)
      gamma[keep] <- 1L
    }

    current <- log_post_gamma(gamma)
    samples <- matrix(0L, nrow = n_saved, ncol = p)
    save_index <- 0L
    accepted <- 0L
    proposal_count <- setNames(integer(length(proposal_probs)), names(proposal_probs))
    accepted_count <- setNames(integer(length(proposal_probs)), names(proposal_probs))
    trace_size <- integer(n_saved)
    trace_log_post <- numeric(n_saved)

    for (iter in seq_len(n_iter)) {
      proposal <- make_group_aware_proposal(
        gamma = gamma,
        group_id = group_id,
        max_size = max_size,
        proposal_probs = proposal_probs
      )
      proposal_count[proposal$type] <- proposal_count[proposal$type] + 1L

      if (proposal$valid) {
        proposed <- log_post_gamma(proposal$gamma)
        if (is.finite(proposed) && log(runif(1)) < proposed - current) {
          gamma <- proposal$gamma
          current <- proposed
          accepted <- accepted + 1L
          accepted_count[proposal$type] <- accepted_count[proposal$type] + 1L
        }
      }

      if (iter > burn && ((iter - burn) %% thin == 0L)) {
        save_index <- save_index + 1L
        samples[save_index, ] <- gamma
        trace_size[save_index] <- sum(gamma)
        trace_log_post[save_index] <- current
      }
    }

    samples <- samples[seq_len(save_index), , drop = FALSE]
    samples_by_chain[[chain]] <- samples
    all_samples <- rbind(all_samples, samples)
    group_pips[chain, ] <- group_pip_from_samples(samples, group_id)

    trace_rows[[chain]] <- data.frame(
      chain = chain,
      iteration = seq_len(save_index),
      size = trace_size[seq_len(save_index)],
      log_posterior = trace_log_post[seq_len(save_index)]
    )
    chain_rows[[chain]] <- data.frame(
      chain = chain,
      accept_rate = accepted / n_iter,
      unique_count = length(unique(apply(samples, 1, gamma_key))),
      mean_size = mean(rowSums(samples)),
      log_posterior_mean = mean(trace_log_post[seq_len(save_index)]),
      log_posterior_sd = stats::sd(trace_log_post[seq_len(save_index)])
    )
    proposal_rows[[chain]] <- data.frame(
      chain = chain,
      proposal_type = names(proposal_count),
      proposed = as.numeric(proposal_count),
      accepted = as.numeric(accepted_count),
      accept_rate = ifelse(proposal_count > 0,
                           accepted_count / proposal_count,
                           NA_real_)
    )
  }

  combined <- combine_sampled_supports(all_samples)
  unique_supports <- combined$supports
  weights <- combined$posterior
  chain_diagnostics <- do.call(rbind, chain_rows)
  proposal_diagnostics <- do.call(rbind, proposal_rows)
  traces <- do.call(rbind, trace_rows)
  group_pip_range <- apply(group_pips, 2, function(z) max(z) - min(z))

  list(
    supports = unique_supports,
    posterior = weights,
    samples = all_samples,
    samples_by_chain = samples_by_chain,
    chain_diagnostics = chain_diagnostics,
    proposal_diagnostics = proposal_diagnostics,
    traces = traces,
    group_pips = group_pips,
    group_pip_max_range = max(group_pip_range),
    group_pip_mean_range = mean(group_pip_range),
    accept_rate = mean(chain_diagnostics$accept_rate),
    unique_count = length(weights),
    mean_size = mean(rowSums(all_samples)),
    log_posterior_mean = mean(chain_diagnostics$log_posterior_mean),
    log_posterior_sd = mean(chain_diagnostics$log_posterior_sd),
    theta = theta,
    tau2 = tau2,
    a0 = a0,
    b0 = b0,
    n = nrow(X),
    p = p,
    prior_type = prior_type,
    max_size = max_size,
    n_iter = n_iter,
    burn = burn,
    thin = thin,
    n_chains = n_chains
  )
}

proposal_rate <- function(fit, type) {
  diagnostics <- fit$proposal_diagnostics
  rows <- diagnostics[diagnostics$proposal_type == type, ]
  if (nrow(rows) == 0L) {
    return(NA_real_)
  }
  sum(rows$accepted) / sum(rows$proposed)
}

safe_true_support_probability <- function(supports, weights, true_support) {
  matches <- which(rowSums(abs(sweep(supports, 2, true_support, "-"))) == 0)
  if (length(matches) == 0L) {
    return(0)
  }
  sum(weights[matches])
}

sampled_oracle_retained_mass <- function(samples, group_id, active_groups) {
  counts <- support_group_counts(samples, group_id)
  active <- counts > 0
  one_rep <- apply(counts, 1, max) <= 1
  target <- rep(FALSE, max(group_id))
  target[active_groups] <- TRUE
  outside <- rowSums(active[, !target, drop = FALSE]) > 0
  mean(one_rep & !outside)
}

sampled_support_summary <- function(data, fit, pred, method, fit_seconds) {
  support_sum <- posterior_summary(fit$posterior)
  counts <- support_group_counts(fit$supports, data$group_id)
  group_pip <- as.numeric(crossprod(fit$posterior / sum(fit$posterior), counts > 0))
  active_group <- seq_len(data$K) %in% data$active_groups
  top_metrics <- top_support_metrics(fit$supports, fit$posterior, data$true_support)

  data.frame(
    rho = data$rho,
    method = method,
    lambda = NA_real_,
    rmse = pred$rmse,
    log_score = pred$mean_log_score,
    support_n95 = support_sum$n_mass,
    support_effective_count = support_sum$effective_count,
    family_n95 = NA_real_,
    family_effective_count = NA_real_,
    retained_mass = NA_real_,
    true_support_posterior = safe_true_support_probability(
      fit$supports,
      fit$posterior,
      data$true_support
    ),
    active_group_pip_mean = mean(group_pip[active_group]),
    inactive_group_pip_max = max(group_pip[!active_group]),
    top_support_size = top_metrics$top_support_size,
    top_support_f1 = top_metrics$top_support_f1,
    fit_seconds = fit_seconds,
    accept_rate = fit$accept_rate,
    flip_accept_rate = proposal_rate(fit, "flip_variable"),
    swap_accept_rate = proposal_rate(fit, "swap_variable"),
    within_group_accept_rate = proposal_rate(fit, "swap_within_group"),
    singleton_group_accept_rate = proposal_rate(fit, "swap_singleton_group"),
    chain_count = fit$n_chains,
    group_pip_max_range = fit$group_pip_max_range,
    group_pip_mean_range = fit$group_pip_mean_range,
    log_posterior_sd = fit$log_posterior_sd,
    unique_count = fit$unique_count,
    mean_size = fit$mean_size
  )
}

candidate_families_from_samples <- function(samples,
                                            group_id,
                                            top_group_count = 10,
                                            max_family_size = 8) {
  counts <- support_group_counts(samples, group_id)
  active_matrix <- counts > 0
  group_pip <- colMeans(active_matrix)
  K <- ncol(active_matrix)
  top_groups <- order(group_pip, decreasing = TRUE)[seq_len(min(top_group_count, K))]

  candidate <- matrix(0L, nrow = 1L, ncol = K)
  G <- length(top_groups)
  ids <- 0:(2^G - 1)
  for (id in ids) {
    mask <- as.integer(intToBits(id))[seq_len(G)]
    if (sum(mask) <= max_family_size) {
      row <- integer(K)
      row[top_groups[mask == 1L]] <- 1L
      candidate <- rbind(candidate, row)
    }
  }

  observed <- unique(active_matrix * 1L)
  observed <- observed[rowSums(observed) <= max_family_size, , drop = FALSE]
  candidate <- unique(rbind(candidate, observed))
  colnames(candidate) <- paste0("g", seq_len(K))
  candidate
}

fit_repfam_from_mcmc <- function(ubma_fit,
                                 group_id,
                                 lambda = 0,
                                 top_group_count = 10,
                                 max_family_size = 8) {
  samples <- ubma_fit$samples
  K <- max(group_id)
  families <- candidate_families_from_samples(
    samples = samples,
    group_id = group_id,
    top_group_count = top_group_count,
    max_family_size = max_family_size
  )

  sample_counts <- support_group_counts(samples, group_id)
  sample_active <- sample_counts > 0
  sample_one_rep <- apply(sample_counts, 1, max) <= 1

  alpha <- numeric(nrow(families))
  for (f in seq_len(nrow(families))) {
    S <- as.logical(families[f, ])
    outside <- if (all(S)) {
      rep(FALSE, nrow(samples))
    } else {
      rowSums(sample_active[, !S, drop = FALSE]) > 0
    }
    alpha[f] <- mean(sample_one_rep & !outside)
  }

  keep <- alpha > 0
  families <- families[keep, , drop = FALSE]
  alpha <- alpha[keep]
  family_size <- rowSums(families)
  log_weight <- log(alpha) - lambda * family_size * log(K)
  family_posterior <- exp(log_weight - log_sum_exp(log_weight))

  unique_counts <- support_group_counts(ubma_fit$supports, group_id)
  unique_active <- unique_counts > 0
  unique_one_rep <- apply(unique_counts, 1, max) <= 1
  aggregated <- numeric(length(ubma_fit$posterior))

  for (f in seq_len(nrow(families))) {
    S <- as.logical(families[f, ])
    outside <- if (all(S)) {
      rep(FALSE, nrow(ubma_fit$supports))
    } else {
      rowSums(unique_active[, !S, drop = FALSE]) > 0
    }
    members <- unique_one_rep & !outside
    member_mass <- sum(ubma_fit$posterior[members])
    if (member_mass > 0) {
      aggregated[members] <- aggregated[members] +
        family_posterior[f] * ubma_fit$posterior[members] / member_mass
    }
  }

  aggregated <- aggregated / sum(aggregated)

  list(
    supports = ubma_fit$supports,
    aggregated_support = aggregated,
    families = families,
    alpha = alpha,
    family_posterior = family_posterior,
    lambda = lambda,
    K = K
  )
}

sampled_repfam_summary <- function(data, ubma_fit, repfam_fit, pred, method, fit_seconds) {
  support_sum <- posterior_summary(repfam_fit$aggregated_support)
  family_sum <- posterior_summary(repfam_fit$family_posterior)
  counts <- support_group_counts(ubma_fit$supports, data$group_id)
  group_pip <- as.numeric(crossprod(
    repfam_fit$aggregated_support / sum(repfam_fit$aggregated_support),
    counts > 0
  ))
  active_group <- seq_len(data$K) %in% data$active_groups
  top_metrics <- top_support_metrics(
    ubma_fit$supports,
    repfam_fit$aggregated_support,
    data$true_support
  )
  retained_mass <- sampled_oracle_retained_mass(
    samples = ubma_fit$samples,
    group_id = data$group_id,
    active_groups = data$active_groups
  )

  data.frame(
    rho = data$rho,
    method = method,
    lambda = repfam_fit$lambda,
    rmse = pred$rmse,
    log_score = pred$mean_log_score,
    support_n95 = support_sum$n_mass,
    support_effective_count = support_sum$effective_count,
    family_n95 = family_sum$n_mass,
    family_effective_count = family_sum$effective_count,
    retained_mass = retained_mass,
    true_support_posterior = safe_true_support_probability(
      ubma_fit$supports,
      repfam_fit$aggregated_support,
      data$true_support
    ),
    active_group_pip_mean = mean(group_pip[active_group]),
    inactive_group_pip_max = max(group_pip[!active_group]),
    top_support_size = top_metrics$top_support_size,
    top_support_f1 = top_metrics$top_support_f1,
    fit_seconds = fit_seconds,
    accept_rate = NA_real_,
    flip_accept_rate = NA_real_,
    swap_accept_rate = NA_real_,
    within_group_accept_rate = NA_real_,
    singleton_group_accept_rate = NA_real_,
    chain_count = ubma_fit$n_chains,
    group_pip_max_range = ubma_fit$group_pip_max_range,
    group_pip_mean_range = ubma_fit$group_pip_mean_range,
    log_posterior_sd = ubma_fit$log_posterior_sd,
    unique_count = length(repfam_fit$aggregated_support),
    mean_size = sum(rowSums(ubma_fit$supports) * repfam_fit$aggregated_support)
  )
}

fit_penalized_baselines <- function(data, seed) {
  rows <- list()

  if (requireNamespace("glmnet", quietly = TRUE)) {
    glmnet_methods <- list(lasso = 1, elastic_net = 0.5)
    for (name in names(glmnet_methods)) {
      set.seed(seed)
      elapsed <- system.time({
        fit <- glmnet::cv.glmnet(
          x = data$X_train,
          y = data$y_train,
          family = "gaussian",
          alpha = glmnet_methods[[name]],
          nfolds = 5,
          intercept = FALSE,
          standardize = FALSE
        )
        pred <- as.numeric(stats::predict(fit, newx = data$X_test, s = "lambda.min"))
      })

      coef_vec <- as.numeric(stats::coef(fit, s = "lambda.min"))[-1]
      rows[[name]] <- penalized_summary_row(
        data = data,
        method = name,
        pred = pred,
        selected = as.integer(abs(coef_vec) > 1e-8),
        fit_seconds = unname(elapsed["elapsed"])
      )
    }
  }

  if (requireNamespace("grpreg", quietly = TRUE)) {
    set.seed(seed + 10)
    elapsed <- try(system.time({
      fit <- grpreg::cv.grpreg(
        X = data$X_train,
        y = data$y_train,
        group = data$group_id,
        penalty = "grLasso",
        family = "gaussian",
        nfolds = 5,
        seed = seed + 10
      )
      pred <- as.numeric(stats::predict(
        fit,
        X = data$X_test,
        type = "response",
        lambda = fit$lambda.min
      ))
      coef_vec <- as.numeric(stats::coef(fit, lambda = fit$lambda.min))[-1]
    }), silent = TRUE)

    if (!inherits(elapsed, "try-error")) {
      rows[["group_lasso"]] <- penalized_summary_row(
        data = data,
        method = "group_lasso",
        pred = pred,
        selected = as.integer(abs(coef_vec) > 1e-8),
        fit_seconds = unname(elapsed["elapsed"])
      )
    }
  }

  if (length(rows) == 0L) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

penalized_summary_row <- function(data, method, pred, selected, fit_seconds) {
  tp <- sum(selected == 1 & data$true_support == 1)
  fp <- sum(selected == 1 & data$true_support == 0)
  fn <- sum(selected == 0 & data$true_support == 1)
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
    NA_real_
  } else {
    2 * precision * recall / (precision + recall)
  }

  selected_groups <- as.numeric(tapply(selected, data$group_id, max))
  active_group <- seq_len(data$K) %in% data$active_groups

  data.frame(
    rho = data$rho,
    method = method,
    lambda = NA_real_,
    rmse = sqrt(mean((data$y_test - pred)^2)),
    log_score = NA_real_,
    support_n95 = NA_real_,
    support_effective_count = NA_real_,
    family_n95 = NA_real_,
    family_effective_count = NA_real_,
    retained_mass = NA_real_,
    true_support_posterior = NA_real_,
    active_group_pip_mean = mean(selected_groups[active_group]),
    inactive_group_pip_max = max(selected_groups[!active_group]),
    top_support_size = sum(selected),
    top_support_f1 = f1,
    fit_seconds = fit_seconds,
    accept_rate = NA_real_,
    flip_accept_rate = NA_real_,
    swap_accept_rate = NA_real_,
    within_group_accept_rate = NA_real_,
    singleton_group_accept_rate = NA_real_,
    chain_count = NA_real_,
    group_pip_max_range = NA_real_,
    group_pip_mean_range = NA_real_,
    log_posterior_sd = NA_real_,
    unique_count = NA_real_,
    mean_size = sum(selected)
  )
}
