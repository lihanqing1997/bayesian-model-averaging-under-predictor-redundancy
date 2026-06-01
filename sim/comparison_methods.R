log_prior_bernoulli <- function(supports, theta) {
  p <- ncol(supports)
  size <- rowSums(supports)
  size * log(theta) + (p - size) * log1p(-theta)
}

support_logdet_correlation <- function(supports, X, jitter = 1e-8) {
  R <- stats::cor(X)
  R[!is.finite(R)] <- 0
  diag(R) <- 1

  out <- numeric(nrow(supports))
  for (i in seq_len(nrow(supports))) {
    active <- which(supports[i, ] == 1)
    q <- length(active)
    if (q <= 1) {
      out[i] <- 0
    } else {
      Rg <- R[active, active, drop = FALSE] + diag(jitter, q)
      chol_Rg <- safe_chol(Rg)
      out[i] <- 2 * sum(log(diag(chol_Rg)))
    }
  }

  out
}

fit_bma_from_log_prior <- function(supports,
                                   log_marginal,
                                   log_prior,
                                   theta,
                                   tau2,
                                   a0,
                                   b0,
                                   n) {
  log_unnormalized <- log_marginal + log_prior
  log_evidence <- log_sum_exp(log_unnormalized)
  posterior <- exp(log_unnormalized - log_evidence)

  list(
    supports = supports,
    log_marginal = log_marginal,
    log_prior = log_prior,
    log_posterior = log(posterior),
    posterior = posterior,
    log_evidence = log_evidence,
    theta = theta,
    tau2 = tau2,
    a0 = a0,
    b0 = b0,
    n = n,
    p = ncol(supports)
  )
}

fit_bma_variants <- function(X,
                             y,
                             supports,
                             theta = 0.12,
                             tau2 = 4,
                             a0 = 1,
                             b0 = 1) {
  log_marginal <- numeric(nrow(supports))
  for (i in seq_len(nrow(supports))) {
    log_marginal[i] <- log_marginal_support(
      gamma = supports[i, ],
      X = X,
      y = y,
      tau2 = tau2,
      a0 = a0,
      b0 = b0
    )
  }

  log_prior_base <- log_prior_bernoulli(supports, theta)
  log_det <- support_logdet_correlation(supports, X)

  list(
    unrestricted_bma = fit_bma_from_log_prior(
      supports = supports,
      log_marginal = log_marginal,
      log_prior = log_prior_base,
      theta = theta,
      tau2 = tau2,
      a0 = a0,
      b0 = b0,
      n = nrow(X)
    ),
    dilution_bma = fit_bma_from_log_prior(
      supports = supports,
      log_marginal = log_marginal,
      log_prior = log_prior_base + 0.5 * log_det,
      theta = theta,
      tau2 = tau2,
      a0 = a0,
      b0 = b0,
      n = nrow(X)
    ),
    dpp_bma = fit_bma_from_log_prior(
      supports = supports,
      log_marginal = log_marginal,
      log_prior = log_prior_base + log_det,
      theta = theta,
      tau2 = tau2,
      a0 = a0,
      b0 = b0,
      n = nrow(X)
    )
  )
}

group_posterior_inclusion <- function(supports, weights, group_id) {
  counts <- support_group_counts(supports, group_id)
  as.numeric(crossprod(weights / sum(weights), counts > 0))
}

oracle_group_pattern_mass <- function(supports, weights, group_id, active_groups) {
  counts <- support_group_counts(supports, group_id)
  pattern <- counts > 0
  target <- rep(FALSE, max(group_id))
  target[active_groups] <- TRUE
  sum(weights[rowSums(abs(sweep(pattern, 2, target, "-"))) == 0])
}

top_support_metrics <- function(supports, weights, true_support) {
  top <- supports[which.max(weights), ]
  tp <- sum(top == 1 & true_support == 1)
  fp <- sum(top == 1 & true_support == 0)
  fn <- sum(top == 0 & true_support == 1)
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
    NA_real_
  } else {
    2 * precision * recall / (precision + recall)
  }

  list(
    top_support_size = sum(top),
    top_support_f1 = f1
  )
}

summarize_bma_method <- function(data, fit, pred, method, fit_seconds = NA_real_) {
  support_sum <- posterior_summary(fit$posterior)
  true_support_idx <- support_index(fit$supports, data$true_support)
  group_pip <- group_posterior_inclusion(fit$supports, fit$posterior, data$group_id)
  top_metrics <- top_support_metrics(fit$supports, fit$posterior, data$true_support)
  active_group <- seq_len(data$K) %in% data$active_groups

  data.frame(
    rho = data$rho,
    method = method,
    lambda = NA_real_,
    rmse = pred$rmse,
    log_score = pred$mean_log_score,
    support_entropy = support_sum$entropy,
    support_effective_count = support_sum$effective_count,
    support_n95 = support_sum$n_mass,
    family_entropy = NA_real_,
    family_effective_count = NA_real_,
    family_n95 = NA_real_,
    one_rep_mass = NA_real_,
    multi_hit_mass = NA_real_,
    oracle_retained_mass = NA_real_,
    oracle_family_posterior = NA_real_,
    true_support_posterior = fit$posterior[true_support_idx],
    oracle_group_pattern_mass = oracle_group_pattern_mass(
      fit$supports,
      fit$posterior,
      data$group_id,
      data$active_groups
    ),
    active_group_pip_mean = mean(group_pip[active_group]),
    inactive_group_pip_max = max(group_pip[!active_group]),
    top_support_size = top_metrics$top_support_size,
    top_support_f1 = top_metrics$top_support_f1,
    fit_seconds = fit_seconds
  )
}

summarize_repfam_method <- function(data,
                                    exact_fit,
                                    repfam_fit,
                                    pred,
                                    method,
                                    fit_seconds = NA_real_) {
  support_sum <- posterior_summary(repfam_fit$aggregated_support)
  family_sum <- posterior_summary(repfam_fit$family_posterior)
  true_support_idx <- support_index(exact_fit$supports, data$true_support)
  oracle_family_idx <- family_index(repfam_fit$families, data$active_groups)
  group_pip <- group_posterior_inclusion(
    exact_fit$supports,
    repfam_fit$aggregated_support,
    data$group_id
  )
  top_metrics <- top_support_metrics(
    exact_fit$supports,
    repfam_fit$aggregated_support,
    data$true_support
  )
  active_group <- seq_len(data$K) %in% data$active_groups

  data.frame(
    rho = data$rho,
    method = method,
    lambda = repfam_fit$lambda,
    rmse = pred$rmse,
    log_score = pred$mean_log_score,
    support_entropy = support_sum$entropy,
    support_effective_count = support_sum$effective_count,
    support_n95 = support_sum$n_mass,
    family_entropy = family_sum$entropy,
    family_effective_count = family_sum$effective_count,
    family_n95 = family_sum$n_mass,
    one_rep_mass = sum(exact_fit$posterior[repfam_fit$one_rep]),
    multi_hit_mass = 1 - sum(exact_fit$posterior[repfam_fit$one_rep]),
    oracle_retained_mass = repfam_fit$alpha[oracle_family_idx],
    oracle_family_posterior = repfam_fit$family_posterior[oracle_family_idx],
    true_support_posterior = repfam_fit$aggregated_support[true_support_idx],
    oracle_group_pattern_mass = oracle_group_pattern_mass(
      exact_fit$supports,
      repfam_fit$aggregated_support,
      data$group_id,
      data$active_groups
    ),
    active_group_pip_mean = mean(group_pip[active_group]),
    inactive_group_pip_max = max(group_pip[!active_group]),
    top_support_size = top_metrics$top_support_size,
    top_support_f1 = top_metrics$top_support_f1,
    fit_seconds = fit_seconds
  )
}

fit_glmnet_baselines <- function(data, seed) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    return(data.frame())
  }

  methods <- list(
    lasso = 1,
    elastic_net = 0.5
  )
  rows <- list()

  for (name in names(methods)) {
    set.seed(seed)
    elapsed <- system.time({
      fit <- glmnet::cv.glmnet(
        x = data$X_train,
        y = data$y_train,
        family = "gaussian",
        alpha = methods[[name]],
        nfolds = 5,
        intercept = FALSE,
        standardize = FALSE
      )
      pred <- as.numeric(stats::predict(fit, newx = data$X_test, s = "lambda.min"))
    })

    coef_vec <- as.numeric(stats::coef(fit, s = "lambda.min"))[-1]
    selected <- as.integer(abs(coef_vec) > 1e-8)
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
    selected_groups <- tapply(selected, data$group_id, max)
    active_group <- seq_len(data$K) %in% data$active_groups

    rows[[name]] <- data.frame(
      rho = data$rho,
      method = name,
      lambda = NA_real_,
      rmse = sqrt(mean((data$y_test - pred)^2)),
      log_score = NA_real_,
      support_entropy = NA_real_,
      support_effective_count = NA_real_,
      support_n95 = NA_real_,
      family_entropy = NA_real_,
      family_effective_count = NA_real_,
      family_n95 = NA_real_,
      one_rep_mass = NA_real_,
      multi_hit_mass = NA_real_,
      oracle_retained_mass = NA_real_,
      oracle_family_posterior = NA_real_,
      true_support_posterior = NA_real_,
      oracle_group_pattern_mass = NA_real_,
      active_group_pip_mean = mean(selected_groups[active_group]),
      inactive_group_pip_max = max(selected_groups[!active_group]),
      top_support_size = sum(selected),
      top_support_f1 = f1,
      fit_seconds = unname(elapsed["elapsed"])
    )
  }

  do.call(rbind, rows)
}

summarize_repeated_results <- function(results) {
  results$lambda_for_summary <- ifelse(is.na(results$lambda), -1, results$lambda)
  out <- aggregate(
    cbind(
      rmse,
      log_score,
      support_n95,
      family_n95,
      support_effective_count,
      family_effective_count,
      one_rep_mass,
      multi_hit_mass,
      oracle_retained_mass,
      oracle_family_posterior,
      true_support_posterior,
      oracle_group_pattern_mass,
      active_group_pip_mean,
      inactive_group_pip_max,
      top_support_size,
      top_support_f1,
      fit_seconds
    ) ~ rho + method + lambda_for_summary,
    data = results,
    FUN = function(x) mean(x, na.rm = TRUE),
    na.action = na.pass
  )
  names(out)[names(out) == "lambda_for_summary"] <- "lambda"
  out$lambda[out$lambda < 0] <- NA_real_
  out[] <- lapply(out, function(column) {
    if (is.numeric(column)) {
      column[is.nan(column)] <- NA_real_
    }
    column
  })
  out
}
