support_index <- function(supports, target_support) {
  matches <- which(rowSums(abs(sweep(supports, 2, target_support, "-"))) == 0)
  if (length(matches) != 1) {
    stop("Could not identify target support")
  }
  matches
}

summarize_trial <- function(data,
                            exact_fit,
                            repfam_fit,
                            pred_bma,
                            pred_repfam) {
  support_sum <- posterior_summary(exact_fit$posterior)
  family_sum <- posterior_summary(repfam_fit$family_posterior)
  agg_sum <- posterior_summary(repfam_fit$aggregated_support)
  true_support_idx <- support_index(exact_fit$supports, data$true_support)
  oracle_family_idx <- family_index(repfam_fit$families, data$active_groups)
  one_rep_mass <- sum(exact_fit$posterior[repfam_fit$one_rep])
  multi_hit_mass <- 1 - one_rep_mass

  data.frame(
    rho = data$rho,
    support_entropy = support_sum$entropy,
    support_effective_count = support_sum$effective_count,
    support_n95 = support_sum$n_mass,
    family_entropy = family_sum$entropy,
    family_effective_count = family_sum$effective_count,
    family_n95 = family_sum$n_mass,
    aggregated_support_entropy = agg_sum$entropy,
    aggregated_support_n95 = agg_sum$n_mass,
    one_rep_mass = one_rep_mass,
    multi_hit_mass = multi_hit_mass,
    oracle_retained_mass = repfam_fit$alpha[oracle_family_idx],
    oracle_family_posterior = repfam_fit$family_posterior[oracle_family_idx],
    true_support_posterior = exact_fit$posterior[true_support_idx],
    bma_rmse = pred_bma$rmse,
    repfam_rmse = pred_repfam$rmse,
    bma_log_score = pred_bma$mean_log_score,
    repfam_log_score = pred_repfam$mean_log_score
  )
}
