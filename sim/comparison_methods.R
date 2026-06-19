log_prior_bernoulli <- function(supports, theta) {
  supports <- as.matrix(supports)
  size <- rowSums(supports)
  p <- ncol(supports)
  size * log(theta) + (p - size) * log1p(-theta)
}

support_logdet_correlation <- function(supports, X, jitter = 1e-8) {
  supports <- as.matrix(supports)
  corr <- stats::cor(X)
  corr[!is.finite(corr)] <- 0
  diag(corr) <- 1

  vapply(seq_len(nrow(supports)), function(i) {
    active <- which(supports[i, ] == 1L)
    q <- length(active)
    if (q <= 1L) {
      return(0)
    }
    R <- corr[active, active, drop = FALSE]
    out <- tryCatch({
      chol_R <- chol(R + diag(jitter, q))
      2 * sum(log(diag(chol_R)))
    }, error = function(e) -Inf)
    out
  }, numeric(1))
}

fit_bma_from_log_prior <- function(supports,
                                   log_marginal,
                                   log_prior,
                                   theta,
                                   tau2,
                                   a0,
                                   b0,
                                   n) {
  supports <- as.matrix(supports)
  log_unnormalized <- as.numeric(log_marginal) + as.numeric(log_prior)
  log_evidence <- log_sum_exp(log_unnormalized)
  posterior <- exp(log_unnormalized - log_evidence)
  posterior[!is.finite(posterior)] <- 0
  posterior <- posterior / sum(posterior)

  list(
    supports = supports,
    log_marginal = as.numeric(log_marginal),
    log_prior = as.numeric(log_prior),
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
