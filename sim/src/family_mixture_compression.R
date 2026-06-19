source(file.path("sim", "src", "family_dictionary.R"))

normalize_weights <- function(w) {
  w <- as.numeric(w)
  if (any(!is.finite(w)) || sum(w) <= 0) {
    stop("weights must be finite and have positive sum")
  }
  w / sum(w)
}

estimate_family_alpha <- function(membership, base_weights, alpha_floor = 1e-8) {
  w <- normalize_weights(base_weights)
  alpha <- as.numeric(crossprod(w, membership))
  list(
    alpha = alpha,
    alpha_truncated = pmax(alpha, alpha_floor),
    alpha_floor = alpha_floor,
    floor_hit = alpha < alpha_floor
  )
}

family_feature_matrix <- function(membership, alpha) {
  A <- sweep(membership * 1, 2, alpha, "/")
  bad <- !is.finite(colSums(A))
  if (any(bad)) {
    A[, bad] <- 0
  }
  A
}

mixture_h <- function(q, membership, alpha) {
  A <- family_feature_matrix(membership, alpha)
  as.numeric(A %*% q)
}

mixture_distortions <- function(q,
                                membership,
                                alpha,
                                base_weights = NULL,
                                h_floor = 1e-12) {
  if (is.null(base_weights)) {
    base_weights <- rep(1 / nrow(membership), nrow(membership))
  }
  w <- normalize_weights(base_weights)
  q <- normalize_weights(q)
  h <- pmax(mixture_h(q, membership, alpha), h_floor)
  h_raw <- mixture_h(q, membership, alpha)
  mean_h <- sum(w * h_raw)
  kl_compressed_to_base <- sum(w * h_raw * ifelse(h_raw > 0, log(pmax(h_raw, h_floor)), 0))
  tv <- 0.5 * sum(w * abs(h_raw - 1))
  kl_base_to_compressed <- sum(w * (-log(h)))
  list(
    mean_h = mean_h,
    kl_compressed_to_base = kl_compressed_to_base,
    kl_base_to_compressed = kl_base_to_compressed,
    tv = tv,
    h = h_raw
  )
}

softmax_normalize <- function(log_q) {
  z <- log_q - max(log_q)
  q <- exp(z)
  q / sum(q)
}

enforce_simplex_safety <- function(q, q0_min = 0, q0_max = 1, safety_index = 1L) {
  q <- as.numeric(q)
  q[!is.finite(q) | q < 0] <- 0
  if (sum(q) <= 0) {
    q <- rep(1 / length(q), length(q))
  } else {
    q <- q / sum(q)
  }
  q0_min <- min(max(as.numeric(q0_min), 0), 1)
  q0_max <- min(max(as.numeric(q0_max), q0_min), 1)
  safety_index <- as.integer(safety_index)
  if (q0_min > 0 && length(q) > 1 && q[safety_index] < q0_min) {
    rest_idx <- setdiff(seq_along(q), safety_index)
    rest <- q[rest_idx]
    if (sum(rest) > 0) {
      rest <- (1 - q0_min) * rest / sum(rest)
    } else {
      rest <- rep((1 - q0_min) / length(rest_idx), length(rest_idx))
    }
    q[rest_idx] <- rest
    q[safety_index] <- q0_min
  }
  if (q0_max < 1 && length(q) > 1 && q[safety_index] > q0_max) {
    rest_idx <- setdiff(seq_along(q), safety_index)
    rest <- q[rest_idx]
    if (sum(rest) > 0) {
      rest <- (1 - q0_max) * rest / sum(rest)
    } else {
      rest <- rep((1 - q0_max) / length(rest_idx), length(rest_idx))
    }
    q[rest_idx] <- rest
    q[safety_index] <- q0_max
  }
  if (q0_min >= 1) {
    q[] <- 0
    q[safety_index] <- 1
  }
  q / sum(q)
}

simplex_directional_kkt_residual <- function(gradient,
                                             q,
                                             q0_min = 0,
                                             q0_max = 1,
                                             safety_index = 1L,
                                             active_tol = 1e-8) {
  gradient <- as.numeric(gradient)
  q <- normalize_weights(q)
  lower <- rep(0, length(q))
  lower[as.integer(safety_index)] <- min(max(as.numeric(q0_min), 0), 1)
  upper <- rep(1, length(q))
  upper[as.integer(safety_index)] <- min(max(as.numeric(q0_max), lower[as.integer(safety_index)]), 1)
  donors <- q > lower + active_tol
  receivers <- q < upper - active_tol
  if (!any(donors)) {
    return(0)
  }
  if (!any(receivers)) {
    return(0)
  }
  max(0, max(gradient[donors], na.rm = TRUE) - min(gradient[receivers], na.rm = TRUE))
}

evaluate_family_mixture_fit <- function(membership,
                                        alpha,
                                        costs,
                                        q,
                                        base_weights = NULL,
                                        beta = 0,
                                        tau = 1e-3,
                                        distortion = c("fkl", "tv", "rkl"),
                                        h_floor = 1e-10,
                                        q0_min = 0,
                                        q0_max = 1,
                                        safety_index = 1L,
                                        status = "feasible_evaluation") {
  distortion <- match.arg(distortion)
  if (is.null(base_weights)) {
    base_weights <- rep(1 / nrow(membership), nrow(membership))
  }
  w <- normalize_weights(base_weights)
  costs <- as.numeric(costs)
  q <- enforce_simplex_safety(q, q0_min = q0_min, q0_max = q0_max, safety_index = safety_index)
  h_floor_eff <- if (distortion == "fkl" && q0_min > 0) min(q0_min, 1) else h_floor
  A <- family_feature_matrix(membership, alpha)
  d <- mixture_distortions(q, membership, alpha, w, h_floor_eff)
  entropy_term <- if (tau > 0) tau * sum(ifelse(q > 0, q * log(q), 0)) else 0
  dist_value <- switch(
    distortion,
    fkl = d$kl_base_to_compressed,
    tv = d$tv,
    rkl = d$kl_compressed_to_base
  )
  h_final <- pmax(as.numeric(A %*% q), h_floor_eff)
  grad_dist <- switch(
    distortion,
    fkl = -as.numeric(crossprod(w / h_final, A)),
    tv = 0.5 * as.numeric(crossprod(w * sign(h_final - 1), A)),
    rkl = as.numeric(crossprod(w * (log(h_final) + 1), A))
  )
  grad_entropy <- if (tau > 0) tau * (log(pmax(q, h_floor)) + 1) else rep(0, length(q))
  grad <- grad_dist + beta * costs + grad_entropy
  entropy <- -sum(ifelse(q > 0, q * log(q), 0))
  list(
    q = q,
    objective = dist_value + beta * sum(q * costs) + entropy_term,
    distortions = d,
    expected_cost = sum(q * costs),
    entropy = entropy,
    effective_families = exp(entropy),
    n95 = which(cumsum(sort(q, decreasing = TRUE)) >= 0.95)[1],
    iterations = 0L,
    converged = TRUE,
    status = status,
    tol = NA_real_,
    kkt_tol = NA_real_,
    objective_change = 0,
    objective_trace = dist_value + beta * sum(q * costs) + entropy_term,
    line_search_failures = 0L,
    polished = FALSE,
    polish_convergence = NA_integer_,
    beta = beta,
    tau = tau,
    distortion = distortion,
    q0_min = q0_min,
    safety_index = safety_index,
    kkt_residual = simplex_directional_kkt_residual(
      gradient = grad,
      q = q,
      q0_min = q0_min,
      q0_max = q0_max,
      safety_index = safety_index
    )
  )
}

optimize_family_mixture <- function(membership,
                                    alpha,
                                    costs,
                                    base_weights = NULL,
                                    beta = 0,
                                    tau = 1e-3,
                                    distortion = c("fkl", "tv", "rkl"),
                                    max_iter = 80L,
                                    step = 0.1,
                                    h_floor = 1e-10,
                                    tol = 1e-6,
                                    kkt_tol = 1e-4,
                                    polish = FALSE,
                                    polish_maxit = 200L,
                                    q_init = NULL,
                                    q0_min = 0,
                                    q0_max = 1,
                                    safety_index = 1L) {
  distortion <- match.arg(distortion)
  if (is.null(base_weights)) {
    base_weights <- rep(1 / nrow(membership), nrow(membership))
  }
  w <- normalize_weights(base_weights)
  costs <- as.numeric(costs)
  A <- family_feature_matrix(membership, alpha)
  M <- ncol(A)
  if (is.null(q_init)) {
    q <- rep(1 / M, M)
  } else {
    q <- normalize_weights(q_init)
  }
  q <- enforce_simplex_safety(q, q0_min = q0_min, q0_max = q0_max, safety_index = safety_index)
  h_floor_eff <- if (distortion == "fkl" && q0_min > 0) min(q0_min, 1) else h_floor

  objective_value <- function(q) {
    q <- enforce_simplex_safety(q, q0_min = q0_min, q0_max = q0_max, safety_index = safety_index)
    d <- mixture_distortions(q, membership, alpha, w, h_floor_eff)
    entropy_term <- if (tau > 0) tau * sum(ifelse(q > 0, q * log(q), 0)) else 0
    dist_value <- switch(
      distortion,
      fkl = d$kl_base_to_compressed,
      tv = d$tv,
      rkl = d$kl_compressed_to_base
    )
    dist_value + beta * sum(q * costs) + entropy_term
  }

  gradient_value <- function(q) {
    q <- enforce_simplex_safety(q, q0_min = q0_min, q0_max = q0_max, safety_index = safety_index)
    h <- pmax(as.numeric(A %*% q), h_floor_eff)
    grad_dist <- switch(
      distortion,
      fkl = -as.numeric(crossprod(w / h, A)),
      tv = 0.5 * as.numeric(crossprod(w * sign(h - 1), A)),
      rkl = as.numeric(crossprod(w * (log(h) + 1), A))
    )
    grad_entropy <- if (tau > 0) tau * (log(pmax(q, h_floor)) + 1) else rep(0, M)
    grad_dist + beta * costs + grad_entropy
  }

  obj_old <- objective_value(q)
  obj_trace <- obj_old
  converged <- FALSE
  objective_change <- Inf
  line_search_failures <- 0L
  for (iter in seq_len(max_iter)) {
    grad <- gradient_value(q)
    q_new <- softmax_normalize(log(pmax(q, h_floor)) - step * (grad - mean(grad)))
    q_new <- enforce_simplex_safety(q_new, q0_min = q0_min, q0_max = q0_max, safety_index = safety_index)
    obj_new <- objective_value(q_new)
    local_step <- step
    while (obj_new > obj_old && local_step > 1e-6) {
      local_step <- local_step / 2
      q_new <- softmax_normalize(log(pmax(q, h_floor)) - local_step * (grad - mean(grad)))
      q_new <- enforce_simplex_safety(q_new, q0_min = q0_min, q0_max = q0_max, safety_index = safety_index)
      obj_new <- objective_value(q_new)
    }
    if (obj_new > obj_old) {
      line_search_failures <- line_search_failures + 1L
    }
    objective_change <- abs(obj_old - obj_new)
    obj_trace <- c(obj_trace, obj_new)
    if (objective_change <= tol * (1 + abs(obj_old))) {
      q <- q_new
      obj_old <- obj_new
      converged <- TRUE
      break
    }
    q <- q_new
    obj_old <- obj_new
    step <- min(0.5, local_step * 1.05)
  }

  polished <- FALSE
  polish_convergence <- NA_integer_
  if (isTRUE(polish) && M > 1L && q0_min < 1) {
    lower <- rep(0, M)
    lower[safety_index] <- min(max(as.numeric(q0_min), 0), 1)
    scale <- 1 - sum(lower)
    q_to_theta <- function(q) {
      r <- pmax((q - lower) / scale, h_floor)
      r <- r / sum(r)
      log(r[-M] / r[M])
    }
    theta_to_q <- function(theta) {
      z <- c(theta, 0)
      z <- z - max(z)
      r <- exp(z)
      r <- r / sum(r)
      lower + scale * r
    }
    theta0 <- q_to_theta(q)
    obj_theta <- function(theta) objective_value(theta_to_q(theta))
    grad_theta <- function(theta) {
      z <- c(theta, 0)
      z <- z - max(z)
      r <- exp(z)
      r <- r / sum(r)
      q_theta <- lower + scale * r
      grad_q <- gradient_value(q_theta)
      scale * r[-M] * (grad_q[-M] - sum(r * grad_q))
    }
    opt <- tryCatch(
      stats::optim(
        par = theta0,
        fn = obj_theta,
        gr = grad_theta,
        method = "BFGS",
        control = list(maxit = as.integer(polish_maxit), reltol = tol)
      ),
      error = function(e) e
    )
    if (!inherits(opt, "error")) {
      q_polished <- enforce_simplex_safety(theta_to_q(opt$par), q0_min = q0_min, q0_max = q0_max, safety_index = safety_index)
      obj_polished <- objective_value(q_polished)
      polish_convergence <- opt$convergence
      if (is.finite(obj_polished) && obj_polished <= obj_old + tol * (1 + abs(obj_old))) {
        objective_change <- abs(obj_old - obj_polished)
        q <- q_polished
        obj_old <- obj_polished
        obj_trace <- c(obj_trace, obj_polished)
        polished <- TRUE
        converged <- converged || identical(opt$convergence, 0L)
      }
    }
  }

  d <- mixture_distortions(q, membership, alpha, w, h_floor_eff)
  h_final <- pmax(as.numeric(A %*% q), h_floor_eff)
  grad_dist_final <- switch(
    distortion,
    fkl = -as.numeric(crossprod(w / h_final, A)),
    tv = 0.5 * as.numeric(crossprod(w * sign(h_final - 1), A)),
    rkl = as.numeric(crossprod(w * (log(h_final) + 1), A))
  )
  grad_entropy_final <- if (tau > 0) tau * (log(pmax(q, h_floor)) + 1) else rep(0, M)
  grad_final <- grad_dist_final + beta * costs + grad_entropy_final
  kkt_residual <- simplex_directional_kkt_residual(
    gradient = grad_final,
    q = q,
    q0_min = q0_min,
    q0_max = q0_max,
    safety_index = safety_index
  )
  status <- if (converged && kkt_residual <= kkt_tol) {
    "converged"
  } else if (converged) {
    "objective_converged_kkt_warning"
  } else {
    "max_iter_reached"
  }
  list(
    q = q,
    objective = obj_old,
    distortions = d,
    expected_cost = sum(q * costs),
    entropy = -sum(ifelse(q > 0, q * log(q), 0)),
    effective_families = exp(-sum(ifelse(q > 0, q * log(q), 0))),
    n95 = which(cumsum(sort(q, decreasing = TRUE)) >= 0.95)[1],
    iterations = iter,
    converged = converged,
    status = status,
    tol = tol,
    kkt_tol = kkt_tol,
    objective_change = objective_change,
    objective_trace = obj_trace,
    line_search_failures = line_search_failures,
    polished = polished,
    polish_convergence = polish_convergence,
    beta = beta,
    tau = tau,
    distortion = distortion,
    q0_min = q0_min,
    q0_max = q0_max,
    safety_index = safety_index,
    kkt_residual = kkt_residual
  )
}

rate_distortion_grid <- function(membership,
                                 alpha,
                                 costs,
                                 base_weights = NULL,
                                 beta_grid = exp(seq(log(0.001), log(10), length.out = 16)),
                                 tau = 1e-3,
                                 distortion = "fkl",
                                 safety_index = 1L,
                                 q0_min = 0,
                                 q0_max = 1) {
  rows <- list()
  q_start <- NULL
  for (b in beta_grid) {
    fit <- optimize_family_mixture(
      membership = membership,
      alpha = alpha,
      costs = costs,
      base_weights = base_weights,
      beta = b,
      tau = tau,
      distortion = distortion,
      q_init = q_start,
      q0_min = q0_min,
      q0_max = q0_max,
      safety_index = safety_index
    )
    q_start <- fit$q
    rows[[length(rows) + 1L]] <- data.frame(
      beta = b,
      tau = tau,
      distortion = distortion,
      tv = fit$distortions$tv,
      kl_base_to_compressed = fit$distortions$kl_base_to_compressed,
      kl_compressed_to_base = fit$distortions$kl_compressed_to_base,
      mean_h = fit$distortions$mean_h,
      expected_cost = fit$expected_cost,
      effective_families = fit$effective_families,
      family_n95 = fit$n95,
      q0 = fit$q[safety_index],
      objective = fit$objective,
      iterations = fit$iterations,
      converged = fit$converged,
      optimizer_status = fit$status,
      objective_change = fit$objective_change,
      kkt_residual = fit$kkt_residual,
      q0_min = fit$q0_min,
      q0_max = fit$q0_max,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

verify_mixture_identities_exact <- function(base_weights,
                                            membership,
                                            alpha,
                                            q,
                                            h_floor = 1e-12) {
  w <- normalize_weights(base_weights)
  h <- mixture_h(q, membership, alpha)
  pi_bar <- w * h
  pi_bar <- pi_bar / sum(pi_bar)
  d <- mixture_distortions(q, membership, alpha, w, h_floor)
  direct_rkl <- sum(ifelse(pi_bar > 0, pi_bar * log(pi_bar / w), 0))
  direct_tv <- 0.5 * sum(abs(pi_bar - w))
  direct_fkl <- if (all(pi_bar > 0)) sum(w * log(w / pi_bar)) else Inf
  data.frame(
    mean_h_error = abs(d$mean_h - 1),
    rkl_identity_error = abs(d$kl_compressed_to_base - direct_rkl),
    tv_identity_error = abs(d$tv - direct_tv),
    fkl_identity_error = abs(d$kl_base_to_compressed - direct_fkl),
    stringsAsFactors = FALSE
  )
}

posterior_functional_error <- function(phi,
                                       q,
                                       membership,
                                       alpha,
                                       base_weights,
                                       h_floor = 1e-12) {
  phi <- as.numeric(phi)
  membership <- as.matrix(membership)
  if (length(phi) != nrow(membership)) {
    stop("phi must have one value per posterior support row")
  }
  weights <- normalize_weights(base_weights)
  q <- normalize_weights(q)
  h <- mixture_h(q, membership, alpha)
  compressed_weights <- weights * pmax(h, h_floor)
  compressed_weights <- compressed_weights / sum(compressed_weights)
  d <- mixture_distortions(q, membership, alpha, weights, h_floor = h_floor)
  b <- max(abs(phi), na.rm = TRUE)
  data.frame(
    reference_mean = sum(weights * phi),
    compressed_mean = sum(compressed_weights * phi),
    absolute_error = abs(sum(compressed_weights * phi) - sum(weights * phi)),
    tv = d$tv,
    tv_bound = 2 * b * d$tv,
    bound_holds = abs(sum(compressed_weights * phi) - sum(weights * phi)) <= 2 * b * d$tv + 1e-10,
    stringsAsFactors = FALSE
  )
}

blocked_mcse <- function(x, block_length) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  block_length <- as.integer(block_length)
  if (block_length <= 0) stop("block_length must be positive")
  n_blocks <- floor(length(x) / block_length)
  if (n_blocks < 2L) {
    return(data.frame(
      n = length(x),
      block_length = block_length,
      n_blocks = n_blocks,
      estimate = mean(x),
      mcse = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  x <- x[seq_len(n_blocks * block_length)]
  block_id <- rep(seq_len(n_blocks), each = block_length)
  block_means <- as.numeric(tapply(x, block_id, mean))
  data.frame(
    n = length(x),
    block_length = block_length,
    n_blocks = n_blocks,
    estimate = mean(block_means),
    mcse = stats::sd(block_means) / sqrt(n_blocks),
    stringsAsFactors = FALSE
  )
}

blocked_validation_summary <- function(values, block_lengths = c(5L, 10L, 20L, 50L)) {
  do.call(rbind, lapply(block_lengths, function(b) blocked_mcse(values, b)))
}
