source(file.path("sim", "src", "support_kernel_benchmark.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit)) sub(pat, "", hit[[1]]) else default
}

mode <- tolower(get_arg("mode", "smoke"))
if (!mode %in% c("smoke", "full")) stop("--mode must be smoke or full")

table_dir <- file.path("sim", "output", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

config <- list(
  rep_count = if (mode == "full") 3L else 1L,
  beta_grid = if (mode == "full") c(0.02, 0.08) else 0.02,
  topm_grid = if (mode == "full") c(16L, 32L) else 16L,
  cluster_grid = if (mode == "full") c(6L, 10L) else 6L,
  credible_grid = if (mode == "full") c(0.90, 0.95) else 0.90
)

stable_loglik_logistic <- function(eta, y) {
  sum(y * eta - ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta))))
}

logistic_laplace_marginal <- function(gamma, X, y, tau2 = 4, intercept_var = 100) {
  active <- which(gamma == 1)
  Z <- if (length(active)) cbind(1, X[, active, drop = FALSE]) else matrix(1, nrow(X), 1)
  d <- ncol(Z)
  prior_var <- c(intercept_var, rep(tau2, d - 1L))
  objective <- function(beta) {
    eta <- as.numeric(Z %*% beta)
    loglik <- stable_loglik_logistic(eta, y)
    logprior <- -0.5 * sum(beta^2 / prior_var) - 0.5 * sum(log(2 * pi * prior_var))
    -(loglik + logprior)
  }
  opt <- tryCatch(stats::optim(rep(0, d), objective, method = "BFGS", control = list(maxit = 120)), error = function(e) e)
  if (inherits(opt, "error") || !all(is.finite(opt$par))) return(-Inf)
  beta <- opt$par
  eta <- as.numeric(Z %*% beta)
  prob <- pmin(pmax(stats::plogis(eta), 1e-8), 1 - 1e-8)
  W <- prob * (1 - prob)
  H <- crossprod(Z * sqrt(W)) + diag(1 / prior_var, d)
  chol_H <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(chol_H)) return(-Inf)
  log_integrand <- -objective(beta)
  as.numeric(log_integrand + 0.5 * d * log(2 * pi) - sum(log(diag(chol_H))))
}

fit_logistic_reference <- function(X, y, theta = 0.25, tau2 = 4) {
  supports <- enumerate_supports(ncol(X))
  log_marginal <- numeric(nrow(supports))
  for (i in seq_len(nrow(supports))) {
    log_marginal[i] <- logistic_laplace_marginal(supports[i, ], X = X, y = y, tau2 = tau2)
  }
  log_prior <- log_prior_bernoulli(supports, theta)
  log_post <- log_marginal + log_prior
  log_post <- log_post - log_sum_exp(log_post)
  list(
    supports = supports,
    log_marginal = log_marginal,
    log_prior = log_prior,
    posterior = exp(log_post),
    log_posterior = log_post,
    theta = theta,
    tau2 = tau2,
    a0 = NA_real_,
    b0 = NA_real_,
    n = nrow(X),
    p = ncol(X)
  )
}

simulate_logistic_cell <- function(rep_id) {
  set.seed(20260531L + rep_id)
  K <- 4L
  m <- 3L
  p <- K * m
  n <- 110L
  group_id <- make_group_id(K, m)
  X_raw <- simulate_block_design(n, K, m, rho = 0.90)
  X <- scale(X_raw)
  X <- as.matrix(X)
  beta <- rep(0, p)
  beta[c(1, 7, 11)] <- c(1.05, -0.95, 0.80)
  eta <- -0.20 + as.numeric(X %*% beta)
  prob <- stats::plogis(eta)
  y <- stats::rbinom(n, 1, prob)
  list(X_train = X, y_train = y, group_id = group_id, K = K, m = m, scenario = "logistic_redundancy")
}

kernel_row <- function(method, fit, W, alpha, q, costs, beta, tuning, family, stored_atoms = NA_real_) {
  q <- normalize_weights(q)
  d <- mixture_distortions(q, W, alpha, fit$posterior, h_floor = max(q[1], 1e-10))
  entropy <- -sum(ifelse(q > 0, q * log(q), 0))
  data.frame(
    method = method,
    method_family = family,
    tuning = tuning,
    beta = beta,
    tv = d$tv,
    fkl = d$kl_base_to_compressed,
    rkl = d$kl_compressed_to_base,
    q0 = q[1],
    expected_code = sum(q * costs),
    active_kernels_001 = sum(q > 0.001),
    q_effective_kernels = exp(entropy),
    n_kernels = ncol(W),
    stored_atoms = stored_atoms,
    stringsAsFactors = FALSE
  )
}

run_cell <- function(rep_id) {
  dat <- simulate_logistic_cell(rep_id)
  ref <- fit_logistic_reference(dat$X_train, dat$y_train, theta = 3 / ncol(dat$X_train), tau2 = 4)
  rows <- list()
  add <- function(x) rows[[length(rows) + 1L]] <<- x
  for (beta in config$beta_grid) {
    ask <- fit_askpc_pooled_pruned(ref$supports, ref$posterior, group_id = dat$group_id, X = dat$X_train, beta = beta, mode = "smoke", max_iter = 90L)
    add(kernel_row("ASK-PC pooled-pruned 99%", ref, ask$pruned$W, ask$pruned$alpha, ask$pruned$fit$q, ask$pruned$costs, beta, paste0("beta=", beta), "support-kernel compression"))
    fixed <- fit_fixed_hard_benchmark(dat, ref, beta = beta, capacity = 1L)
    add(kernel_row("Fixed hard dictionary", ref, fixed$W, fixed$alpha, fixed$q, fixed$costs, beta, paste0("beta=", beta), "fixed region dictionary"))
    for (k in config$cluster_grid) {
      cl <- fit_posterior_clustering_benchmark(dat, ref, beta = beta, n_clusters = k, seed = 99L + rep_id + k)
      add(kernel_row("Posterior clustering", ref, cl$W, cl$alpha, cl$q, cl$costs, beta, paste0("clusters=", k, ", beta=", beta), "posterior clustering"))
    }
    for (top_m in config$topm_grid) {
      top <- fit_topm_benchmark(ref, beta = beta, top_m = top_m)
      add(kernel_row("Top-M support atoms", ref, top$W, top$alpha, top$q, top$costs, beta, paste0("M=", top_m, ", beta=", beta), "atom truncation", stored_atoms = top_m))
    }
    for (coverage in config$credible_grid) {
      cred <- fit_credible_support_benchmark(ref, beta = beta, coverage = coverage)
      add(kernel_row("Credible support set", ref, cred$W, cred$alpha, cred$q, cred$costs, beta, paste0("coverage=", coverage, ", beta=", beta), "credible support summary", stored_atoms = cred$stored_atoms))
    }
  }
  out <- do.call(rbind, rows)
  out$replication_id <- rep_id
  out$mode <- mode
  out$reference <- "logistic Laplace support posterior"
  out$reference_support_n95 <- posterior_summary(ref$posterior)$n_mass
  out
}

se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

detail <- do.call(rbind, lapply(seq_len(config$rep_count), run_cell))
write.csv(detail, file.path(table_dir, "table_logistic_support_kernel_check_detail.csv"), row.names = FALSE)

best_parts <- split(detail, interaction(detail$replication_id, detail$method, drop = TRUE), drop = TRUE)
best <- do.call(rbind, lapply(best_parts, function(d) d[which.min(d$fkl), , drop = FALSE]))
summary <- do.call(rbind, lapply(split(best, best$method, drop = TRUE), function(d) {
  data.frame(
    method = d$method[1],
    tv_mean = mean(d$tv),
    tv_se = se(d$tv),
    fkl_mean = mean(d$fkl),
    fkl_se = se(d$fkl),
    q0_mean = mean(d$q0),
    q0_se = se(d$q0),
    expected_code_mean = mean(d$expected_code),
    expected_code_se = se(d$expected_code),
    active_kernels_mean = mean(d$active_kernels_001),
    active_kernels_se = se(d$active_kernels_001),
    n_replications = length(unique(d$replication_id)),
    mode = mode,
    stringsAsFactors = FALSE
  )
}))
method_order <- c("ASK-PC pooled-pruned 99%", "Posterior clustering", "Top-M support atoms", "Credible support set", "Fixed hard dictionary")
summary <- summary[order(match(summary$method, method_order)), ]
write.csv(summary, file.path(table_dir, "table_logistic_support_kernel_check.csv"), row.names = FALSE)

fmt <- function(x, digits = 3) {
  if (!is.finite(x)) return("--")
  formatC(x, format = "f", digits = digits)
}
fmt_mean <- function(mu, se_val, digits = 3) {
  if (!is.finite(mu)) return("--")
  if (!is.finite(se_val)) return(fmt(mu, digits))
  paste0(fmt(mu, digits), " (", fmt(se_val, digits), ")")
}
labels <- c(
  "ASK-PC pooled-pruned 99%" = "Pooled-pruned",
  "Posterior clustering" = "Cluster kernels",
  "Top-M support atoms" = "Top-M atoms",
  "Credible support set" = "Credible set",
  "Fixed hard dictionary" = "Fixed regions"
)
tex <- c(
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "Method & TV & FKL & $q_0$ & Storage & Active \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(summary))) {
  tex <- c(
    tex,
    sprintf(
      "%s & %s & %s & %s & %s & %s \\\\",
      labels[summary$method[i]],
      fmt_mean(summary$tv_mean[i], summary$tv_se[i]),
      fmt_mean(summary$fkl_mean[i], summary$fkl_se[i]),
      fmt_mean(summary$q0_mean[i], summary$q0_se[i]),
      fmt_mean(summary$expected_code_mean[i], summary$expected_code_se[i], 2),
      fmt_mean(summary$active_kernels_mean[i], summary$active_kernels_se[i], 1)
    )
  )
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(table_dir, "table_logistic_support_kernel_check.tex"))

cat("Wrote logistic support-kernel check outputs\n")
