source(file.path("sim", "src", "support_kernel_benchmark.R"))
source(file.path("sim", "mcmc_bma.R"))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit)) sub(pat, "", hit[[1]]) else default
}

has_flag <- function(name) {
  paste0("--", name) %in% args
}

mode_arg <- get_arg("mode", "smoke")
if (!mode_arg %in% c("smoke", "pilot", "medium", "full", "summarize")) {
  stop("--mode must be smoke, pilot, medium, full, or summarize")
}

run_tag <- get_arg("run-tag", "")
if (nzchar(run_tag) && !grepl("^[A-Za-z0-9_.-]+$", run_tag)) {
  stop("--run-tag may only contain letters, numbers, dots, underscores, and hyphens")
}

tag_file <- function(filename, tag = run_tag) {
  if (!nzchar(tag)) {
    return(filename)
  }
  sub("(\\.[^.]+)$", paste0("_", tag, "\\1"), filename)
}

table_dir <- file.path("sim", "output", "tables")
object_dir <- file.path("sim", "output", "large_end_to_end")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

large_mode_config <- function(mode) {
  stopifnot(mode %in% c("smoke", "pilot", "medium", "full"))
  list(
    mode = mode,
    p = switch(mode, smoke = 60L, pilot = 100L, medium = 100L, full = 100L),
    n_train = switch(mode, smoke = 80L, pilot = 100L, medium = 110L, full = 120L),
    n_test = switch(mode, smoke = 160L, pilot = 240L, medium = 300L, full = 400L),
    group_size = 5L,
    rep_count = switch(mode, smoke = 1L, pilot = 1L, medium = 2L, full = 3L),
    rho_grid = switch(mode, smoke = 0.85, pilot = 0.85, medium = c(0.7, 0.9), full = c(0.7, 0.9)),
    scenario_grid = switch(
      mode,
      smoke = c("one_representative"),
      pilot = c("one_representative", "ordered_interval"),
      medium = c("one_representative", "multi_representative", "weak_signal", "noisy_grouping"),
      full = c(
        "one_representative",
        "multi_representative",
        "weak_signal",
        "noisy_grouping",
        "ordered_interval",
        "graph_community"
      )
    ),
    beta_grid = switch(mode, smoke = 0.02, pilot = c(0.01, 0.04), medium = c(0.005, 0.02, 0.08), full = c(0.005, 0.02, 0.08)),
    topm_grid = switch(mode, smoke = 24L, pilot = c(24L, 48L), medium = c(24L, 48L, 96L), full = c(24L, 48L, 96L)),
    credible_grid = switch(mode, smoke = 0.90, pilot = c(0.90, 0.95), medium = c(0.80, 0.90, 0.95), full = c(0.80, 0.90, 0.95)),
    cluster_grid = switch(mode, smoke = 8L, pilot = c(8L, 12L), medium = c(8L, 12L, 16L), full = c(8L, 12L, 16L)),
    dilution_power_grid = switch(mode, smoke = 0.5, pilot = c(0.5), medium = c(0.25, 0.5, 1), full = c(0.25, 0.5, 1)),
    dpp_power_grid = switch(mode, smoke = 1, pilot = 1, medium = c(0.5, 1), full = c(0.5, 1)),
    mcmc = list(
      n_iter = switch(mode, smoke = 700L, pilot = 2200L, medium = 8000L, full = 20000L),
      burn = switch(mode, smoke = 250L, pilot = 800L, medium = 2500L, full = 5000L),
      thin = switch(mode, smoke = 5L, pilot = 5L, medium = 5L, full = 5L),
      n_chains = switch(mode, smoke = 2L, pilot = 3L, medium = 4L, full = 6L),
      max_size = switch(mode, smoke = 10L, pilot = 12L, medium = 14L, full = 16L),
      start_size = switch(mode, smoke = 5L, pilot = 6L, medium = 7L, full = 8L)
    )
  )
}

override_config <- function(config) {
  scalar_int <- function(name, current) {
    val <- get_arg(name, NULL)
    if (is.null(val)) current else as.integer(val)
  }
  scalar_num_vec <- function(name, current) {
    val <- get_arg(name, NULL)
    if (is.null(val)) current else as.numeric(strsplit(val, ",", fixed = TRUE)[[1]])
  }
  scalar_chr_vec <- function(name, current) {
    val <- get_arg(name, NULL)
    if (is.null(val)) current else strsplit(val, ",", fixed = TRUE)[[1]]
  }
  config$p <- scalar_int("p", config$p)
  config$n_train <- scalar_int("n-train", config$n_train)
  config$n_test <- scalar_int("n-test", config$n_test)
  config$rep_count <- scalar_int("replications", config$rep_count)
  config$rho_grid <- scalar_num_vec("rho", config$rho_grid)
  config$scenario_grid <- scalar_chr_vec("scenarios", config$scenario_grid)
  config$beta_grid <- scalar_num_vec("beta", config$beta_grid)
  config$mcmc$n_iter <- scalar_int("n-iter", config$mcmc$n_iter)
  config$mcmc$burn <- scalar_int("burn", config$mcmc$burn)
  config$mcmc$thin <- scalar_int("thin", config$mcmc$thin)
  config$mcmc$n_chains <- scalar_int("chains", config$mcmc$n_chains)
  config$mcmc$max_size <- scalar_int("max-size", config$mcmc$max_size)
  config$mcmc$start_size <- scalar_int("start-size", config$mcmc$start_size)
  config$run_tag <- run_tag
  config
}

simulate_large_askpc_data <- function(scenario, rho, seed, config) {
  set.seed(seed)
  p <- as.integer(config$p)
  m <- as.integer(config$group_size)
  K <- ceiling(p / m)
  p <- K * m
  group_id <- make_group_id(K, m)
  n_train <- as.integer(config$n_train)
  n_test <- as.integer(config$n_test)

  X_train_raw <- simulate_block_design(n_train, K, m, rho)
  X_test_raw <- simulate_block_design(n_test, K, m, rho)
  beta <- rep(0, p)
  sigma <- 1.25

  if (identical(scenario, "ordered_interval")) {
    z_train <- as.numeric(stats::filter(rnorm(n_train + p), filter = rep(1, 5) / 5, sides = 1))
    z_test <- as.numeric(stats::filter(rnorm(n_test + p), filter = rep(1, 5) / 5, sides = 1))
    z_train[!is.finite(z_train)] <- rnorm(sum(!is.finite(z_train)))
    z_test[!is.finite(z_test)] <- rnorm(sum(!is.finite(z_test)))
    X_train_raw <- matrix(0, n_train, p)
    X_test_raw <- matrix(0, n_test, p)
    for (j in seq_len(p)) {
      X_train_raw[, j] <- sqrt(rho) * z_train[j:(j + n_train - 1L)] + sqrt(1 - rho) * rnorm(n_train)
      X_test_raw[, j] <- sqrt(rho) * z_test[j:(j + n_test - 1L)] + sqrt(1 - rho) * rnorm(n_test)
    }
    center <- floor(p / 2)
    band <- (center - 2L):(center + 2L)
    beta[band] <- c(0.6, 0.85, 1.05, 0.85, 0.6)
  } else {
    active_groups <- c(2L, 5L, 9L, 13L, 17L)
    active_groups <- active_groups[active_groups <= K]
    rep1 <- (active_groups - 1L) * m + 1L
    if (identical(scenario, "multi_representative")) {
      beta[rep1] <- c(1.15, 1.05, 0.9, 0.8, 0.7)[seq_along(rep1)]
      second <- pmin(rep1 + 1L, p)
      beta[second] <- beta[second] + c(0.55, 0.5, 0.45, 0.4, 0.35)[seq_along(second)]
    } else if (identical(scenario, "weak_signal")) {
      beta[rep1] <- c(0.75, 0.68, 0.6, 0.55, 0.5)[seq_along(rep1)]
      sigma <- 1.7
    } else if (identical(scenario, "graph_community")) {
      beta[rep1] <- c(0.9, 0.85, 0.75, 0.65, 0.6)[seq_along(rep1)]
      beta[pmin(rep1 + 1L, p)] <- beta[pmin(rep1 + 1L, p)] + 0.35
    } else {
      beta[rep1] <- c(1.25, 1.1, 0.95, 0.85, 0.75)[seq_along(rep1)]
    }
  }

  fit_group_id <- group_id
  if (identical(scenario, "noisy_grouping")) {
    swap <- sample(seq_along(fit_group_id), size = max(1L, floor(0.20 * length(fit_group_id))))
    fit_group_id[swap] <- sample(fit_group_id[swap])
  }

  z <- standardize_with_training(X_train_raw, X_test_raw)
  y_train_raw <- as.numeric(z$X_train %*% beta + rnorm(n_train, sd = sigma))
  y_test_raw <- as.numeric(z$X_test %*% beta + rnorm(n_test, sd = sigma))
  y_center <- mean(y_train_raw)
  list(
    X_train = z$X_train,
    X_test = z$X_test,
    y_train = y_train_raw - y_center,
    y_test = y_test_raw - y_center,
    beta = beta,
    true_support = as.integer(beta != 0),
    group_id = as.integer(fit_group_id),
    true_group_id = as.integer(group_id),
    active_groups = sort(unique(group_id[beta != 0])),
    K = max(fit_group_id),
    m = m,
    rho = rho,
    sigma = sigma,
    scenario = scenario
  )
}

enrich_mcmc_reference_fit <- function(fit, X, y) {
  supports <- as.matrix(fit$supports)
  log_marginal <- numeric(nrow(supports))
  for (i in seq_len(nrow(supports))) {
    log_marginal[i] <- log_marginal_support(
      gamma = supports[i, ],
      X = X,
      y = y,
      tau2 = fit$tau2,
      a0 = fit$a0,
      b0 = fit$b0
    )
  }
  log_prior <- log_prior_bernoulli(supports, fit$theta)
  fit$log_marginal <- log_marginal
  fit$log_prior <- log_prior
  fit$log_posterior <- log(pmax(fit$posterior, 1e-300))
  fit$log_evidence <- NA_real_
  fit
}

ess_one <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4L || stats::var(x) == 0) return(n)
  ac <- as.numeric(stats::acf(x, plot = FALSE, lag.max = min(100L, n - 1L))$acf)[-1L]
  first_negative <- which(ac < 0)
  if (length(first_negative)) ac <- ac[seq_len(max(1L, first_negative[[1]] - 1L))]
  tau <- 1 + 2 * sum(pmax(ac, 0), na.rm = TRUE)
  max(1, min(n, n / tau))
}

split_rhat_one <- function(chains) {
  parts <- unlist(lapply(chains, function(x) {
    n <- length(x)
    h <- floor(n / 2L)
    if (h < 2L) return(list())
    list(as.numeric(x[seq_len(h)]), as.numeric(x[(h + 1L):(2L * h)]))
  }), recursive = FALSE)
  if (length(parts) < 2L) return(NA_real_)
  n <- min(vapply(parts, length, integer(1)))
  mat <- do.call(rbind, lapply(parts, function(x) x[seq_len(n)]))
  W <- mean(apply(mat, 1, stats::var), na.rm = TRUE)
  B <- n * stats::var(rowMeans(mat), na.rm = TRUE)
  if (!is.finite(W) || W <= 0) return(NA_real_)
  sqrt((((n - 1) / n) * W + B / n) / W)
}

mcmc_reference_diagnostics <- function(fit, group_id, runtime_sec) {
  chains <- fit$samples_by_chain
  group_chains <- lapply(chains, function(s) support_group_counts(s, group_id) > 0)
  all_group <- do.call(rbind, group_chains)
  pip <- colMeans(all_group)
  top_groups <- order(pip, decreasing = TRUE)[seq_len(min(20L, length(pip)))]
  size_chains <- lapply(chains, rowSums)
  logp_chains <- split(fit$traces$log_posterior, fit$traces$chain)
  rhat_vals <- c(
    split_rhat_one(size_chains),
    split_rhat_one(logp_chains),
    vapply(top_groups, function(g) split_rhat_one(lapply(group_chains, function(z) z[, g] * 1)), numeric(1))
  )
  ess_vals <- c(
    ess_one(unlist(size_chains)),
    ess_one(unlist(logp_chains)),
    vapply(top_groups, function(g) ess_one(unlist(lapply(group_chains, function(z) z[, g] * 1))), numeric(1))
  )
  pip_mcse <- vapply(top_groups, function(g) {
    z <- unlist(lapply(group_chains, function(x) x[, g] * 1))
    sqrt(stats::var(z) / max(1, ess_one(z)))
  }, numeric(1))
  diagnostic_warning <- character(0)
  if (is.finite(fit$accept_rate) && fit$accept_rate < 0.03) {
    diagnostic_warning <- c(diagnostic_warning, "low acceptance")
  }
  status <- "pass"
  if (!is.finite(max(rhat_vals, na.rm = TRUE)) ||
      max(rhat_vals, na.rm = TRUE) > 1.25 ||
      min(ess_vals, na.rm = TRUE) < 50 ||
      fit$group_pip_max_range > 0.45) {
    status <- "stress-test"
  }
  data.frame(
    reference_mode = "mcmc_empirical_visited_supports",
    sampler = if (!is.null(fit$sampler_label)) fit$sampler_label else "multi_chain_mc3",
    n_iter = if (!is.null(fit$n_iter)) fit$n_iter else NA_integer_,
    burn = if (!is.null(fit$burn)) fit$burn else NA_integer_,
    thin = if (!is.null(fit$thin)) fit$thin else NA_integer_,
    n_chains = fit$n_chains,
    retained_draws = nrow(fit$samples),
    unique_supports = fit$unique_count,
    unique_group_states = length(unique(apply(all_group, 1, paste0, collapse = ""))),
    top_support_mass = max(fit$posterior),
    effective_support_count = posterior_summary(fit$posterior)$effective_count,
    mean_size = fit$mean_size,
    accept_rate = fit$accept_rate,
    group_pip_max_range = fit$group_pip_max_range,
    group_pip_mean_range = fit$group_pip_mean_range,
    split_rhat_max = suppressWarnings(max(rhat_vals, na.rm = TRUE)),
    split_rhat_median = suppressWarnings(stats::median(rhat_vals, na.rm = TRUE)),
    ess_min = suppressWarnings(min(ess_vals, na.rm = TRUE)),
    ess_median = suppressWarnings(stats::median(ess_vals, na.rm = TRUE)),
    group_pip_mcse_max = suppressWarnings(max(pip_mcse, na.rm = TRUE)),
    runtime_sec_reference = runtime_sec,
    diagnostic_status = status,
    diagnostic_warning = if (length(diagnostic_warning)) paste(diagnostic_warning, collapse = "; ") else "",
    stringsAsFactors = FALSE
  )
}

normalize_reference_diagnostics <- function(detail) {
  needed <- c("split_rhat_max", "ess_min", "group_pip_max_range")
  if (!all(needed %in% names(detail))) {
    return(detail)
  }
  num <- function(x) suppressWarnings(as.numeric(x))
  core_fail <- !is.finite(num(detail$split_rhat_max)) |
    num(detail$split_rhat_max) > 1.25 |
    !is.finite(num(detail$ess_min)) |
    num(detail$ess_min) < 50 |
    !is.finite(num(detail$group_pip_max_range)) |
    num(detail$group_pip_max_range) > 0.45
  detail$diagnostic_status <- ifelse(core_fail, "stress-test", "pass")
  if ("accept_rate" %in% names(detail)) {
    warning <- ifelse(is.finite(num(detail$accept_rate)) & num(detail$accept_rate) < 0.03,
                      "low acceptance", "")
    if ("diagnostic_warning" %in% names(detail)) {
      existing <- as.character(detail$diagnostic_warning)
      existing[is.na(existing)] <- ""
      detail$diagnostic_warning <- ifelse(nzchar(existing), existing, warning)
    } else {
      detail$diagnostic_warning <- warning
    }
  }
  detail
}

cell_key <- function(scenario, rho, replication_id) {
  paste(scenario, sprintf("%.6f", rho), replication_id, sep = "|")
}

benchmark_grid <- function(config) {
  expand.grid(
    scenario = config$scenario_grid,
    rho = config$rho_grid,
    replication_id = seq_len(config$rep_count),
    stringsAsFactors = FALSE
  )
}

fit_reference_posterior <- function(dat, config, seed) {
  cfg <- config$mcmc
  elapsed <- system.time({
    fit <- fit_mcmc_bma(
      X = dat$X_train,
      y = dat$y_train,
      group_id = dat$group_id,
      theta = 0.04,
      tau2 = 4,
      a0 = 1,
      b0 = 1,
      prior_type = "unrestricted",
      n_iter = cfg$n_iter,
      burn = cfg$burn,
      thin = cfg$thin,
      max_size = cfg$max_size,
      start_size = cfg$start_size,
      n_chains = cfg$n_chains,
      seed = seed
    )
  })[["elapsed"]]
  fit$sampler_label <- if (cfg$n_iter >= 10000L) {
    "long multi-chain MC3"
  } else {
    "multi-chain MC3"
  }
  fit <- enrich_mcmc_reference_fit(fit, dat$X_train, dat$y_train)
  fit$runtime_sec_reference <- elapsed
  fit$diagnostics <- mcmc_reference_diagnostics(fit, dat$group_id, elapsed)
  fit
}

add_reference_columns <- function(rows, diag) {
  for (nm in names(diag)) {
    rows[[nm]] <- diag[[nm]][1]
  }
  rows
}

run_large_cell <- function(scenario, rho, replication_id, config) {
  seed <- 20260528L + replication_id + round(1000 * rho) + 10000L * match(scenario, config$scenario_grid)
  dat <- simulate_large_askpc_data(scenario, rho, seed, config)
  ref <- fit_reference_posterior(dat, config, seed + 1000L)
  ref_path <- file.path(
    object_dir,
    tag_file(sprintf("%s_%s_rho%s_rep%03d_reference.rds",
                     config$mode, scenario, gsub("\\.", "p", sprintf("%.2f", rho)), replication_id),
             tag = config$run_tag)
  )
  saveRDS(list(data = dat, reference_fit = ref), ref_path)
  reference_pred <- predict_mixture(
    ref,
    X_train = dat$X_train,
    y_train = dat$y_train,
    X_test = dat$X_test,
    y_test = dat$y_test,
    weights = ref$posterior
  )
  rows <- list()
  add_row <- function(x) rows[[length(rows) + 1L]] <<- x
  add_failure <- function(method, method_family, tuning, beta, error) {
    add_row(failure_summary(method, method_family, tuning, beta, error))
  }

  ask_mode <- if (identical(config$mode, "full")) "medium" else config$mode
  for (beta in config$beta_grid) {
    tuning_beta <- paste0("beta=", signif(beta, 3))
    ask <- tryCatch(
      fit_askpc_pooled_pruned(
        supports = ref$supports,
        weights = ref$posterior,
        group_id = dat$group_id,
        X = dat$X_train,
        beta = beta,
        tau = 1e-3,
        q0_min = 1e-3,
        prune_mass = 0.99,
        mode = ask_mode,
        max_iter = if (identical(config$mode, "smoke")) 80L else 140L,
        fit_mode = "exact_weighted"
      ),
      error = function(e) e
    )
    if (inherits(ask, "error")) {
      add_failure("ASK-PC pooled union", "support-kernel compression", tuning_beta, beta, ask)
      add_failure("ASK-PC pooled-pruned 99%", "support-kernel compression", tuning_beta, beta, ask)
    } else {
      add_row(evaluate_kernel_summary(
        method = "ASK-PC pooled union",
        reference_fit = ref,
        dat = dat,
        W = ask$pooled$W,
        alpha = ask$pooled$alpha,
        q = ask$pooled$fit$q,
        costs = ask$pooled$costs,
        reference_pred = reference_pred,
        beta = beta,
        tuning = tuning_beta,
        method_family = "support-kernel compression",
        objective = ask$pooled$fit$objective,
        kkt_residual = ask$pooled$fit$kkt_residual,
        optimizer_status = ask$pooled$fit$status,
        optimizer_iterations = ask$pooled$fit$iterations,
        objective_change = ask$pooled$fit$objective_change,
        raw_candidates = ask$pool_diagnostics$raw_candidates,
        key_deduped_candidates = ask$pool_diagnostics$key_deduped_candidates,
        near_deduped_candidates = ask$pool_diagnostics$near_deduped_candidates,
        objective_consistent = ask$pool_diagnostics$objective_consistent
      ))
      add_row(evaluate_kernel_summary(
        method = "ASK-PC pooled-pruned 99%",
        reference_fit = ref,
        dat = dat,
        W = ask$pruned$W,
        alpha = ask$pruned$alpha,
        q = ask$pruned$fit$q,
        costs = ask$pruned$costs,
        reference_pred = reference_pred,
        beta = beta,
        tuning = tuning_beta,
        method_family = "support-kernel compression",
        objective = ask$pruned$fit$objective,
        kkt_residual = ask$pruned$fit$kkt_residual,
        optimizer_status = ask$pruned$fit$status,
        optimizer_iterations = ask$pruned$fit$iterations,
        objective_change = ask$pruned$fit$objective_change,
        raw_candidates = ask$pool_diagnostics$raw_candidates,
        key_deduped_candidates = ask$pool_diagnostics$key_deduped_candidates,
        near_deduped_candidates = ask$pool_diagnostics$near_deduped_candidates,
        objective_consistent = ask$pool_diagnostics$objective_consistent
      ))
    }

    fixed <- tryCatch(fit_fixed_hard_benchmark(dat, ref, beta = beta, capacity = 1L), error = function(e) e)
    if (inherits(fixed, "error")) {
      add_failure("Fixed hard dictionary", "fixed region dictionary", tuning_beta, beta, fixed)
    } else {
      add_row(evaluate_kernel_summary(
        method = "Fixed hard dictionary",
        reference_fit = ref,
        dat = dat,
        W = fixed$W,
        alpha = fixed$alpha,
        q = fixed$q,
        costs = fixed$costs,
        reference_pred = reference_pred,
        beta = beta,
        tuning = tuning_beta,
        method_family = "fixed region dictionary",
        objective = fixed$fit$objective,
        kkt_residual = fixed$fit$kkt_residual,
        optimizer_status = fixed$fit$status,
        optimizer_iterations = fixed$fit$iterations,
        objective_change = fixed$fit$objective_change
      ))
    }

    for (top_m in config$topm_grid) {
      tuning_topm <- paste0("M=", top_m, ", beta=", signif(beta, 3))
      topm <- tryCatch(fit_topm_benchmark(ref, beta = beta, top_m = top_m), error = function(e) e)
      if (inherits(topm, "error")) {
        add_failure("Top-M support atoms", "atom truncation", tuning_topm, beta, topm)
      } else {
        z <- evaluate_kernel_summary(
          method = "Top-M support atoms",
          reference_fit = ref,
          dat = dat,
          W = topm$W,
          alpha = topm$alpha,
          q = topm$q,
          costs = topm$costs,
          reference_pred = reference_pred,
          beta = beta,
          tuning = tuning_topm,
          method_family = "atom truncation",
          objective = topm$fit$objective,
          kkt_residual = topm$fit$kkt_residual,
          optimizer_status = topm$fit$status,
          optimizer_iterations = topm$fit$iterations,
          objective_change = topm$fit$objective_change
        )
        z$stored_atoms <- top_m
        add_row(z)
      }
    }

    for (coverage in config$credible_grid) {
      tuning_credible <- paste0("coverage=", signif(coverage, 3), ", beta=", signif(beta, 3))
      credible <- tryCatch(fit_credible_support_benchmark(ref, beta = beta, coverage = coverage), error = function(e) e)
      if (inherits(credible, "error")) {
        add_failure("Credible support set", "credible support summary", tuning_credible, beta, credible)
      } else {
        z <- evaluate_kernel_summary(
          method = "Credible support set",
          reference_fit = ref,
          dat = dat,
          W = credible$W,
          alpha = credible$alpha,
          q = credible$q,
          costs = credible$costs,
          reference_pred = reference_pred,
          beta = beta,
          tuning = tuning_credible,
          method_family = "credible support summary",
          objective = credible$fit$objective,
          kkt_residual = credible$fit$kkt_residual,
          optimizer_status = credible$fit$status,
          optimizer_iterations = credible$fit$iterations,
          objective_change = credible$fit$objective_change
        )
        z$stored_atoms <- credible$stored_atoms
        add_row(z)
      }
    }

    for (k in config$cluster_grid) {
      tuning_cluster <- paste0("clusters=", k, ", beta=", signif(beta, 3))
      cl <- tryCatch(
        fit_posterior_clustering_benchmark(dat, ref, beta = beta, n_clusters = k, seed = seed + k),
        error = function(e) e
      )
      if (inherits(cl, "error")) {
        add_failure("Posterior clustering", "posterior clustering", tuning_cluster, beta, cl)
      } else {
        add_row(evaluate_kernel_summary(
          method = "Posterior clustering",
          reference_fit = ref,
          dat = dat,
          W = cl$W,
          alpha = cl$alpha,
          q = cl$q,
          costs = cl$costs,
          reference_pred = reference_pred,
          beta = beta,
          tuning = tuning_cluster,
          method_family = "posterior clustering",
          objective = cl$fit$objective,
          kkt_residual = cl$fit$kkt_residual,
          optimizer_status = cl$fit$status,
          optimizer_iterations = cl$fit$iterations,
          objective_change = cl$fit$objective_change
        ))
      }
    }
  }

  for (power in config$dilution_power_grid) {
    tuning_dilution <- paste0("power=", signif(power, 3), " over visited supports")
    dil <- tryCatch(fit_dilution_benchmark(dat, ref, power = power), error = function(e) e)
    if (inherits(dil, "error")) {
      add_failure("Dilution-prior BMA", "redundancy-aware prior", tuning_dilution, NA_real_, dil)
    } else {
      support_sum <- posterior_summary(dil$posterior)
      add_row(evaluate_weight_summary(
        method = "Dilution-prior BMA",
        reference_fit = ref,
        approx_weights = dil$posterior,
        dat = dat,
        reference_pred = reference_pred,
        expected_code = support_sum$n_mass,
        n_kernels = NA_real_,
        stored_atoms = support_sum$n_mass,
        objective = NA_real_,
        tuning = tuning_dilution,
        method_family = "redundancy-aware prior"
      ))
    }
  }

  for (power in config$dpp_power_grid) {
    tuning_dpp <- paste0("power=", signif(power, 3), " over visited supports")
    dpp <- tryCatch(fit_dpp_benchmark(dat, ref, power = power), error = function(e) e)
    if (inherits(dpp, "error")) {
      add_failure("DPP-prior BMA", "redundancy-aware prior", tuning_dpp, NA_real_, dpp)
    } else {
      support_sum <- posterior_summary(dpp$posterior)
      add_row(evaluate_weight_summary(
        method = "DPP-prior BMA",
        reference_fit = ref,
        approx_weights = dpp$posterior,
        dat = dat,
        reference_pred = reference_pred,
        expected_code = support_sum$n_mass,
        n_kernels = NA_real_,
        stored_atoms = support_sum$n_mass,
        objective = NA_real_,
        tuning = tuning_dpp,
        method_family = "redundancy-aware prior"
      ))
    }
  }

  predictive_specs <- list(
    list(method = "Lasso", alpha = 1),
    list(method = "Elastic net", alpha = 0.5)
  )
  for (spec in predictive_specs) {
    pred_fit <- tryCatch(
      fit_glmnet_predictive_benchmark(dat, alpha = spec$alpha, method = spec$method, seed = seed + round(100 * spec$alpha)),
      error = function(e) e
    )
    if (inherits(pred_fit, "error")) {
      add_failure(spec$method, "predictive selection", paste0("alpha=", spec$alpha), NA_real_, pred_fit)
    } else {
      add_row(evaluate_predictive_summary(
        method = spec$method,
        dat = dat,
        reference_pred = reference_pred,
        pred = pred_fit$pred,
        selected_count = pred_fit$selected_count,
        runtime_sec = pred_fit$runtime_sec,
        tuning = pred_fit$tuning
      ))
    }
  }
  group_pred <- tryCatch(fit_group_lasso_predictive_benchmark(dat, seed = seed + 333L), error = function(e) e)
  if (inherits(group_pred, "error")) {
    add_failure("Group lasso", "predictive selection", "group lasso", NA_real_, group_pred)
  } else {
    add_row(evaluate_predictive_summary(
      method = "Group lasso",
      dat = dat,
      reference_pred = reference_pred,
      pred = group_pred$pred,
      selected_count = group_pred$selected_count,
      runtime_sec = group_pred$runtime_sec,
      tuning = group_pred$tuning
    ))
  }

  out <- do.call(rbind, rows)
  out$scenario <- scenario
  out$rho <- rho
  out$replication_id <- replication_id
  out$mode <- config$mode
  out$p <- ncol(dat$X_train)
  out$n_train <- nrow(dat$X_train)
  out$n_test <- nrow(dat$X_test)
  out$reference_object <- ref_path
  out$reference_support_entropy <- posterior_summary(ref$posterior)$entropy
  out$reference_support_n95 <- posterior_summary(ref$posterior)$n_mass
  add_reference_columns(out, ref$diagnostics)
}

append_csv <- function(path, rows) {
  write.table(
    rows,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(path),
    append = file.exists(path),
    quote = TRUE
  )
}

summarize_large_detail <- function(detail) {
  metrics <- c(
    "tv", "fkl", "rkl", "q0", "expected_code", "active_kernels_001",
    "q_effective_kernels", "n_kernels", "stored_atoms", "objective",
    "rmse_gap", "logscore_gap", "kkt_residual", "runtime_sec_reference",
    "optimizer_iterations", "objective_change", "raw_candidates",
    "key_deduped_candidates", "near_deduped_candidates",
    "reference_support_n95", "unique_supports", "ess_min", "split_rhat_max"
  )
  metrics <- intersect(metrics, names(detail))
  groups <- c("scenario", "rho", "method", "method_family", "tuning", "mode", "p", "reference_mode")
  groups <- intersect(groups, names(detail))
  se <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  }
  parts <- split(detail, interaction(detail[groups], drop = TRUE), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(d) {
    key <- d[1, groups, drop = FALSE]
    vals <- unlist(lapply(metrics, function(m) c(mean = mean(d[[m]], na.rm = TRUE), se = se(d[[m]]))))
    names(vals) <- paste0(rep(metrics, each = 2L), rep(c("_mean", "_se"), length(metrics)))
    data.frame(key, as.list(vals), n_replications = length(unique(d$replication_id)), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

best_fkl_summary <- function(summary) {
  ok <- is.finite(summary$fkl_mean)
  z <- summary[ok, , drop = FALSE]
  if (!nrow(z)) return(z)
  parts <- split(z, interaction(z$scenario, z$method, drop = TRUE), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(d) d[which.min(d$fkl_mean), , drop = FALSE]))
  rownames(out) <- NULL
  out[order(out$scenario, out$fkl_mean), , drop = FALSE]
}

write_outputs <- function(detail, mode, publish = FALSE) {
  detail <- normalize_reference_diagnostics(detail)
  detail_path <- file.path(table_dir, tag_file("table_large_end_to_end_askpc_detail.csv"))
  summary_path <- file.path(table_dir, tag_file("table_large_end_to_end_askpc_summary.csv"))
  best_path <- file.path(table_dir, tag_file("table_large_end_to_end_askpc_best_fkl.csv"))
  diag_path <- file.path(table_dir, tag_file("table_large_end_to_end_reference_diagnostics.csv"))
  write.csv(detail, detail_path, row.names = FALSE)
  summary <- summarize_large_detail(detail)
  write.csv(summary, summary_path, row.names = FALSE)
  write.csv(best_fkl_summary(summary), best_path, row.names = FALSE)
  diag_cols <- c(
    "scenario", "rho", "replication_id", "mode", "p", "n_train", "n_test",
    "reference_object", "reference_mode", "sampler", "n_iter", "burn", "thin", "n_chains",
    "retained_draws", "unique_supports", "unique_group_states",
    "top_support_mass", "effective_support_count", "mean_size",
    "accept_rate", "group_pip_max_range", "group_pip_mean_range",
    "split_rhat_max", "split_rhat_median", "ess_min", "ess_median",
    "group_pip_mcse_max", "runtime_sec_reference", "diagnostic_status",
    "diagnostic_warning"
  )
  diag_cols <- intersect(diag_cols, names(detail))
  diag <- unique(detail[, diag_cols, drop = FALSE])
  write.csv(diag, diag_path, row.names = FALSE)
  if (publish && nzchar(run_tag)) {
    write.csv(detail, file.path(table_dir, "table_large_end_to_end_askpc_detail.csv"), row.names = FALSE)
    write.csv(summary, file.path(table_dir, "table_large_end_to_end_askpc_summary.csv"), row.names = FALSE)
    write.csv(best_fkl_summary(summary), file.path(table_dir, "table_large_end_to_end_askpc_best_fkl.csv"), row.names = FALSE)
    write.csv(diag, file.path(table_dir, "table_large_end_to_end_reference_diagnostics.csv"), row.names = FALSE)
  }
  list(detail = detail, summary = summary, best = best_fkl_summary(summary), diagnostics = diag)
}

run_large_benchmark <- function(mode, checkpoint = TRUE, resume = TRUE) {
  config <- override_config(large_mode_config(mode))
  shard_id <- as.integer(get_arg("shard-id", "1"))
  shard_total <- as.integer(get_arg("shard-total", "1"))
  if (is.na(shard_id) || is.na(shard_total) || shard_id < 1L || shard_total < 1L || shard_id > shard_total) {
    stop("Shard arguments must satisfy 1 <= shard-id <= shard-total")
  }
  checkpoint_name <- if (shard_total > 1L) {
    sprintf("table_large_end_to_end_askpc_detail_%s_shard%02dof%02d_checkpoint.csv", mode, shard_id, shard_total)
  } else {
    sprintf("table_large_end_to_end_askpc_detail_%s_checkpoint.csv", mode)
  }
  checkpoint_name <- tag_file(checkpoint_name)
  checkpoint_path <- file.path(table_dir, checkpoint_name)
  cells <- benchmark_grid(config)
  cells$key <- cell_key(cells$scenario, cells$rho, cells$replication_id)
  cells <- cells[((seq_len(nrow(cells)) - 1L) %% shard_total) + 1L == shard_id, , drop = FALSE]
  done <- character(0)
  if (checkpoint && resume && file.exists(checkpoint_path)) {
    old <- read.csv(checkpoint_path, stringsAsFactors = FALSE)
    done <- unique(cell_key(old$scenario, old$rho, old$replication_id))
    message(sprintf("resuming from %s with %d completed cells", checkpoint_path, length(done)))
  }
  remaining <- cells[!(cells$key %in% done), , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(remaining))) {
    z <- remaining[i, , drop = FALSE]
    message(sprintf(
      "large ASK-PC cell %d of %d: %s rho=%.2f rep=%d p=%d",
      i, nrow(remaining), z$scenario, z$rho, z$replication_id, config$p
    ))
    cell <- tryCatch(
      run_large_cell(z$scenario, z$rho, z$replication_id, config),
      error = function(e) {
        out <- cell_failure_result(z$scenario, z$rho, z$replication_id, config, e)
        out$p <- config$p
        out$n_train <- config$n_train
        out$n_test <- config$n_test
        out$reference_mode <- "mcmc_empirical_visited_supports"
        out$diagnostic_status <- "cell failed"
        out
      }
    )
    rows[[length(rows) + 1L]] <- cell
    if (checkpoint) append_csv(checkpoint_path, cell)
  }
  if (checkpoint && file.exists(checkpoint_path)) {
    detail <- read.csv(checkpoint_path, stringsAsFactors = FALSE)
    if (mode != "full" || shard_total == 1L) {
      write_outputs(detail, mode, publish = has_flag("publish"))
    }
    return(list(detail = detail))
  }
  detail <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write_outputs(detail, mode, publish = has_flag("publish"))
}

if (identical(mode_arg, "summarize")) {
  detail_path <- file.path(table_dir, tag_file("table_large_end_to_end_askpc_detail.csv"))
  checkpoint_path <- get_arg("checkpoint-path", NULL)
  if (!is.null(checkpoint_path)) {
    detail_path <- checkpoint_path
  }
  if (is.null(checkpoint_path)) {
    full_pattern <- if (nzchar(run_tag)) {
      sprintf("table_large_end_to_end_askpc_detail_full_shard*_checkpoint_%s.csv", run_tag)
    } else {
      "table_large_end_to_end_askpc_detail_full_shard*_checkpoint.csv"
    }
    full_shards <- sort(Sys.glob(file.path(table_dir, full_pattern)))
    candidates <- if (length(full_shards)) {
      full_shards
    } else if (!file.exists(detail_path)) {
      fallback_pattern <- if (nzchar(run_tag)) {
        sprintf("table_large_end_to_end_askpc_detail_*_checkpoint_%s.csv", run_tag)
      } else {
        "table_large_end_to_end_askpc_detail_*_checkpoint.csv"
      }
      sort(Sys.glob(file.path(table_dir, fallback_pattern)))
    } else {
      character()
    }
  } else {
    candidates <- character()
  }
  if (length(candidates)) {
    if (!length(candidates)) stop("No large end-to-end detail or checkpoint file found")
    detail <- do.call(rbind, lapply(candidates, function(path) read.csv(path, stringsAsFactors = FALSE)))
    out <- write_outputs(detail, mode = "summarize", publish = has_flag("publish"))
    message(sprintf("Summarized %d rows from %d checkpoint files", nrow(detail), length(candidates)))
    message(sprintf("Summary rows: %d", nrow(out$summary)))
    quit(save = "no")
  }
  detail <- read.csv(detail_path, stringsAsFactors = FALSE)
  out <- write_outputs(detail, mode = "summarize", publish = has_flag("publish"))
  message(sprintf("Summarized %d rows from %s", nrow(detail), detail_path))
  message(sprintf("Summary rows: %d", nrow(out$summary)))
  quit(save = "no")
}

checkpoint <- !has_flag("no-checkpoint")
resume <- !has_flag("no-resume")
out <- run_large_benchmark(mode_arg, checkpoint = checkpoint, resume = resume)
message("Large end-to-end ASK-PC benchmark wrote:")
message("  sim/output/tables/table_large_end_to_end_askpc_detail.csv")
message("  sim/output/tables/table_large_end_to_end_askpc_summary.csv")
message("  sim/output/tables/table_large_end_to_end_askpc_best_fkl.csv")
message("  sim/output/tables/table_large_end_to_end_reference_diagnostics.csv")
message(sprintf("Rows currently available: %d", nrow(out$detail)))
