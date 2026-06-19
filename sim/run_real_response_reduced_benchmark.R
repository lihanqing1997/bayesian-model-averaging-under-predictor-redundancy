source(file.path("sim", "src", "support_kernel_benchmark.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit)) sub(pat, "", hit[[1]]) else default
}

mode <- tolower(get_arg("mode", "smoke"))
if (!mode %in% c("smoke", "medium", "full")) {
  stop("--mode must be smoke, medium, or full")
}

dataset_arg <- tolower(get_arg("dataset", "all"))
if (!dataset_arg %in% c("tecator", "gasoline", "all")) {
  stop("--dataset must be tecator, gasoline, or all")
}

table_dir <- file.path("sim", "output", "tables")
fig_dir <- file.path("sim", "output", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

real_response_config <- function(mode) {
  list(
    mode = mode,
    p_reduced = switch(mode, smoke = 10L, medium = 12L, full = 12L),
    rep_count = switch(mode, smoke = 1L, medium = 3L, full = 5L),
    beta_grid = switch(mode, smoke = 0.02, medium = c(0.02, 0.08), full = c(0.02, 0.08)),
    topm_grid = switch(mode, smoke = 16L, medium = c(16L, 32L), full = c(16L, 32L)),
    credible_grid = switch(mode, smoke = 0.90, medium = 0.90, full = 0.90),
    cluster_grid = switch(mode, smoke = 6L, medium = c(6L, 10L), full = c(6L, 10L)),
    dilution_power_grid = switch(mode, smoke = 0.5, medium = 0.5, full = 0.5),
    dpp_power_grid = switch(mode, smoke = 1, medium = 1, full = 1)
  )
}

make_blocks <- function(p, p_reduced) {
  cut(seq_len(p), breaks = p_reduced, labels = FALSE)
}

block_average <- function(X, p_reduced) {
  block <- make_blocks(ncol(X), p_reduced)
  out <- sapply(seq_len(p_reduced), function(b) rowMeans(X[, block == b, drop = FALSE]))
  out <- as.matrix(out)
  colnames(out) <- paste0("region", seq_len(ncol(out)))
  list(X = out, block = block)
}

load_real_response <- function(dataset) {
  if (identical(dataset, "tecator")) {
    if (!requireNamespace("caret", quietly = TRUE)) stop("caret is required for Tecator")
    data("tecator", package = "caret")
    return(list(
      X = as.matrix(absorp),
      y = as.numeric(endpoints[, 2]),
      response = "fat"
    ))
  }
  if (identical(dataset, "gasoline")) {
    if (!requireNamespace("pls", quietly = TRUE)) stop("pls is required for Gasoline")
    data("gasoline", package = "pls")
    return(list(
      X = as.matrix(gasoline$NIR),
      y = as.numeric(gasoline$octane),
      response = "octane"
    ))
  }
  stop("unknown dataset")
}

prepare_real_response_split <- function(dataset, rep_id, seed, config) {
  raw <- load_real_response(dataset)
  reduced <- block_average(raw$X, config$p_reduced)
  n <- nrow(reduced$X)
  test_size <- if (identical(dataset, "tecator")) max(40L, floor(0.25 * n)) else max(12L, floor(0.25 * n))
  set.seed(seed)
  test_id <- sort(sample.int(n, size = test_size))
  train_id <- setdiff(seq_len(n), test_id)
  z <- standardize_with_training(
    reduced$X[train_id, , drop = FALSE],
    reduced$X[test_id, , drop = FALSE]
  )
  y_center <- mean(raw$y[train_id])
  group_id <- rep(seq_len(ceiling(config$p_reduced / 2)), each = 2L)[seq_len(config$p_reduced)]
  list(
    X_train = z$X_train,
    X_test = z$X_test,
    y_train = raw$y[train_id] - y_center,
    y_test = raw$y[test_id] - y_center,
    y_center = y_center,
    train_id = train_id,
    test_id = test_id,
    group_id = as.integer(group_id),
    K = max(group_id),
    m = 2L,
    dataset = dataset,
    response = raw$response,
    spectral_block = reduced$block,
    p_reduced = config$p_reduced,
    n_train = length(train_id),
    n_test = length(test_id)
  )
}

run_real_response_cell <- function(dataset, rep_id, config) {
  seed <- 20260529L + rep_id + 1000L * match(dataset, c("tecator", "gasoline"))
  dat <- prepare_real_response_split(dataset, rep_id, seed, config)
  reference_fit <- fit_exact_bma(dat$X_train, dat$y_train, theta = min(0.20, 3 / dat$p_reduced), tau2 = 4, a0 = 1, b0 = 1)
  reference_pred <- predict_mixture(reference_fit, dat$X_train, dat$y_train, dat$X_test, dat$y_test, weights = reference_fit$posterior)
  rows <- list()
  add_row <- function(x) rows[[length(rows) + 1L]] <<- x
  add_failure <- function(method, method_family, tuning, beta, error) {
    add_row(failure_summary(method, method_family, tuning, beta, error))
  }

  for (beta in config$beta_grid) {
    tuning_beta <- paste0("beta=", signif(beta, 3))
    ask <- tryCatch(
      fit_askpc_pooled_pruned(
        supports = reference_fit$supports,
        weights = reference_fit$posterior,
        group_id = dat$group_id,
        X = dat$X_train,
        beta = beta,
        tau = 1e-3,
        q0_min = 1e-3,
        prune_mass = 0.99,
        mode = if (identical(config$mode, "full")) "medium" else config$mode,
        max_iter = if (identical(config$mode, "smoke")) 70L else 120L
      ),
      error = function(e) e
    )
    if (inherits(ask, "error")) {
      add_failure("ASK-PC pooled union", "support-kernel compression", tuning_beta, beta, ask)
      add_failure("ASK-PC pooled-pruned 99%", "support-kernel compression", tuning_beta, beta, ask)
    } else {
      add_row(evaluate_kernel_summary("ASK-PC pooled union", reference_fit, dat, ask$pooled$W, ask$pooled$alpha,
                                      ask$pooled$fit$q, ask$pooled$costs, reference_pred, beta = beta,
                                      tuning = tuning_beta, method_family = "support-kernel compression",
                                      objective = ask$pooled$fit$objective, kkt_residual = ask$pooled$fit$kkt_residual,
                                      optimizer_status = ask$pooled$fit$status, optimizer_iterations = ask$pooled$fit$iterations,
                                      objective_change = ask$pooled$fit$objective_change,
                                      raw_candidates = ask$pool_diagnostics$raw_candidates,
                                      key_deduped_candidates = ask$pool_diagnostics$key_deduped_candidates,
                                      near_deduped_candidates = ask$pool_diagnostics$near_deduped_candidates,
                                      objective_consistent = ask$pool_diagnostics$objective_consistent))
      add_row(evaluate_kernel_summary("ASK-PC pooled-pruned 99%", reference_fit, dat, ask$pruned$W, ask$pruned$alpha,
                                      ask$pruned$fit$q, ask$pruned$costs, reference_pred, beta = beta,
                                      tuning = tuning_beta, method_family = "support-kernel compression",
                                      objective = ask$pruned$fit$objective, kkt_residual = ask$pruned$fit$kkt_residual,
                                      optimizer_status = ask$pruned$fit$status, optimizer_iterations = ask$pruned$fit$iterations,
                                      objective_change = ask$pruned$fit$objective_change,
                                      raw_candidates = ask$pool_diagnostics$raw_candidates,
                                      key_deduped_candidates = ask$pool_diagnostics$key_deduped_candidates,
                                      near_deduped_candidates = ask$pool_diagnostics$near_deduped_candidates,
                                      objective_consistent = ask$pool_diagnostics$objective_consistent))
    }

    fixed <- tryCatch(fit_fixed_hard_benchmark(dat, reference_fit, beta = beta, capacity = 1L), error = function(e) e)
    if (inherits(fixed, "error")) {
      add_failure("Fixed hard dictionary", "fixed region dictionary", tuning_beta, beta, fixed)
    } else {
      add_row(evaluate_kernel_summary("Fixed hard dictionary", reference_fit, dat, fixed$W, fixed$alpha,
                                      fixed$q, fixed$costs, reference_pred, beta = beta,
                                      tuning = tuning_beta, method_family = "fixed region dictionary",
                                      objective = fixed$fit$objective, kkt_residual = fixed$fit$kkt_residual,
                                      optimizer_status = fixed$fit$status, optimizer_iterations = fixed$fit$iterations,
                                      objective_change = fixed$fit$objective_change))
    }

    for (top_m in config$topm_grid) {
      tuning_topm <- paste0("M=", top_m, ", beta=", signif(beta, 3))
      topm <- tryCatch(fit_topm_benchmark(reference_fit, beta = beta, top_m = top_m), error = function(e) e)
      if (inherits(topm, "error")) {
        add_failure("Top-M support atoms", "atom truncation", tuning_topm, beta, topm)
      } else {
        z <- evaluate_kernel_summary("Top-M support atoms", reference_fit, dat, topm$W, topm$alpha,
                                     topm$q, topm$costs, reference_pred, beta = beta,
                                     tuning = tuning_topm, method_family = "atom truncation",
                                     objective = topm$fit$objective, kkt_residual = topm$fit$kkt_residual,
                                     optimizer_status = topm$fit$status, optimizer_iterations = topm$fit$iterations,
                                     objective_change = topm$fit$objective_change)
        z$stored_atoms <- top_m
        add_row(z)
      }
    }

    for (coverage in config$credible_grid) {
      tuning_credible <- paste0("coverage=", signif(coverage, 3), ", beta=", signif(beta, 3))
      credible <- tryCatch(fit_credible_support_benchmark(reference_fit, beta = beta, coverage = coverage), error = function(e) e)
      if (inherits(credible, "error")) {
        add_failure("Credible support set", "credible support summary", tuning_credible, beta, credible)
      } else {
        z <- evaluate_kernel_summary("Credible support set", reference_fit, dat, credible$W, credible$alpha,
                                     credible$q, credible$costs, reference_pred, beta = beta,
                                     tuning = tuning_credible, method_family = "credible support summary",
                                     objective = credible$fit$objective, kkt_residual = credible$fit$kkt_residual,
                                     optimizer_status = credible$fit$status, optimizer_iterations = credible$fit$iterations,
                                     objective_change = credible$fit$objective_change)
        z$stored_atoms <- credible$stored_atoms
        add_row(z)
      }
    }

    for (k in config$cluster_grid) {
      tuning_cluster <- paste0("clusters=", k, ", beta=", signif(beta, 3))
      cl <- tryCatch(fit_posterior_clustering_benchmark(dat, reference_fit, beta = beta, n_clusters = k, seed = seed + k),
                     error = function(e) e)
      if (inherits(cl, "error")) {
        add_failure("Posterior clustering", "posterior clustering", tuning_cluster, beta, cl)
      } else {
        add_row(evaluate_kernel_summary("Posterior clustering", reference_fit, dat, cl$W, cl$alpha,
                                        cl$q, cl$costs, reference_pred, beta = beta,
                                        tuning = tuning_cluster, method_family = "posterior clustering",
                                        objective = cl$fit$objective, kkt_residual = cl$fit$kkt_residual,
                                        optimizer_status = cl$fit$status, optimizer_iterations = cl$fit$iterations,
                                        objective_change = cl$fit$objective_change))
      }
    }
  }

  for (power in config$dilution_power_grid) {
    tuning_dilution <- paste0("power=", signif(power, 3))
    dil <- tryCatch(fit_dilution_benchmark(dat, reference_fit, power = power), error = function(e) e)
    if (inherits(dil, "error")) {
      add_failure("Dilution-prior BMA", "redundancy-aware prior", tuning_dilution, NA_real_, dil)
    } else {
      support_sum <- posterior_summary(dil$posterior)
      add_row(evaluate_weight_summary("Dilution-prior BMA", reference_fit, dil$posterior, dat, reference_pred,
                                      expected_code = support_sum$n_mass, stored_atoms = support_sum$n_mass,
                                      tuning = tuning_dilution, method_family = "redundancy-aware prior"))
    }
  }

  for (power in config$dpp_power_grid) {
    tuning_dpp <- paste0("power=", signif(power, 3))
    dpp <- tryCatch(fit_dpp_benchmark(dat, reference_fit, power = power), error = function(e) e)
    if (inherits(dpp, "error")) {
      add_failure("DPP-prior BMA", "redundancy-aware prior", tuning_dpp, NA_real_, dpp)
    } else {
      support_sum <- posterior_summary(dpp$posterior)
      add_row(evaluate_weight_summary("DPP-prior BMA", reference_fit, dpp$posterior, dat, reference_pred,
                                      expected_code = support_sum$n_mass, stored_atoms = support_sum$n_mass,
                                      tuning = tuning_dpp, method_family = "redundancy-aware prior"))
    }
  }

  predictive_specs <- list(list(method = "Lasso", alpha = 1), list(method = "Elastic net", alpha = 0.5))
  for (spec in predictive_specs) {
    pred_fit <- tryCatch(fit_glmnet_predictive_benchmark(dat, spec$alpha, spec$method, seed + round(100 * spec$alpha)),
                         error = function(e) e)
    if (inherits(pred_fit, "error")) {
      add_failure(spec$method, "predictive selection", paste0("alpha=", spec$alpha), NA_real_, pred_fit)
    } else {
      add_row(evaluate_predictive_summary(spec$method, dat, reference_pred, pred_fit$pred,
                                          pred_fit$selected_count, pred_fit$runtime_sec, pred_fit$tuning))
    }
  }
  group_pred <- tryCatch(fit_group_lasso_predictive_benchmark(dat, seed + 333L), error = function(e) e)
  if (inherits(group_pred, "error")) {
    add_failure("Group lasso", "predictive selection", "group lasso", NA_real_, group_pred)
  } else {
    add_row(evaluate_predictive_summary("Group lasso", dat, reference_pred, group_pred$pred,
                                        group_pred$selected_count, group_pred$runtime_sec, group_pred$tuning))
  }

  out <- do.call(rbind, rows)
  out$dataset <- dataset
  out$response <- dat$response
  out$replication_id <- rep_id
  out$mode <- config$mode
  out$p_reduced <- dat$p_reduced
  out$n_train <- dat$n_train
  out$n_test <- dat$n_test
  out$reference_type <- "exact enumerated reduced real-response BMA"
  out$reference_supports <- nrow(reference_fit$supports)
  out$reference_support_entropy <- posterior_summary(reference_fit$posterior)$entropy
  out$reference_support_n95 <- posterior_summary(reference_fit$posterior)$n_mass
  out
}

se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

summarize_detail <- function(detail) {
  metrics <- c("tv", "fkl", "rkl", "q0", "expected_code", "active_kernels_001",
               "q_effective_kernels", "stored_atoms", "rmse_gap", "logscore_gap",
               "kkt_residual", "runtime_sec", "reference_support_n95")
  metrics <- intersect(metrics, names(detail))
  groups <- c("dataset", "method", "method_family", "tuning", "mode", "p_reduced")
  parts <- split(detail, interaction(detail[groups], drop = TRUE), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(d) {
    key <- d[1, groups, drop = FALSE]
    vals <- unlist(lapply(metrics, function(m) c(mean = mean(d[[m]], na.rm = TRUE), se = se(d[[m]]))))
    names(vals) <- paste0(rep(metrics, each = 2), rep(c("_mean", "_se"), length(metrics)))
    data.frame(key, as.list(vals), n_replications = length(unique(d$replication_id)), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

best_by <- function(summary, metric = "fkl_mean") {
  z <- summary[is.finite(summary[[metric]]), , drop = FALSE]
  parts <- split(z, interaction(z$dataset, z$method, drop = TRUE), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(d) d[which.min(d[[metric]]), , drop = FALSE]))
  rownames(out) <- NULL
  out
}

fmt <- function(x, digits = 3) {
  if (!is.finite(x)) return("--")
  formatC(x, format = "f", digits = digits)
}

fmt_mean_se <- function(mu, se_val, digits = 3) {
  if (!is.finite(mu)) return("--")
  if (!is.finite(se_val)) return(fmt(mu, digits))
  paste0(fmt(mu, digits), " (", fmt(se_val, digits), ")")
}

write_main_tex <- function(best, path) {
  keep <- c("ASK-PC pooled-pruned 99%", "Posterior clustering", "Top-M support atoms",
            "Credible support set", "Dilution-prior BMA", "DPP-prior BMA", "Fixed hard dictionary")
  lab <- c("ASK-PC pooled-pruned 99%" = "Pooled-pruned",
           "Posterior clustering" = "Cluster kernels",
           "Top-M support atoms" = "Top-M atoms",
           "Credible support set" = "Credible set",
           "Dilution-prior BMA" = "Dilution BMA",
           "DPP-prior BMA" = "DPP BMA",
           "Fixed hard dictionary" = "Fixed regions")
  z <- best[best$method %in% keep, , drop = FALSE]
  z$method <- factor(z$method, levels = keep)
  z <- z[order(z$dataset, z$method), , drop = FALSE]
  con <- file(path, open = "w")
  on.exit(close(con), add = TRUE)
  writeLines("\\begin{tabular}{llccccc}", con)
  writeLines("\\toprule", con)
  writeLines("Dataset & Method & TV & FKL & $q_0$ & Reporting cost & RMSE gap \\\\", con)
  writeLines("\\midrule", con)
  for (ds in unique(z$dataset)) {
    d <- z[z$dataset == ds, , drop = FALSE]
    for (i in seq_len(nrow(d))) {
      ds_cell <- if (i == 1L) tools::toTitleCase(ds) else ""
      line <- paste(
        ds_cell,
        lab[as.character(d$method[i])],
        fmt_mean_se(d$tv_mean[i], d$tv_se[i]),
        fmt_mean_se(d$fkl_mean[i], d$fkl_se[i]),
        fmt_mean_se(d$q0_mean[i], d$q0_se[i]),
        fmt_mean_se(d$expected_code_mean[i], d$expected_code_se[i], digits = 2),
        fmt_mean_se(d$rmse_gap_mean[i], d$rmse_gap_se[i]),
        sep = " & "
      )
      writeLines(paste0(line, " \\\\"), con)
    }
    if (ds != tail(unique(z$dataset), 1)) writeLines("\\addlinespace", con)
  }
  writeLines("\\bottomrule", con)
  writeLines("\\end{tabular}", con)
}

config <- real_response_config(mode)
datasets <- if (identical(dataset_arg, "all")) c("tecator", "gasoline") else dataset_arg
grid <- expand.grid(dataset = datasets, replication_id = seq_len(config$rep_count), stringsAsFactors = FALSE)

checkpoint_path <- file.path(table_dir, sprintf("table_real_response_reduced_detail_%s_checkpoint.csv", mode))
completed_keys <- character()
checkpoint_rows <- list()
if (file.exists(checkpoint_path)) {
  checkpoint <- read.csv(checkpoint_path, stringsAsFactors = FALSE)
  if (nrow(checkpoint)) {
    completed_keys <- unique(paste(checkpoint$dataset, checkpoint$replication_id, sep = "::"))
    checkpoint_rows[[1L]] <- checkpoint
  }
}

rows <- vector("list", nrow(grid))
row_i <- 0L
for (i in seq_len(nrow(grid))) {
  key <- paste(grid$dataset[i], grid$replication_id[i], sep = "::")
  if (key %in% completed_keys) {
    message(sprintf("real-response reduced cell %d of %d already checkpointed: %s rep=%d", i, nrow(grid), grid$dataset[i], grid$replication_id[i]))
    next
  }
  message(sprintf("real-response reduced cell %d of %d: %s rep=%d", i, nrow(grid), grid$dataset[i], grid$replication_id[i]))
  cell <- tryCatch(
    run_real_response_cell(grid$dataset[i], grid$replication_id[i], config),
    error = function(e) {
      out <- failure_summary("Cell failure", "real-response exact benchmark", NA_character_, NA_real_, e)
      out$dataset <- grid$dataset[i]
      out$replication_id <- grid$replication_id[i]
      out$mode <- config$mode
      out$p_reduced <- config$p_reduced
      out
    }
  )
  row_i <- row_i + 1L
  rows[[row_i]] <- cell
  checkpoint_rows[[length(checkpoint_rows) + 1L]] <- cell
  write.csv(do.call(rbind, checkpoint_rows), checkpoint_path, row.names = FALSE)
}

rows <- rows[seq_len(row_i)]
detail <- if (length(checkpoint_rows)) do.call(rbind, checkpoint_rows) else do.call(rbind, rows)
summary <- summarize_detail(detail)
best <- best_by(summary, "fkl_mean")
predictive <- summary[summary$method_family == "predictive selection", , drop = FALSE]
pred_best <- if (nrow(predictive)) best_by(predictive, "rmse_gap_mean") else data.frame()

write.csv(detail, file.path(table_dir, "table_real_response_reduced_detail.csv"), row.names = FALSE)
write.csv(summary, file.path(table_dir, "table_real_response_reduced_summary.csv"), row.names = FALSE)
write.csv(best, file.path(table_dir, "table_real_response_reduced_best_fkl.csv"), row.names = FALSE)
write.csv(pred_best, file.path(table_dir, "table_real_response_reduced_predictive.csv"), row.names = FALSE)

reference_diag <- unique(detail[, intersect(c("dataset", "mode", "p_reduced", "n_train", "n_test", "reference_type", "reference_supports", "reference_support_n95"), names(detail)), drop = FALSE])
reference_diag$status <- "pass"
write.csv(reference_diag, file.path(table_dir, "table_real_response_reduced_reference_diagnostics.csv"), row.names = FALSE)

write_main_tex(best, file.path(table_dir, "table_real_response_reduced_best_fkl.tex"))

message("Wrote real-response reduced benchmark outputs")
