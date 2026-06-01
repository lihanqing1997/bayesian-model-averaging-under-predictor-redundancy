source(file.path("sim", "src", "support_kernel_benchmark.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit)) sub(pat, "", hit[[1]]) else default
}

run_tag <- get_arg("run-tag", "reliable_full_v2")
out_tag <- get_arg("out-tag", paste0(run_tag, "_augmented"))
parse_grid <- function(name, default) {
  x <- get_arg(name, default)
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}
beta_grid_arg <- parse_grid("beta-grid", "0.02")
coverage_grid_arg <- parse_grid("coverage-grid", "0.9")
dpp_power_grid_arg <- parse_grid("dpp-power-grid", "0.5,1")

table_dir <- file.path("sim", "output", "tables")
object_dir <- file.path("sim", "output", "large_end_to_end")
base_detail_path <- file.path(table_dir, sprintf("table_large_end_to_end_askpc_detail_%s.csv", run_tag))
if (!file.exists(base_detail_path)) {
  stop("Missing base detail file: ", base_detail_path)
}
base_detail <- read.csv(base_detail_path, stringsAsFactors = FALSE)

drop_methods <- c("Credible support set", "DPP-prior BMA", "Lasso", "Elastic net", "Group lasso")
base_detail <- base_detail[!(base_detail$method %in% drop_methods), , drop = FALSE]

reference_paths <- sort(list.files(object_dir, pattern = paste0("_", run_tag, "\\.rds$"), full.names = TRUE))
if (!length(reference_paths)) {
  stop("No saved reference objects found for run tag: ", run_tag)
}

parse_reference_path <- function(path) {
  b <- basename(path)
  m <- regexec("^full_(.*)_rho([0-9]p[0-9]+)_rep([0-9]+)_reference_", b)
  r <- regmatches(b, m)[[1]]
  if (length(r) < 4L) stop("Cannot parse reference filename: ", b)
  data.frame(
    scenario = r[2],
    rho = as.numeric(gsub("p", ".", r[3], fixed = TRUE)),
    replication_id = as.integer(r[4]),
    stringsAsFactors = FALSE
  )
}

align_to_template <- function(row, template) {
  for (nm in names(template)) {
    if (!nm %in% names(row)) {
      row[[nm]] <- template[[nm]][1]
    }
  }
  for (nm in names(row)) {
    if (!nm %in% names(template)) {
      template[[nm]] <- NA
    }
  }
  row[, names(template), drop = FALSE]
}

checkpoint_path <- file.path(table_dir, sprintf("table_large_end_to_end_augmented_rows_%s.csv", out_tag))
new_rows <- list()
rid <- 0L
completed_keys <- character()
if (file.exists(checkpoint_path)) {
  ck <- read.csv(checkpoint_path, stringsAsFactors = FALSE)
  if (nrow(ck)) {
    new_rows[[1L]] <- ck
    rid <- 1L
    completed_keys <- unique(paste(ck$scenario, ck$rho, ck$replication_id, sep = "::"))
  }
}
for (path in reference_paths) {
  meta <- parse_reference_path(path)
  meta_key <- paste(meta$scenario, meta$rho, meta$replication_id, sep = "::")
  if (meta_key %in% completed_keys) {
    message("already checkpointed ", basename(path))
    next
  }
  obj <- readRDS(path)
  dat <- obj$data
  ref <- obj$reference_fit
  template <- base_detail[
    base_detail$scenario == meta$scenario &
      abs(base_detail$rho - meta$rho) < 1e-8 &
      base_detail$replication_id == meta$replication_id,
    ,
    drop = FALSE
  ]
  if (!nrow(template)) {
    warning("No template row for ", basename(path))
    next
  }
  template <- template[1, , drop = FALSE]
  reference_pred <- predict_mixture(ref, dat$X_train, dat$y_train, dat$X_test, dat$y_test, weights = ref$posterior)
  beta_grid <- beta_grid_arg
  path_rows <- list()
  path_rid <- 0L

  for (beta in beta_grid) {
    for (coverage in coverage_grid_arg) {
      tuning <- paste0("coverage=", signif(coverage, 3), ", beta=", signif(beta, 3))
      fit <- tryCatch(fit_credible_support_benchmark(ref, beta = beta, coverage = coverage), error = function(e) e)
      if (inherits(fit, "error")) {
        row <- failure_summary("Credible support set", "credible support summary", tuning, beta, fit)
      } else {
        row <- evaluate_kernel_summary(
          method = "Credible support set",
          reference_fit = ref,
          dat = dat,
          W = fit$W,
          alpha = fit$alpha,
          q = fit$q,
          costs = fit$costs,
          reference_pred = reference_pred,
          beta = beta,
          tuning = tuning,
          method_family = "credible support summary",
          objective = fit$fit$objective,
          kkt_residual = fit$fit$kkt_residual,
          optimizer_status = fit$fit$status,
          optimizer_iterations = fit$fit$iterations,
          objective_change = fit$fit$objective_change
        )
        row$stored_atoms <- fit$stored_atoms
      }
      path_rid <- path_rid + 1L
      path_rows[[path_rid]] <- align_to_template(row, template)
    }
  }

  for (power in dpp_power_grid_arg) {
    tuning <- paste0("power=", signif(power, 3), " over visited supports")
    fit <- tryCatch(fit_dpp_benchmark(dat, ref, power = power), error = function(e) e)
    if (inherits(fit, "error")) {
      row <- failure_summary("DPP-prior BMA", "redundancy-aware prior", tuning, NA_real_, fit)
    } else {
      support_sum <- posterior_summary(fit$posterior)
      row <- evaluate_weight_summary(
        method = "DPP-prior BMA",
        reference_fit = ref,
        approx_weights = fit$posterior,
        dat = dat,
        reference_pred = reference_pred,
        expected_code = support_sum$n_mass,
        stored_atoms = support_sum$n_mass,
        tuning = tuning,
        method_family = "redundancy-aware prior"
      )
    }
    path_rid <- path_rid + 1L
    path_rows[[path_rid]] <- align_to_template(row, template)
  }

  predictive_specs <- list(list(method = "Lasso", alpha = 1), list(method = "Elastic net", alpha = 0.5))
  for (spec in predictive_specs) {
    fit <- tryCatch(fit_glmnet_predictive_benchmark(dat, spec$alpha, spec$method, seed = 20260529L + meta$replication_id), error = function(e) e)
    if (inherits(fit, "error")) {
      row <- failure_summary(spec$method, "predictive selection", paste0("alpha=", spec$alpha), NA_real_, fit)
    } else {
      row <- evaluate_predictive_summary(spec$method, dat, reference_pred, fit$pred, fit$selected_count, fit$runtime_sec, fit$tuning)
    }
    path_rid <- path_rid + 1L
    path_rows[[path_rid]] <- align_to_template(row, template)
  }

  fit <- tryCatch(fit_group_lasso_predictive_benchmark(dat, seed = 20260529L + meta$replication_id + 333L), error = function(e) e)
  if (inherits(fit, "error")) {
    row <- failure_summary("Group lasso", "predictive selection", "group lasso", NA_real_, fit)
  } else {
    row <- evaluate_predictive_summary("Group lasso", dat, reference_pred, fit$pred, fit$selected_count, fit$runtime_sec, fit$tuning)
  }
  path_rid <- path_rid + 1L
  path_rows[[path_rid]] <- align_to_template(row, template)

  rid <- rid + 1L
  new_rows[[rid]] <- do.call(rbind, path_rows)
  write.csv(do.call(rbind, new_rows), checkpoint_path, row.names = FALSE)
  message("augmented ", basename(path))
}

augmented <- rbind(base_detail, do.call(rbind, new_rows))
out_path <- file.path(table_dir, sprintf("table_large_end_to_end_askpc_detail_%s.csv", out_tag))
write.csv(augmented, out_path, row.names = FALSE)
message("Wrote augmented large benchmark detail: ", out_path)
