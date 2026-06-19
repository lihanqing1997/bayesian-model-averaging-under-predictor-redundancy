source(file.path("sim", "src", "support_kernel_benchmark.R"))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit)) sub(pat, "", hit[[1]]) else default
}

mode <- get_arg("mode", "full")
if (!mode %in% c("smoke", "full")) stop("--mode must be smoke or full")

qmax_grid <- as.numeric(strsplit(get_arg("qmax", "0.05,0.10,0.25,0.50"), ",", fixed = TRUE)[[1]])
beta_grid <- as.numeric(strsplit(get_arg("beta", "0.02"), ",", fixed = TRUE)[[1]])
eta <- as.numeric(get_arg("eta", "0.001"))
alpha_floor <- as.numeric(get_arg("alpha-floor", "1e-8"))
max_cells <- as.integer(get_arg("max-cells", if (mode == "smoke") "2" else "0"))
max_iter <- as.integer(get_arg("max-iter", if (mode == "smoke") "30" else "80"))

table_dir <- file.path("sim", "output", "tables")
figure_dir <- file.path("sim", "output", "figures")
result_dir <- file.path("sim", "output", "results")
object_dir <- file.path("sim", "output", "large_end_to_end")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

fmt <- function(x, digits = 3) {
  if (!is.finite(x)) return("--")
  formatC(x, format = "f", digits = digits)
}

fmt_pm <- function(mu, err, digits = 3) {
  if (!is.finite(mu)) return("--")
  if (!is.finite(err)) return(fmt(mu, digits))
  paste0(fmt(mu, digits), " (", fmt(err, digits), ")")
}

cell_from_path <- function(path) {
  base <- basename(path)
  m <- regexec("^full_(.*)_rho([0-9]p[0-9]+)_rep([0-9]+)_reference_reliable_full_v2\\.rds$", base)
  z <- regmatches(base, m)[[1]]
  if (length(z) != 4L) return(NULL)
  list(
    scenario = z[2],
    rho = as.numeric(gsub("p", ".", z[3], fixed = TRUE)),
    replication_id = as.integer(z[4])
  )
}

reference_files <- list.files(
  object_dir,
  pattern = "^full_.*_reference_reliable_full_v2\\.rds$",
  full.names = TRUE
)
reference_files <- reference_files[vapply(reference_files, function(x) !is.null(cell_from_path(x)), logical(1))]
reference_files <- sort(reference_files)
if (max_cells > 0L) reference_files <- head(reference_files, max_cells)
if (!length(reference_files)) stop("No large-p reliable reference objects found")

report_metrics <- function(q, costs, eta = 1e-3, safety_index = 1L) {
  q <- normalize_weights(q)
  nonsafety <- setdiff(seq_along(q), safety_index)
  active <- nonsafety[q[nonsafety] > eta]
  q_ns_sum <- sum(q[nonsafety])
  if (q_ns_sum > 0) {
    qtilde <- q[nonsafety] / q_ns_sum
    keff <- 1 / sum(qtilde^2)
  } else {
    keff <- 0
  }
  list(
    cexp = sum(q * costs),
    clist = if (length(active)) sum(costs[active]) else 0,
    klist = length(active),
    keff = keff,
    active_nonfallback = length(active)
  )
}

floor_metrics <- function(q, W, alpha_used, weights, alpha_floor, dict = NULL,
                          sample_supports = NULL, group_id = NULL,
                          eta = 1e-3, safety_index = 1L) {
  q <- normalize_weights(q)
  weights <- normalize_weights(weights)
  alpha_emp <- as.numeric(crossprod(weights, W))
  floor_hit <- alpha_emp < alpha_floor
  active <- q > eta
  active[safety_index] <- FALSE
  delta <- sum(q[floor_hit] * pmax(0, 1 - alpha_emp[floor_hit] / alpha_floor))
  alpha_mcse <- rep(NA_real_, length(alpha_emp))
  if (!is.null(dict) && !is.null(sample_supports)) {
    Ws <- support_kernel_weight_matrix(as.matrix(sample_supports), dict, group_id = group_id)
    alpha_draw <- colMeans(Ws)
    alpha_mcse <- sqrt(pmax(alpha_draw * (1 - alpha_draw), 0) / nrow(Ws))
  }
  list(
    min_alpha_empirical = min(alpha_emp),
    floor_q_mass = sum(q[floor_hit]),
    floor_active_q_mass = sum(q[floor_hit & active]),
    delta_a = delta,
    active_floored_kernels = sum(floor_hit & active),
    max_active_reciprocal_uncertainty = suppressWarnings(max(alpha_mcse[active] / alpha_floor^2, na.rm = TRUE)),
    final_mass_mode = if (any(floor_hit & active)) "floored" else "unfloored empirical"
  )
}

fit_lean_pooled_pruned <- function(supports, weights, group_id, X,
                                   beta = 0.02, q0_max = 1,
                                   alpha_floor = 1e-8,
                                   prune_mass = 0.99,
                                   max_iter = 80L,
                                   eta = 1e-3) {
  supports <- as.matrix(supports)
  weights <- normalize_weights(weights)
  p <- ncol(supports)
  kernels <- list(safety_kernel(cost = max(10, 4 * log(max(p, 2)))))
  kernels <- c(
    kernels,
    generate_active_set_kernel_pool(
      supports = supports,
      group_id = group_id,
      weights = weights,
      max_sets = 12L,
      capacities = c(1L, 2L),
      rho_grid = c(1)
    ),
    generate_posterior_cluster_kernels(
      supports = supports,
      weights = weights,
      group_id = group_id,
      top_medoids = 10L,
      rho_grid = c(0.5, 1),
      group_level = TRUE
    ),
    generate_interval_kernel_pool(
      p = p,
      lengths = unique(pmin(c(3, 6, 12), p)),
      rho_grid = c(1),
      capacities = c(1L, 2L)
    )
  )
  seed_dict <- new_support_kernel_dictionary(dedupe_kernels(kernels))
  W_seed <- support_kernel_weight_matrix(supports, seed_dict, group_id = group_id)
  alpha_seed <- estimate_kernel_alpha(W_seed, weights, alpha_floor = alpha_floor)$alpha_truncated
  q_seed_fit <- optimize_family_mixture(
    W_seed,
    alpha_seed,
    kernel_costs(seed_dict),
    weights,
    beta = beta,
    tau = 1e-3,
    distortion = "fkl",
    q0_min = 1e-3,
    q0_max = q0_max,
    max_iter = max_iter,
    tol = 1e-7,
    polish = TRUE,
    polish_maxit = 400L
  )
  h_seed <- pmax(mixture_h(q_seed_fit$q, W_seed, alpha_seed), max(1e-3, 1e-10))
  residual_pool <- generate_residual_cover_kernels(
    supports = supports,
    residual_weights = weights / h_seed,
    group_id = group_id,
    top_medoids = 8L,
    rho_grid = c(0.5, 1),
    group_level = TRUE
  )
  kernels <- dedupe_kernels(c(seed_dict$kernels, residual_pool))
  dict <- new_support_kernel_dictionary(kernels)
  W <- support_kernel_weight_matrix(supports, dict, group_id = group_id)
  keep <- near_duplicate_kernel_keep(W, weights, alpha_floor = alpha_floor, tol = 1e-10, always_keep = 1L)
  kernels <- kernels[keep]
  dict <- new_support_kernel_dictionary(kernels)
  W <- support_kernel_weight_matrix(supports, dict, group_id = group_id)
  alpha <- estimate_kernel_alpha(W, weights, alpha_floor = alpha_floor)$alpha_truncated
  costs <- kernel_costs(dict)
  pooled_fit <- optimize_family_mixture(
    W,
    alpha,
    costs,
    weights,
    beta = beta,
    tau = 1e-3,
    distortion = "fkl",
    q0_min = 1e-3,
    q0_max = q0_max,
    max_iter = max_iter,
    tol = 1e-7,
    polish = TRUE,
    polish_maxit = 400L
  )
  keep_idx <- prune_kernel_indices_by_q_mass(pooled_fit$q, safety_index = 1L, mass = prune_mass)
  pruned_kernels <- dict$kernels[keep_idx]
  pruned_dict <- new_support_kernel_dictionary(pruned_kernels)
  Wp <- support_kernel_weight_matrix(supports, pruned_dict, group_id = group_id)
  alphap <- estimate_kernel_alpha(Wp, weights, alpha_floor = alpha_floor)$alpha_truncated
  costsp <- kernel_costs(pruned_dict)
  q_init <- pooled_fit$q[keep_idx]
  q_init <- q_init / sum(q_init)
  pruned_fit <- optimize_family_mixture(
    Wp,
    alphap,
    costsp,
    weights,
    beta = beta,
    tau = 1e-3,
    distortion = "fkl",
    q_init = q_init,
    q0_min = 1e-3,
    q0_max = q0_max,
    max_iter = max_iter,
    tol = 1e-7,
    polish = TRUE,
    polish_maxit = 400L
  )
  list(
    pooled = list(dict = dict, W = W, alpha = alpha, costs = costs, fit = pooled_fit),
    pruned = list(dict = pruned_dict, W = Wp, alpha = alphap, costs = costsp, fit = pruned_fit),
    raw_candidates = length(kernels),
    retained_q_mass = sum(pooled_fit$q[keep_idx])
  )
}

unfloored_distortions <- function(q, W, weights, q0_min) {
  alpha_emp <- as.numeric(crossprod(normalize_weights(weights), W))
  active <- q > 1e-12
  if (any(alpha_emp[active] <= 0)) {
    return(list(tv = NA_real_, fkl = NA_real_, rkl = NA_real_, mean_h_error = NA_real_))
  }
  d <- mixture_distortions(q, W, alpha_emp, weights, h_floor = max(q0_min, 1e-12))
  list(
    tv = d$tv,
    fkl = d$kl_base_to_compressed,
    rkl = d$kl_compressed_to_base,
    mean_h_error = abs(d$mean_h - 1)
  )
}

make_row <- function(cell, method, method_family, tuning, beta, qmax, stage, fit_obj,
                     W, alpha, costs, dict = NULL, ref, dat, runtime_sec,
                     eta = 1e-3, alpha_floor = 1e-8) {
  q <- fit_obj$fit$q
  rep <- report_metrics(q, costs, eta = eta)
  fl <- floor_metrics(
    q = q,
    W = W,
    alpha_used = alpha,
    weights = ref$posterior,
    alpha_floor = alpha_floor,
    dict = dict,
    sample_supports = ref$samples,
    group_id = dat$group_id,
    eta = eta
  )
  unf <- unfloored_distortions(q, W, ref$posterior, q0_min = 1e-3)
  data.frame(
    scenario = cell$scenario,
    rho = cell$rho,
    replication_id = cell$replication_id,
    method = method,
    method_family = method_family,
    tuning = tuning,
    beta = beta,
    tau = 1e-3,
    epsilon0 = 1e-3,
    qmax = qmax,
    q0_constraint_active = abs(q[1] - qmax) <= 5e-4,
    stage = stage,
    tv = unf$tv,
    fkl = unf$fkl,
    rkl = unf$rkl,
    mean_h_error = unf$mean_h_error,
    q0 = q[1],
    cexp = rep$cexp,
    clist = rep$clist,
    klist = rep$klist,
    keff = rep$keff,
    active_nonfallback = rep$active_nonfallback,
    n_kernels = ncol(W),
    objective = fit_obj$fit$objective,
    kkt_residual = fit_obj$fit$kkt_residual,
    optimizer_status = fit_obj$fit$status,
    optimizer_iterations = fit_obj$fit$iterations,
    objective_change = fit_obj$fit$objective_change,
    min_alpha_empirical = fl$min_alpha_empirical,
    floor_q_mass = fl$floor_q_mass,
    floor_active_q_mass = fl$floor_active_q_mass,
    delta_a = fl$delta_a,
    active_floored_kernels = fl$active_floored_kernels,
    max_active_reciprocal_uncertainty = fl$max_active_reciprocal_uncertainty,
    final_mass_mode = fl$final_mass_mode,
    runtime_sec = runtime_sec,
    stringsAsFactors = FALSE
  )
}

hard_ball_matrix <- function(supports, centers, radii, group_id = NULL, group_level = FALSE) {
  supports <- as.matrix(supports)
  if (group_level) {
    S <- support_group_counts(supports, group_id) > 0
    C <- support_group_counts(centers, group_id) > 0
  } else {
    S <- supports > 0
    C <- centers > 0
  }
  W <- matrix(0, nrow(S), nrow(C) * length(radii))
  names <- character(ncol(W))
  k <- 0L
  for (i in seq_len(nrow(C))) {
    dist <- rowSums(sweep(S, 2, C[i, ], `!=`))
    for (r in radii) {
      k <- k + 1L
      W[, k] <- as.numeric(dist <= r)
      names[k] <- paste0(if (group_level) "group" else "support", "_ball_", i, "_r", r)
    }
  }
  colnames(W) <- names
  W
}

fit_ball_baseline <- function(ref, dat, beta, group_level = FALSE,
                              n_centers = 5L, radii = 0:5, q0_max = 1) {
  weights <- normalize_weights(ref$posterior)
  ord <- order(weights, decreasing = TRUE)
  centers <- unique(ref$supports[ord[seq_len(min(length(ord), n_centers * 4L))], , drop = FALSE])
  centers <- centers[seq_len(min(nrow(centers), n_centers)), , drop = FALSE]
  W_ball <- hard_ball_matrix(ref$supports, centers, radii, group_id = dat$group_id, group_level = group_level)
  W <- cbind(1, W_ball)
  alpha <- estimate_family_alpha(W, weights, alpha_floor = alpha_floor)$alpha_truncated
  costs <- c(max(10, 4 * log(max(ncol(ref$supports), 2))), rep(log(ncol(ref$supports) + 1), ncol(W_ball)) + rep(radii, each = nrow(centers)))
  fit <- optimize_family_mixture(
    W,
    alpha,
    costs,
    weights,
    beta = beta,
    tau = 1e-3,
    distortion = "fkl",
    q0_min = 1e-3,
    q0_max = q0_max,
    max_iter = 120L,
    tol = 1e-7
  )
  list(W = W, alpha = alpha, costs = costs, fit = fit, dict = NULL)
}

q0_rows <- list()
baseline_rows <- list()
start_time <- Sys.time()

for (path in reference_files) {
  cell <- cell_from_path(path)
  message(sprintf("Cell %s rho %.2f rep %03d", cell$scenario, cell$rho, cell$replication_id))
  obj <- readRDS(path)
  dat <- obj$data
  ref <- obj$reference_fit
  for (beta in beta_grid) {
    for (qmax in qmax_grid) {
      label <- sprintf("beta=%.3g, qmax=%.2f", beta, qmax)
      t0 <- proc.time()[["elapsed"]]
      ask <- fit_lean_pooled_pruned(
        supports = ref$supports,
        weights = ref$posterior,
        group_id = dat$group_id,
        X = dat$X_train,
        beta = beta,
        alpha_floor = alpha_floor,
        q0_max = qmax,
        prune_mass = 0.99,
        max_iter = max_iter,
        eta = eta
      )
      rt <- proc.time()[["elapsed"]] - t0
      q0_rows[[length(q0_rows) + 1L]] <- make_row(cell, "Pooled union", "support-kernel compression", label, beta, qmax, "pre-prune", ask$pooled, ask$pooled$W, ask$pooled$alpha, ask$pooled$costs, ask$pooled$dict, ref, dat, rt, eta, alpha_floor)
      q0_rows[[length(q0_rows) + 1L]] <- make_row(cell, "Pooled-pruned", "support-kernel compression", label, beta, qmax, "post-prune", ask$pruned, ask$pruned$W, ask$pruned$alpha, ask$pruned$costs, ask$pruned$dict, ref, dat, rt, eta, alpha_floor)
    }
  }

  beta0 <- beta_grid[1]
  t0 <- proc.time()[["elapsed"]]
  hb <- fit_ball_baseline(ref, dat, beta = beta0, group_level = FALSE, n_centers = 5L, radii = 0:5)
  baseline_rows[[length(baseline_rows) + 1L]] <- make_row(cell, "Hamming balls", "same-target credible region", "top centers, r=0:5", beta0, 1, "fit", hb, hb$W, hb$alpha, hb$costs, NULL, ref, dat, proc.time()[["elapsed"]] - t0, eta, alpha_floor)

  t0 <- proc.time()[["elapsed"]]
  ghb <- fit_ball_baseline(ref, dat, beta = beta0, group_level = TRUE, n_centers = 5L, radii = 0:3)
  baseline_rows[[length(baseline_rows) + 1L]] <- make_row(cell, "Group-Hamming balls", "same-target credible region", "top group centers, r=0:3", beta0, 1, "fit", ghb, ghb$W, ghb$alpha, ghb$costs, NULL, ref, dat, proc.time()[["elapsed"]] - t0, eta, alpha_floor)
}

q0_detail <- do.call(rbind, q0_rows)
baseline_detail <- do.call(rbind, baseline_rows)

write.csv(q0_detail, file.path(result_dir, "q0_constrained_largep.csv"), row.names = FALSE)
write.csv(q0_detail, file.path(result_dir, "floor_diagnostics_largep.csv"), row.names = FALSE)
write.csv(q0_detail, file.path(result_dir, "report_length_diagnostics.csv"), row.names = FALSE)
write.csv(baseline_detail, file.path(result_dir, "same_target_baselines_largep.csv"), row.names = FALSE)

summarize_by <- function(d, group_cols) {
  split_d <- split(d, interaction(d[group_cols], drop = TRUE), drop = TRUE)
  out <- lapply(split_d, function(z) {
    head <- z[1, group_cols, drop = FALSE]
    nums <- c("tv", "fkl", "rkl", "q0", "cexp", "clist", "klist", "keff", "active_nonfallback",
              "objective", "kkt_residual", "floor_q_mass", "floor_active_q_mass", "delta_a",
              "active_floored_kernels", "runtime_sec")
    for (nm in nums) {
      head[[paste0(nm, "_mean")]] <- mean(z[[nm]], na.rm = TRUE)
      head[[paste0(nm, "_se")]] <- se(z[[nm]])
    }
    head$n_cells <- length(unique(paste(z$scenario, z$rho, z$replication_id, sep = "::")))
    head
  })
  do.call(rbind, out)
}

q0_summary <- summarize_by(q0_detail[q0_detail$stage == "post-prune", , drop = FALSE], c("method", "qmax", "beta", "stage"))
write.csv(q0_summary, file.path(result_dir, "q0_constrained_largep_summary.csv"), row.names = FALSE)

floor_summary <- summarize_by(q0_detail[q0_detail$stage == "post-prune", , drop = FALSE], c("method", "qmax", "beta", "stage"))
write.csv(floor_summary, file.path(result_dir, "floor_diagnostics_largep_summary.csv"), row.names = FALSE)

baseline_summary <- summarize_by(baseline_detail, c("method", "beta", "stage"))
write.csv(baseline_summary, file.path(result_dir, "same_target_baselines_largep_summary.csv"), row.names = FALSE)

tex_q0 <- c(
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  "$q_0$ bound & TV & FKL & $q_0$ & $C_{\\mathrm{exp}}$ & $C_{\\mathrm{list}}$ & KKT \\\\",
  "\\midrule"
)
q0_show <- q0_summary[order(q0_summary$qmax), , drop = FALSE]
for (i in seq_len(nrow(q0_show))) {
  tex_q0 <- c(tex_q0, sprintf(
    "$q_0\\le%s$ & %s & %s & %s & %s & %s & %s \\\\",
    fmt(q0_show$qmax[i], 2),
    fmt_pm(q0_show$tv_mean[i], q0_show$tv_se[i]),
    fmt_pm(q0_show$fkl_mean[i], q0_show$fkl_se[i]),
    fmt_pm(q0_show$q0_mean[i], q0_show$q0_se[i]),
    fmt_pm(q0_show$cexp_mean[i], q0_show$cexp_se[i], 2),
    fmt_pm(q0_show$clist_mean[i], q0_show$clist_se[i], 1),
    fmt_pm(q0_show$kkt_residual_mean[i], q0_show$kkt_residual_se[i], 3)
  ))
}
tex_q0 <- c(tex_q0, "\\bottomrule", "\\end{tabular}")
writeLines(tex_q0, file.path(table_dir, "table_q0_constrained_largep.tex"))

tex_floor <- c(
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "$q_0$ bound & Floor $q$-mass & Active floor $q$ & $\\Delta_a(q)$ & Active floored & Mass mode \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(q0_show))) {
  mode_tab <- if (q0_show$floor_active_q_mass_mean[i] > 1e-8) "floored" else "unfloored"
  tex_floor <- c(tex_floor, sprintf(
    "$q_0\\le%s$ & %s & %s & %s & %s & %s \\\\",
    fmt(q0_show$qmax[i], 2),
    fmt_pm(q0_show$floor_q_mass_mean[i], q0_show$floor_q_mass_se[i], 4),
    fmt_pm(q0_show$floor_active_q_mass_mean[i], q0_show$floor_active_q_mass_se[i], 4),
    fmt_pm(q0_show$delta_a_mean[i], q0_show$delta_a_se[i], 4),
    fmt_pm(q0_show$active_floored_kernels_mean[i], q0_show$active_floored_kernels_se[i], 1),
    mode_tab
  ))
}
tex_floor <- c(tex_floor, "\\bottomrule", "\\end{tabular}")
writeLines(tex_floor, file.path(table_dir, "table_floor_diagnostics_largep.tex"))

tex_base <- c(
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  "Method & TV & FKL & $q_0$ & $C_{\\mathrm{exp}}$ & $C_{\\mathrm{list}}$ & $K_{\\mathrm{list}}$ \\\\",
  "\\midrule"
)
baseline_show <- baseline_summary[order(baseline_summary$method), , drop = FALSE]
for (i in seq_len(nrow(baseline_show))) {
  tex_base <- c(tex_base, sprintf(
    "%s & %s & %s & %s & %s & %s & %s \\\\",
    baseline_show$method[i],
    fmt_pm(baseline_show$tv_mean[i], baseline_show$tv_se[i]),
    fmt_pm(baseline_show$fkl_mean[i], baseline_show$fkl_se[i]),
    fmt_pm(baseline_show$q0_mean[i], baseline_show$q0_se[i]),
    fmt_pm(baseline_show$cexp_mean[i], baseline_show$cexp_se[i], 2),
    fmt_pm(baseline_show$clist_mean[i], baseline_show$clist_se[i], 1),
    fmt_pm(baseline_show$klist_mean[i], baseline_show$klist_se[i], 1)
  ))
}
tex_base <- c(tex_base, "\\bottomrule", "\\end{tabular}")
writeLines(tex_base, file.path(table_dir, "table_same_target_baselines.tex"))

report_summary <- rbind(
  q0_summary[, c("method", "qmax", "beta", "cexp_mean", "cexp_se", "clist_mean", "clist_se", "klist_mean", "klist_se", "keff_mean", "keff_se"), drop = FALSE],
  transform(baseline_summary[, c("method", "beta", "cexp_mean", "cexp_se", "clist_mean", "clist_se", "klist_mean", "klist_se", "keff_mean", "keff_se"), drop = FALSE], qmax = NA_real_)
)
write.csv(report_summary, file.path(result_dir, "report_length_diagnostics_summary.csv"), row.names = FALSE)

tex_report <- c(
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "Method & $q_0$ bound & $C_{\\mathrm{exp}}$ & $C_{\\mathrm{list}}$ & $K_{\\mathrm{list}}$ & $K_{\\mathrm{eff}}$ \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(report_summary))) {
  bound <- if (is.finite(report_summary$qmax[i])) paste0("$q_0\\le", fmt(report_summary$qmax[i], 2), "$") else "--"
  tex_report <- c(tex_report, sprintf(
    "%s & %s & %s & %s & %s & %s \\\\",
    report_summary$method[i],
    bound,
    fmt_pm(report_summary$cexp_mean[i], report_summary$cexp_se[i], 2),
    fmt_pm(report_summary$clist_mean[i], report_summary$clist_se[i], 1),
    fmt_pm(report_summary$klist_mean[i], report_summary$klist_se[i], 1),
    fmt_pm(report_summary$keff_mean[i], report_summary$keff_se[i], 1)
  ))
}
tex_report <- c(tex_report, "\\bottomrule", "\\end{tabular}")
writeLines(tex_report, file.path(table_dir, "table_report_length_diagnostics.tex"))

q0_plot <- q0_summary[q0_summary$stage == "post-prune", , drop = FALSE]
plot_data <- data.frame(
  method = sprintf("Pooled-pruned q0<=%.2f", q0_plot$qmax),
  family = "Pooled-pruned",
  fkl = q0_plot$fkl_mean,
  cexp = q0_plot$cexp_mean,
  clist = q0_plot$clist_mean,
  q0 = q0_plot$q0_mean,
  stringsAsFactors = FALSE
)
if (nrow(baseline_summary)) {
  plot_data <- rbind(
    plot_data,
    data.frame(
      method = baseline_summary$method,
      family = baseline_summary$method,
      fkl = baseline_summary$fkl_mean,
      cexp = baseline_summary$cexp_mean,
      clist = baseline_summary$clist_mean,
      q0 = baseline_summary$q0_mean,
      stringsAsFactors = FALSE
    )
  )
}
write.csv(plot_data, file.path(result_dir, "largep_frontiers_clean_data.csv"), row.names = FALSE)

draw_frontier <- function(device, file) {
  device(file, width = 10.2, height = 3.8)
  old <- par(mfrow = c(1, 3), mar = c(3.45, 4.15, 1.35, 0.65), oma = c(0, 0, 0, 0), las = 1, pty = "s")
  on.exit({ par(old); dev.off() }, add = TRUE)
  cols <- c("Pooled-pruned" = "#1b9e77",
            "Group-Hamming balls" = "#d95f02",
            "Hamming balls" = "#7570b3")
  pchs <- c("Pooled-pruned" = 17,
            "Group-Hamming balls" = 19,
            "Hamming balls" = 15)
  panel <- function(x, xlab, main, qline = FALSE, legend_panel = FALSE) {
    xpad <- diff(range(x, finite = TRUE)) * 0.08
    ypad <- diff(range(plot_data$fkl, finite = TRUE)) * 0.10
    if (!is.finite(xpad) || xpad == 0) xpad <- 1
    if (!is.finite(ypad) || ypad == 0) ypad <- 0.01
    plot(x, plot_data$fkl,
         pch = pchs[plot_data$family],
         col = cols[plot_data$family],
         xlab = xlab,
         ylab = "FKL",
         main = main,
         cex = 1.35,
         cex.lab = 1.04,
         cex.axis = 0.92,
         cex.main = 1.08,
         xlim = range(x, finite = TRUE) + c(-xpad, xpad),
         ylim = range(plot_data$fkl, finite = TRUE) + c(-ypad, ypad))
    grid(col = "gray90")
    if (qline) abline(v = 0.25, lty = 2, col = "gray45")
    if (legend_panel) {
      legend("topleft",
             inset = c(0.035, 0.025),
             legend = names(cols),
             col = cols,
             pch = pchs,
             bty = "n",
             cex = 0.64,
             pt.cex = 0.92)
    }
  }
  panel(plot_data$cexp, expression(C[exp]), "A", legend_panel = TRUE)
  panel(plot_data$clist, expression(C[list]), "B")
  panel(plot_data$q0, expression(q[0]), "C", qline = TRUE)
}
draw_frontier(grDevices::pdf, file.path(figure_dir, "largep_frontiers_clean.pdf"))
draw_frontier(function(file, width, height) grDevices::png(file, width = width, height = height, units = "in", res = 300),
              file.path(figure_dir, "largep_frontiers_clean.png"))

log_lines <- c(
  "# Computational Revision Log",
  "",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Mode: ", mode),
  paste0("Reference objects used: ", length(reference_files)),
  "",
  "## Existing scripts and objects",
  "- Manuscript source is maintained separately from the companion code release.",
  "- Core optimization code: sim/src/family_mixture_compression.R and sim/src/adaptive_support_kernel_compression.R.",
  "- Large-p benchmark driver: sim/run_large_end_to_end_askpc_benchmark.R.",
  "- Existing large-p reference posterior objects: sim/output/large_end_to_end/*_reference_reliable_full_v2.rds.",
  "- Existing detailed summaries: sim/output/tables/table_large_end_to_end_askpc_detail_reliable_full_v2_augmented.csv.",
  "",
  "## Generated in this run",
  "- results/q0_constrained_largep.csv from actual refits with epsilon0 <= q0 <= qmax.",
  "- The constrained mixture refits use BFGS polish after projected first-order optimization.",
  "- results/floor_diagnostics_largep.csv from fitted kernels and empirical retained masses.",
  "- results/report_length_diagnostics.csv with Cexp, Clist, Klist, and Keff.",
  "- results/same_target_baselines_largep.csv for Hamming and group-Hamming ball baselines.",
  "- tables/table_q0_constrained_largep.tex.",
  "- tables/table_floor_diagnostics_largep.tex.",
  "- tables/table_report_length_diagnostics.tex.",
  "- tables/table_same_target_baselines.tex.",
  "- figures/largep_frontiers_clean.pdf and figures/largep_frontiers_clean.png.",
  "",
  "## Feasibility notes",
  "- Large-p refits used saved reference posterior objects, not a new MC3 run.",
  "- Exact and Tecator constrained refits were not run by this script.",
  "- The new q0-constrained rows are actual constrained optimizer runs, not screens over a previous grid.",
  paste0("Elapsed seconds: ", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1))
)
writeLines(log_lines, file.path("sim", "output", "results", "largep_report_diagnostics_log.txt"))

cat("Large-p report diagnostics written.\n")
