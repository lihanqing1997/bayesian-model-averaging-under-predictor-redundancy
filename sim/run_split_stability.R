source(file.path("sim", "src", "support_kernel_benchmark.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  pat <- paste0("^--", name, "=")
  hit <- grep(pat, args, value = TRUE)
  if (length(hit)) sub(pat, "", hit[[1]]) else default
}

mode <- get_arg("mode", "smoke")
if (!mode %in% c("smoke", "screen", "full")) stop("--mode must be smoke, screen, or full")
run_tag <- get_arg("run-tag", "reliable_full_v2")
splits <- as.integer(get_arg("splits", switch(mode, smoke = "5", screen = "10", full = "20")))
max_cells <- as.integer(get_arg("max-cells", switch(mode, smoke = "1", screen = "4", full = "0")))
draws_per_split <- as.integer(get_arg("draws-per-split", switch(mode, smoke = "1000", screen = "3000", full = "8000")))
beta <- as.numeric(get_arg("beta", "0.02"))
max_iter <- as.integer(get_arg("max-iter", switch(mode, smoke = "40", screen = "80", full = "140")))

table_dir <- file.path("sim", "output", "tables")
object_dir <- file.path("sim", "output", "large_end_to_end")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

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

tex_escape <- function(x) {
  x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x
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

empirical_reference_from_draws <- function(ref, n, seed) {
  set.seed(seed)
  if (!is.null(ref$samples) && nrow(ref$samples) > 0L) {
    idx <- sample.int(nrow(ref$samples), size = min(n, nrow(ref$samples)), replace = n > nrow(ref$samples))
    draws <- as.matrix(ref$samples[idx, , drop = FALSE])
  } else {
    idx <- sample.int(nrow(ref$supports), size = n, replace = TRUE, prob = normalize_weights(ref$posterior))
    draws <- as.matrix(ref$supports[idx, , drop = FALSE])
  }
  keys <- apply(draws, 1, paste, collapse = "")
  tab <- table(keys)
  first <- match(names(tab), keys)
  list(
    supports = draws[first, , drop = FALSE],
    posterior = as.numeric(tab) / sum(tab)
  )
}

kernel_group_profile <- function(kernels, q, group_id, eta = 1e-3) {
  K <- max(group_id)
  prof <- rep(0, K)
  active <- which(q > eta)
  active <- active[!vapply(kernels[active], function(k) identical(k$type, "safety"), logical(1))]
  for (j in active) {
    k <- kernels[[j]]
    p <- k$params
    groups <- integer(0)
    if (k$type %in% c("hard_active", "soft_capacity")) {
      groups <- which(p$capacity > 0)
    } else if (k$type %in% c("soft_group_hamming", "hard_group_hamming_ball")) {
      groups <- p$center_groups
    } else if (k$type %in% c("soft_hamming", "hard_hamming_ball")) {
      groups <- unique(group_id[which(p$center_support > 0)])
    } else if (k$type %in% c("soft_interval", "soft_graph")) {
      nodes <- if (k$type == "soft_interval") p$interval else p$nodes
      groups <- unique(group_id[nodes])
    }
    if (length(groups)) prof[groups] <- prof[groups] + q[j]
  }
  prof
}

jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  u <- union(a, b)
  if (!length(u)) return(1)
  length(intersect(a, b)) / length(u)
}

pairwise_mean <- function(items, fun) {
  if (length(items) < 2L) return(NA_real_)
  vals <- c()
  for (i in seq_len(length(items) - 1L)) {
    for (j in (i + 1L):length(items)) vals <- c(vals, fun(items[[i]], items[[j]]))
  }
  mean(vals, na.rm = TRUE)
}

reference_paths <- sort(list.files(object_dir, pattern = paste0("_", run_tag, "\\.rds$"), full.names = TRUE))
if (!length(reference_paths)) stop("No reference objects found for run tag: ", run_tag)
if (max_cells > 0L && length(reference_paths) > max_cells) {
  reference_paths <- reference_paths[unique(round(seq(1, length(reference_paths), length.out = max_cells)))]
}

detail_rows <- list()
rid <- 0L
for (path in reference_paths) {
  meta <- parse_reference_path(path)
  obj <- readRDS(path)
  ref_full <- obj$reference_fit
  dat <- obj$data
  for (s in seq_len(splits)) {
    seed <- 20260612L + 1000L * s + meta$replication_id
    emp <- empirical_reference_from_draws(ref_full, draws_per_split, seed)
    fit <- tryCatch(
      fit_askpc_pooled_pruned(
        supports = emp$supports,
        weights = emp$posterior,
        group_id = dat$group_id,
        X = dat$X_train,
        beta = beta,
        tau = 1e-3,
        q0_min = 1e-3,
        prune_mass = 0.99,
        mode = if (mode == "full") "full" else "medium",
        include_metric_balls = TRUE,
        max_iter = max_iter,
        polish_pooled = mode != "smoke"
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      rid <- rid + 1L
      detail_rows[[rid]] <- data.frame(
        scenario = meta$scenario, rho = meta$rho, replication_id = meta$replication_id,
        split_id = s, status = "failure", error = conditionMessage(fit),
        tv = NA_real_, fkl = NA_real_, rkl = NA_real_, q0 = NA_real_, C_exp = NA_real_,
        C_list = NA_real_, K_list = NA_real_, K_eff = NA_real_, active_keys = NA_character_,
        active_families = NA_character_, coverage_profile = NA_character_, stringsAsFactors = FALSE
      )
      next
    }
    dict <- fit$pruned$dict
    W_full <- support_kernel_weight_matrix(ref_full$supports, dict, group_id = dat$group_id)
    alpha_full <- estimate_kernel_alpha(W_full, ref_full$posterior, alpha_floor = 1e-8)$alpha_truncated
    full_fit <- evaluate_family_mixture_fit(
      membership = W_full,
      alpha = alpha_full,
      costs = fit$pruned$costs,
      base_weights = ref_full$posterior,
      q = fit$pruned$fit$q,
      beta = beta,
      tau = 1e-3,
      distortion = "fkl",
      safety_index = 1L
    )
    q <- fit$pruned$fit$q
    active <- which(q > 1e-3)
    active <- active[!vapply(dict$kernels[active], function(k) identical(k$type, "safety"), logical(1))]
    q_ns <- q[-1]
    K_eff <- if (sum(q_ns) > 0) 1 / sum((q_ns / sum(q_ns))^2) else 0
    prof <- kernel_group_profile(dict$kernels, q, dat$group_id)
    rid <- rid + 1L
    detail_rows[[rid]] <- data.frame(
      scenario = meta$scenario,
      rho = meta$rho,
      replication_id = meta$replication_id,
      split_id = s,
      status = "ok",
      error = NA_character_,
      tv = full_fit$distortions$tv,
      fkl = full_fit$distortions$kl_base_to_compressed,
      rkl = full_fit$distortions$kl_compressed_to_base,
      q0 = q[1],
      C_exp = sum(q * fit$pruned$costs),
      C_list = if (length(active)) sum(fit$pruned$costs[active]) else 0,
      K_list = length(active),
      K_eff = K_eff,
      active_keys = paste(vapply(dict$kernels[active], function(k) k$key, character(1)), collapse = "||"),
      active_families = paste(vapply(dict$kernels[active], kernel_origin, character(1)), collapse = "||"),
      coverage_profile = paste(signif(prof, 6), collapse = ","),
      stringsAsFactors = FALSE
    )
  }
}

detail <- do.call(rbind, detail_rows)
detail_path <- file.path(table_dir, sprintf("table_split_stability_detail_%s.csv", mode))
write.csv(detail, detail_path, row.names = FALSE)

coverage_vec <- function(x) as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
summary_rows <- lapply(split(detail[detail$status == "ok", , drop = FALSE], interaction(detail$scenario[detail$status == "ok"], detail$rho[detail$status == "ok"], detail$replication_id[detail$status == "ok"], drop = TRUE)), function(d) {
  keys <- strsplit(d$active_keys, "||", fixed = TRUE)
  fams <- strsplit(d$active_families, "||", fixed = TRUE)
  covs <- lapply(d$coverage_profile, coverage_vec)
  cov_cor <- pairwise_mean(covs, function(a, b) {
    if (stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
    stats::cor(a, b)
  })
  data.frame(
    setting = paste0(d$scenario[1], ", rho=", d$rho[1], ", rep=", d$replication_id[1]),
    splits = nrow(d),
    tv_mean = mean(d$tv), tv_se = se(d$tv),
    fkl_mean = mean(d$fkl), fkl_se = se(d$fkl),
    q0_mean = mean(d$q0), q0_se = se(d$q0),
    C_list_mean = mean(d$C_list), C_list_se = se(d$C_list),
    active_jaccard = pairwise_mean(keys, jaccard),
    family_jaccard = pairwise_mean(fams, jaccard),
    coverage_correlation = cov_cor,
    stable = ifelse(is.finite(cov_cor) && cov_cor >= 0.75 && pairwise_mean(fams, jaccard) >= 0.5, "yes", "warning"),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
summary_path <- file.path(table_dir, sprintf("table_split_stability_summary_%s.csv", mode))
write.csv(summary, summary_path, row.names = FALSE)

tex <- c(
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  "Setting & Splits & TV & FKL & $q_0$ & Active Jaccard & Coverage cor. \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(summary))) {
  z <- summary[i, ]
  tex <- c(tex, paste0(
    tex_escape(z$setting), " & ", z$splits, " & ",
    fmt_pm(z$tv_mean, z$tv_se), " & ",
    fmt_pm(z$fkl_mean, z$fkl_se), " & ",
    fmt_pm(z$q0_mean, z$q0_se), " & ",
    fmt(z$active_jaccard), " & ",
    fmt(z$coverage_correlation), " \\\\"
  ))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(table_dir, sprintf("table_split_stability_%s.tex", mode)))
message("Wrote split-stability diagnostics: ", detail_path)
