table_dir <- file.path("sim", "output", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

fmt <- function(x, digits = 3) {
  y <- as.numeric(x)
  y[is.finite(y) & abs(y) < 0.5 * 10^(-digits)] <- 0
  ifelse(is.na(y), "--", formatC(y, format = "f", digits = digits))
}

label_method <- function(x) {
  map <- c(
    fixed_hard_dictionary = "Fixed hard",
    adaptive_support_kernel = "Adaptive kernel",
    topM_support_atoms = "Top-$M$ atoms",
    full_adaptive = "Full adaptive",
    no_residual_cover = "No residual",
    no_posterior_cluster = "No cluster",
    posterior_cluster_only = "Cluster only",
    residual_cover_only = "Residual only",
    intervals_only = "Intervals only",
    graph_only = "Graph only",
    pooled_candidate_union = "Pooled union",
    pooled_pruned_99 = "Pooled-pruned"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

write_lines <- function(lines, file) {
  con <- file(file, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con)
}

write_ablation_table <- function() {
  path <- file.path(table_dir, "table_adaptive_kernel_ablation.csv")
  if (!file.exists(path)) {
    return(FALSE)
  }
  ab <- read.csv(path)
  ab <- subset(ab, scenario == "one_representative" & rho == 0.95)
  method_order <- c(
    "pooled_candidate_union",
    "pooled_pruned_99",
    "full_adaptive",
    "no_residual_cover",
    "no_posterior_cluster",
    "posterior_cluster_only",
    "residual_cover_only",
    "intervals_only",
    "graph_only",
    "fixed_hard_dictionary",
    "topM_support_atoms"
  )
  ab$order <- match(ab$method, method_order)
  ab <- ab[order(ab$order), ]
  ab$label <- label_method(ab$method)

  lines <- c(
    "\\begin{tabular}{lrrrrrrrrrr}",
    "\\toprule",
    "Method & TV & FKL & Reporting cost & $\\beta c$ & $\\tau\\sum q\\log q$ & Obj. & $q_0$ & Eff. & Active & Sec. \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(ab))) {
    active_val <- if ("active_kernels_001_mean" %in% names(ab)) ab$active_kernels_001_mean[i] else NA_real_
    active <- if (is.finite(active_val)) fmt(active_val, 1) else "--"
    obj_val <- if ("objective_mean" %in% names(ab)) ab$objective_mean[i] else NA_real_
    code_penalty <- if ("code_penalty_mean" %in% names(ab)) ab$code_penalty_mean[i] else NA_real_
    entropy_penalty <- if ("entropy_penalty_mean" %in% names(ab)) ab$entropy_penalty_mean[i] else NA_real_
    obj_check <- if ("objective_check_mean" %in% names(ab)) ab$objective_check_mean[i] else obj_val
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
      ab$label[i],
      fmt(ab$tv_mean[i]),
      fmt(ab$fkl_mean[i]),
      fmt(ab$expected_code_mean[i]),
      fmt(code_penalty),
      fmt(entropy_penalty),
      fmt(obj_check),
      fmt(ab$q0_mean[i]),
      fmt(ab$q_effective_kernels_mean[i], 1),
      active,
      fmt(ab$runtime_sec_mean[i], 2)
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(table_dir, "table_adaptive_kernel_ablation.tex"))
  TRUE
}

write_semisynthetic_table <- function() {
  path <- file.path(table_dir, "table_adaptive_kernel_semisynthetic_realx.csv")
  if (!file.exists(path)) {
    return(FALSE)
  }
  realx <- read.csv(path)
  realx$label <- label_method(realx$method)
  lines <- c(
    "\\begin{tabular}{lrrrrrr}",
    "\\toprule",
    "Method & TV & FKL & $q_0$ & Reporting cost & RMSE gap & Kernels \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(realx))) {
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & %s & %s & %s \\\\",
      realx$label[i],
      fmt(realx$tv_mean[i]),
      fmt(realx$fkl_mean[i]),
      fmt(realx$q0_mean[i]),
      fmt(realx$expected_code_mean[i]),
      fmt(realx$rmse_gap_mean[i], 4),
      fmt(realx$n_kernels_mean[i], 1)
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  write_lines(lines, file.path(table_dir, "table_adaptive_kernel_semisynthetic_realx.tex"))
  TRUE
}

ok <- c(
  ablation = write_ablation_table(),
  semisynthetic = write_semisynthetic_table()
)

cat("Adaptive-kernel LaTeX tables written: ", paste(names(ok)[ok], collapse = ", "), "\n", sep = "")
