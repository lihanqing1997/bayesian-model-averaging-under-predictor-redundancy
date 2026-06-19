table_dir <- file.path("sim", "output", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

escape_tex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x, fixed = TRUE)
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "--", sprintf(paste0("%.", digits, "f"), x))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "--", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

read_best_fkl <- function(path, id_cols) {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[d$status == "ok" & is.finite(d$fkl), , drop = FALSE]
  d$key <- do.call(paste, c(d[id_cols], sep = "|"))
  groups <- split(d, interaction(d$key, d$method, drop = TRUE))
  best <- do.call(rbind, lapply(groups, function(x) x[which.min(x$fkl), , drop = FALSE]))
  rownames(best) <- NULL
  best
}

paired_gap <- function(best, base, comparison, metric) {
  b <- best[best$method == base, c("key", metric), drop = FALSE]
  z <- best[best$method == comparison, c("key", metric), drop = FALSE]
  names(b)[2] <- "base"
  names(z)[2] <- "comparison"
  m <- merge(b, z, by = "key")
  gap <- m$comparison - m$base
  data.frame(
    n_pairs = length(gap),
    mean_gap = mean(gap),
    median_gap = median(gap),
    paired_t_p = tryCatch(t.test(gap)$p.value, error = function(e) NA_real_),
    wilcoxon_p = tryCatch(wilcox.test(gap, exact = FALSE)$p.value, error = function(e) NA_real_),
    stringsAsFactors = FALSE
  )
}

comparison_rows <- list()
add_comparison <- function(evidence, best, comparison, metric, label = NULL) {
  base <- "ASK-PC pooled-pruned 99%"
  if (!comparison %in% best$method) return(invisible(NULL))
  out <- paired_gap(best, base, comparison, metric)
  out$evidence <- evidence
  out$comparison <- if (is.null(label)) comparison else label
  out$metric <- metric
  comparison_rows[[length(comparison_rows) + 1L]] <<- out
}

exact_best <- read_best_fkl(
  file.path(table_dir, "table_support_kernel_competitor_benchmark_detail.csv"),
  c("scenario", "rho", "replication_id")
)
large_best <- read_best_fkl(
  file.path(table_dir, "table_large_end_to_end_askpc_detail_reliable_full_v2_augmented.csv"),
  c("scenario", "rho", "replication_id")
)
real_best <- read_best_fkl(
  file.path(table_dir, "table_real_response_reduced_detail.csv"),
  c("dataset", "replication_id")
)

for (best_name in c("Exact", "Large end-to-end", "Reduced real-response")) {
  best <- switch(
    best_name,
    "Exact" = exact_best,
    "Large end-to-end" = large_best,
    "Reduced real-response" = real_best
  )
  add_comparison(best_name, best, "Top-M support atoms", "fkl", "Top-M atoms")
  add_comparison(best_name, best, "Credible support set", "fkl", "Credible set")
  add_comparison(best_name, best, "Posterior clustering", "fkl", "Posterior clustering")
  add_comparison(best_name, best, "Posterior clustering", "expected_code", "Posterior clustering")
  add_comparison(best_name, best, "Dilution-prior BMA", "fkl", "Dilution-prior BMA")
  add_comparison(best_name, best, "DPP-prior BMA", "fkl", "DPP-prior BMA")
}

paired <- do.call(rbind, comparison_rows)
paired <- paired[, c(
  "evidence", "comparison", "metric", "n_pairs", "mean_gap",
  "median_gap", "paired_t_p", "wilcoxon_p"
)]
write.csv(paired, file.path(table_dir, "table_key_paired_comparisons.csv"), row.names = FALSE)

show_rows <- paired[
  (paired$evidence == "Exact" & paired$comparison %in% c("Top-M atoms", "Posterior clustering", "Dilution-prior BMA")) |
    (paired$evidence == "Large end-to-end" & paired$comparison %in% c("Top-M atoms", "Posterior clustering", "Dilution-prior BMA")) |
    (paired$evidence == "Reduced real-response" & paired$comparison %in% c("Top-M atoms", "Posterior clustering")),
  ,
  drop = FALSE
]

tex <- c(
  "\\begin{tabular}{llcrr}",
  "\\toprule",
  "Evidence & Comparison & metric & mean gap & Wilcoxon $p$ \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(show_rows))) {
  row <- show_rows[i, ]
  metric_label <- if (row$metric == "fkl") "FKL" else "reporting cost"
  tex <- c(
    tex,
    sprintf(
      "%s & %s & %s & %s & %s \\\\",
      escape_tex(row$evidence),
      escape_tex(row$comparison),
      metric_label,
      fmt_num(row$mean_gap, 3),
      fmt_p(row$wilcoxon_p)
    )
  )
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(table_dir, "table_key_paired_comparisons.tex"))

large_detail <- read.csv(
  file.path(table_dir, "table_large_end_to_end_askpc_detail_reliable_full_v2_augmented.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ref <- unique(large_detail[, c(
  "scenario", "rho", "replication_id", "n_iter", "burn", "thin",
  "n_chains", "retained_draws", "unique_supports", "ess_min",
  "split_rhat_max", "group_pip_mcse_max", "runtime_sec_reference",
  "diagnostic_status", "diagnostic_warning"
)])

ref_summary <- data.frame(
  quantity = c(
    "reference runs",
    "chains per run",
    "iterations per chain",
    "burn-in per chain",
    "thinning interval",
    "retained draws per run",
    "reference runtime, seconds",
    "unique supports",
    "minimum ESS",
    "maximum split Rhat",
    "maximum group-PIP MCSE",
    "passing core diagnostic gate",
    "low-acceptance warnings"
  ),
  value = c(
    nrow(ref),
    paste0(unique(ref$n_chains), collapse = ", "),
    paste0(unique(ref$n_iter), collapse = ", "),
    paste0(unique(ref$burn), collapse = ", "),
    paste0(unique(ref$thin), collapse = ", "),
    paste0(unique(ref$retained_draws), collapse = ", "),
    sprintf("%.1f [%.1f, %.1f]", mean(ref$runtime_sec_reference), min(ref$runtime_sec_reference), max(ref$runtime_sec_reference)),
    sprintf("%.1f [%.1f, %.1f]", mean(ref$unique_supports), min(ref$unique_supports), max(ref$unique_supports)),
    sprintf("%.1f [%.1f, %.1f]", mean(ref$ess_min), min(ref$ess_min), max(ref$ess_min)),
    sprintf("%.3f [%.3f, %.3f]", mean(ref$split_rhat_max), min(ref$split_rhat_max), max(ref$split_rhat_max)),
    sprintf("%.3f [%.3f, %.3f]", mean(ref$group_pip_mcse_max), min(ref$group_pip_mcse_max), max(ref$group_pip_mcse_max)),
    sprintf("%d/%d", sum(ref$diagnostic_status == "pass"), nrow(ref)),
    sum(nzchar(ref$diagnostic_warning))
  ),
  stringsAsFactors = FALSE
)
write.csv(ref_summary, file.path(table_dir, "table_reference_computation_cost.csv"), row.names = FALSE)

tex2 <- c(
  "\\begin{tabular}{ll}",
  "\\toprule",
  "Quantity & Value \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(ref_summary))) {
  tex2 <- c(
    tex2,
    sprintf("%s & %s \\\\", escape_tex(ref_summary$quantity[i]), escape_tex(ref_summary$value[i]))
  )
}
tex2 <- c(tex2, "\\bottomrule", "\\end{tabular}")
writeLines(tex2, file.path(table_dir, "table_reference_computation_cost.tex"))

message("Wrote paired comparison and reference-cost outputs.")
