table_dir <- file.path("sim", "output", "tables")

summary_path <- file.path(table_dir, "table_large_end_to_end_askpc_summary.csv")
diag_path <- file.path(table_dir, "table_large_end_to_end_reference_diagnostics.csv")

if (!file.exists(summary_path)) {
  stop("Missing ", summary_path, ". Run sim/run_large_end_to_end_askpc_benchmark.R --mode=summarize first.")
}
if (!file.exists(diag_path)) {
  stop("Missing ", diag_path, ". Run sim/run_large_end_to_end_askpc_benchmark.R --mode=summarize first.")
}

summary <- read.csv(summary_path, stringsAsFactors = FALSE)
diag <- read.csv(diag_path, stringsAsFactors = FALSE)

method_keep <- c(
  "ASK-PC pooled-pruned 99%",
  "Posterior clustering",
  "Top-M support atoms",
  "Credible support set",
  "Dilution-prior BMA",
  "DPP-prior BMA",
  "Fixed hard dictionary"
)

pretty_method <- function(x) {
  out <- x
  out[out == "ASK-PC pooled-pruned 99%"] <- "Pooled-pruned"
  out[out == "Posterior clustering"] <- "Cluster kernels"
  out[out == "Fixed hard dictionary"] <- "Fixed regions"
  out[out == "Top-M support atoms"] <- "Top-M atoms"
  out[out == "Credible support set"] <- "Credible set"
  out[out == "Dilution-prior BMA"] <- "Dilution BMA"
  out[out == "DPP-prior BMA"] <- "DPP BMA"
  out
}

se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

fmt_num <- function(x, digits = 3) {
  if (!is.finite(x)) return("--")
  formatC(x, format = "f", digits = digits)
}

fmt_mean_se <- function(mu, se_val, digits = 3) {
  if (!is.finite(mu)) return("--")
  if (!is.finite(se_val)) return(fmt_num(mu, digits))
  paste0(fmt_num(mu, digits), " (", fmt_num(se_val, digits), ")")
}

best_one <- function(d) {
  d <- d[is.finite(d$fkl_mean), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d[which.min(d$fkl_mean), , drop = FALSE]
}

summary_sub <- summary[summary$method %in% method_keep, , drop = FALSE]
parts <- split(summary_sub, interaction(summary_sub$scenario, summary_sub$rho, summary_sub$method, drop = TRUE), drop = TRUE)
best <- do.call(rbind, Filter(Negate(is.null), lapply(parts, best_one)))
rownames(best) <- NULL

metric_cols <- c("tv_mean", "fkl_mean", "q0_mean", "expected_code_mean", "rmse_gap_mean")
agg_rows <- lapply(method_keep, function(m) {
  z <- best[best$method == m, , drop = FALSE]
  data.frame(
    method = pretty_method(m),
    n_cells = nrow(z),
    tv = mean(z$tv_mean, na.rm = TRUE),
    tv_se = se(z$tv_mean),
    fkl = mean(z$fkl_mean, na.rm = TRUE),
    fkl_se = se(z$fkl_mean),
    q0 = mean(z$q0_mean, na.rm = TRUE),
    q0_se = se(z$q0_mean),
    code = mean(z$expected_code_mean, na.rm = TRUE),
    code_se = se(z$expected_code_mean),
    list_count = mean(ifelse(is.finite(z$stored_atoms_mean), z$stored_atoms_mean, z$active_kernels_001_mean), na.rm = TRUE),
    list_count_se = se(ifelse(is.finite(z$stored_atoms_mean), z$stored_atoms_mean, z$active_kernels_001_mean)),
    rmse_gap = mean(z$rmse_gap_mean, na.rm = TRUE),
    rmse_gap_se = se(z$rmse_gap_mean),
    stringsAsFactors = FALSE
  )
})
agg <- do.call(rbind, agg_rows)
agg$q0[is.nan(agg$q0)] <- NA_real_
agg$q0_se[is.nan(agg$q0_se)] <- NA_real_

write.csv(agg, file.path(table_dir, "table_large_end_to_end_askpc_family_summary.csv"), row.names = FALSE)

predictive <- summary[summary$method_family == "predictive selection", , drop = FALSE]
if (nrow(predictive)) {
  pred_parts <- split(predictive, interaction(predictive$scenario, predictive$method, drop = TRUE), drop = TRUE)
  pred_best <- do.call(rbind, lapply(pred_parts, function(d) d[which.min(d$rmse_gap_mean), , drop = FALSE]))
  pred_agg <- do.call(rbind, lapply(unique(pred_best$method), function(m) {
    z <- pred_best[pred_best$method == m, , drop = FALSE]
    data.frame(
      method = m,
      n_cells = nrow(z),
      rmse_gap = mean(z$rmse_gap_mean, na.rm = TRUE),
      rmse_gap_se = se(z$rmse_gap_mean),
      logscore_gap = mean(z$logscore_gap_mean, na.rm = TRUE),
      logscore_gap_se = se(z$logscore_gap_mean),
      selected = mean(z$expected_code_mean, na.rm = TRUE),
      selected_se = se(z$expected_code_mean),
      stringsAsFactors = FALSE
    )
  }))
} else {
  pred_agg <- data.frame()
}
write.csv(pred_agg, file.path(table_dir, "table_large_end_to_end_predictive_baselines.csv"), row.names = FALSE)

tex_path <- file.path(table_dir, "table_large_end_to_end_askpc_family_summary.tex")
con <- file(tex_path, open = "w")
on.exit(close(con), add = TRUE)
writeLines("\\begin{tabular}{lcccccc}", con)
writeLines("\\toprule", con)
writeLines("Method & TV & FKL & $q_0$ & Storage & List & RMSE gap \\\\", con)
writeLines("\\midrule", con)
for (i in seq_len(nrow(agg))) {
  line <- paste(
    agg$method[i],
    fmt_mean_se(agg$tv[i], agg$tv_se[i]),
    fmt_mean_se(agg$fkl[i], agg$fkl_se[i]),
    fmt_mean_se(agg$q0[i], agg$q0_se[i]),
    fmt_mean_se(agg$code[i], agg$code_se[i], digits = 2),
    fmt_mean_se(agg$list_count[i], agg$list_count_se[i], digits = 1),
    fmt_mean_se(agg$rmse_gap[i], agg$rmse_gap_se[i]),
    sep = " & "
  )
  writeLines(paste0(line, " \\\\"), con)
}
writeLines("\\bottomrule", con)
writeLines("\\end{tabular}", con)
close(con)

pred_tex <- file.path(table_dir, "table_large_end_to_end_predictive_baselines.tex")
conp <- file(pred_tex, open = "w")
writeLines("\\begin{tabular}{lccc}", conp)
writeLines("\\toprule", conp)
writeLines("Method & RMSE gap & Log-score gap & Selected \\\\", conp)
writeLines("\\midrule", conp)
if (nrow(pred_agg)) {
  for (i in seq_len(nrow(pred_agg))) {
    line <- paste(
      pred_agg$method[i],
      fmt_mean_se(pred_agg$rmse_gap[i], pred_agg$rmse_gap_se[i]),
      fmt_mean_se(pred_agg$logscore_gap[i], pred_agg$logscore_gap_se[i]),
      fmt_mean_se(pred_agg$selected[i], pred_agg$selected_se[i], digits = 1),
      sep = " & "
    )
    writeLines(paste0(line, " \\\\"), conp)
  }
}
writeLines("\\bottomrule", conp)
writeLines("\\end{tabular}", conp)
close(conp)

status_counts <- table(diag$diagnostic_status)
get_count <- function(name) {
  if (name %in% names(status_counts)) as.integer(status_counts[[name]]) else 0L
}
range_text <- function(x, digits = 3) {
  x <- x[is.finite(x)]
  paste0(fmt_num(mean(x), digits), " [", fmt_num(min(x), digits), ", ", fmt_num(max(x), digits), "]")
}
diag_summary <- data.frame(
  Diagnostic = c(
    "Reference runs",
    "Passed gate",
    "Stress-test label",
    "Low-acceptance warnings",
    "Retained draws per run",
    "Split Rhat max",
    "Minimum ESS",
    "Max group-PIP MCSE",
    "Unique supports"
  ),
  Value = c(
    as.character(nrow(diag)),
    as.character(get_count("pass")),
    as.character(get_count("stress-test")),
    if ("diagnostic_warning" %in% names(diag)) {
      as.character(sum(grepl("low acceptance", diag$diagnostic_warning, fixed = TRUE), na.rm = TRUE))
    } else {
      "0"
    },
    paste(sort(unique(diag$retained_draws)), collapse = ", "),
    range_text(diag$split_rhat_max),
    range_text(diag$ess_min),
    range_text(diag$group_pip_mcse_max),
    range_text(diag$unique_supports, digits = 1)
  ),
  stringsAsFactors = FALSE
)
write.csv(diag_summary, file.path(table_dir, "table_large_end_to_end_reference_diagnostics_summary.csv"), row.names = FALSE)

diag_tex <- file.path(table_dir, "table_large_end_to_end_reference_diagnostics.tex")
con2 <- file(diag_tex, open = "w")
writeLines("\\begin{tabular}{lc}", con2)
writeLines("\\toprule", con2)
writeLines("Diagnostic & Value \\\\", con2)
writeLines("\\midrule", con2)
for (i in seq_len(nrow(diag_summary))) {
  writeLines(paste0(diag_summary$Diagnostic[i], " & ", diag_summary$Value[i], " \\\\"), con2)
}
writeLines("\\bottomrule", con2)
writeLines("\\end{tabular}", con2)
close(con2)

cat("Wrote large end-to-end ASK-PC LaTeX outputs:\n")
cat("  ", tex_path, "\n", sep = "")
cat("  ", diag_tex, "\n", sep = "")
cat("  ", pred_tex, "\n", sep = "")
