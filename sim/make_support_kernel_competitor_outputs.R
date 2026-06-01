args <- commandArgs(trailingOnly = TRUE)
mode_arg <- args[grepl("^--mode=", args)]
mode <- if (length(mode_arg)) sub("^--mode=", "", mode_arg[[1]]) else NA_character_

table_dir <- file.path("sim", "output", "tables")
fig_dir <- file.path("sim", "output", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

detail_path <- file.path(table_dir, "table_support_kernel_competitor_benchmark_detail.csv")
if (!file.exists(detail_path)) {
  stop("Missing detail CSV: ", detail_path)
}
detail <- read.csv(detail_path, stringsAsFactors = FALSE)
if (!is.na(mode)) {
  detail <- detail[detail$mode == mode, , drop = FALSE]
}
if (!nrow(detail)) {
  stop("No benchmark rows available for requested mode")
}

method_order <- c(
  "ASK-PC pooled-pruned 99%",
  "ASK-PC pooled union",
  "Posterior clustering",
  "Top-M support atoms",
  "Credible support set",
  "Dilution-prior BMA",
  "DPP-prior BMA",
  "Fixed hard dictionary"
)
main_table_methods <- setdiff(method_order, "ASK-PC pooled union")
method_short <- c(
  "ASK-PC pooled-pruned 99%" = "Pooled-pruned",
  "ASK-PC pooled union" = "Pooled union",
  "Posterior clustering" = "Cluster kernels",
  "Top-M support atoms" = "Top-M",
  "Credible support set" = "Credible set",
  "Dilution-prior BMA" = "Dilution",
  "DPP-prior BMA" = "DPP",
  "Fixed hard dictionary" = "Fixed hard"
)
method_cols <- c(
  "ASK-PC pooled-pruned 99%" = "#1F77B4",
  "ASK-PC pooled union" = "#0B4F8A",
  "Posterior clustering" = "#CC6677",
  "Top-M support atoms" = "#117733",
  "Credible support set" = "#44AA99",
  "Dilution-prior BMA" = "#AA4499",
  "DPP-prior BMA" = "#882255",
  "Fixed hard dictionary" = "#DDCC77"
)

select_best_by_rep <- function(d, criterion = "fkl") {
  split_vars <- interaction(d$scenario, d$rho, d$replication_id, d$method, drop = TRUE)
  parts <- split(d, split_vars)
  out <- do.call(rbind, lapply(parts, function(z) {
    if ("status" %in% names(z)) {
      ok <- z$status == "ok" | is.na(z$status)
      z_ok <- z[ok & is.finite(z[[criterion]]), , drop = FALSE]
    } else {
      z_ok <- z[is.finite(z[[criterion]]), , drop = FALSE]
    }
    if (!nrow(z_ok)) {
      return(z[1, , drop = FALSE])
    }
    z_ok[which.min(z_ok[[criterion]]), , drop = FALSE]
  }))
  rownames(out) <- NULL
  out
}

se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

posterior_detail <- detail[detail$method_family != "predictive selection" | is.na(detail$method_family), , drop = FALSE]
best <- select_best_by_rep(posterior_detail, "fkl")
best$method <- factor(best$method, levels = method_order)
best <- best[!is.na(best$method), , drop = FALSE]

summary_parts <- split(best, interaction(best$scenario, best$method, drop = TRUE), drop = TRUE)
compact <- do.call(rbind, lapply(summary_parts, function(d) {
  data.frame(
    scenario = d$scenario[1],
    method = as.character(d$method[1]),
    tv_mean = mean(d$tv, na.rm = TRUE),
    tv_se = se(d$tv),
    fkl_mean = mean(d$fkl, na.rm = TRUE),
    fkl_se = se(d$fkl),
    rkl_mean = mean(d$rkl, na.rm = TRUE),
    rkl_se = se(d$rkl),
    q0_mean = mean(d$q0, na.rm = TRUE),
    q0_se = se(d$q0),
    code_mean = mean(d$expected_code, na.rm = TRUE),
    code_se = se(d$expected_code),
    rmse_gap_mean = mean(d$rmse_gap, na.rm = TRUE),
    rmse_gap_se = se(d$rmse_gap),
    n = nrow(d),
    mode = d$mode[1],
    stringsAsFactors = FALSE
  )
}))
compact <- compact[order(compact$scenario, match(compact$method, method_order)), ]
write.csv(compact, file.path(table_dir, "table_support_kernel_competitor_best_fkl.csv"), row.names = FALSE)
compact_main <- compact[compact$method %in% main_table_methods, , drop = FALSE]

predictive <- detail[detail$method_family == "predictive selection", , drop = FALSE]
if (nrow(predictive)) {
  predictive_parts <- split(predictive, interaction(predictive$scenario, predictive$method, drop = TRUE), drop = TRUE)
  predictive_summary <- do.call(rbind, lapply(predictive_parts, function(d) {
    data.frame(
      scenario = d$scenario[1],
      method = d$method[1],
      rmse_mean = mean(d$rmse, na.rm = TRUE),
      rmse_se = se(d$rmse),
      rmse_gap_mean = mean(d$rmse_gap, na.rm = TRUE),
      rmse_gap_se = se(d$rmse_gap),
      logscore_gap_mean = mean(d$logscore_gap, na.rm = TRUE),
      logscore_gap_se = se(d$logscore_gap),
      selected_mean = mean(d$expected_code, na.rm = TRUE),
      selected_se = se(d$expected_code),
      n = nrow(d),
      mode = d$mode[1],
      stringsAsFactors = FALSE
    )
  }))
  predictive_summary <- predictive_summary[order(predictive_summary$scenario, predictive_summary$method), ]
} else {
  predictive_summary <- data.frame()
}
write.csv(predictive_summary, file.path(table_dir, "table_support_kernel_predictive_baselines.csv"), row.names = FALSE)

fmt <- function(mean, se_val, digits = 3) {
  if (!is.finite(mean)) return("--")
  if (!is.finite(se_val)) return(sprintf(paste0("%.", digits, "f"), mean))
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean, se_val)
}
fmt_code <- function(mean, se_val) fmt(mean, se_val, digits = 2)

tex_path <- file.path(table_dir, "table_support_kernel_competitor_best_fkl.tex")
con <- file(tex_path, open = "wt")
writeLines("\\begin{tabular}{llcccc}", con)
writeLines("\\toprule", con)
writeLines("Scenario & Method & TV & FKL & $q_0$ & Code \\\\", con)
writeLines("\\midrule", con)
for (sc in unique(compact_main$scenario)) {
  z <- compact_main[compact_main$scenario == sc, , drop = FALSE]
  for (i in seq_len(nrow(z))) {
    scenario_cell <- if (i == 1L) gsub("_", " ", sc, fixed = TRUE) else ""
    q0_txt <- if (is.nan(z$q0_mean[i])) "--" else fmt(z$q0_mean[i], z$q0_se[i], digits = 3)
    line <- sprintf(
      "%s & %s & %s & %s & %s & %s \\\\",
      scenario_cell,
      method_short[[z$method[i]]],
      fmt(z$tv_mean[i], z$tv_se[i], digits = 3),
      fmt(z$fkl_mean[i], z$fkl_se[i], digits = 3),
      q0_txt,
      fmt_code(z$code_mean[i], z$code_se[i])
    )
    writeLines(line, con)
  }
  if (sc != tail(unique(compact_main$scenario), 1)) {
    writeLines("\\addlinespace", con)
  }
}
writeLines("\\bottomrule", con)
writeLines("\\end{tabular}", con)
close(con)

pred_tex_path <- file.path(table_dir, "table_support_kernel_predictive_baselines.tex")
pred_con <- file(pred_tex_path, open = "wt")
writeLines("\\begin{tabular}{llccc}", pred_con)
writeLines("\\toprule", pred_con)
writeLines("Scenario & Method & RMSE gap & Log-score gap & Selected \\\\", pred_con)
writeLines("\\midrule", pred_con)
if (nrow(predictive_summary)) {
  for (sc in unique(predictive_summary$scenario)) {
    z <- predictive_summary[predictive_summary$scenario == sc, , drop = FALSE]
    for (i in seq_len(nrow(z))) {
      scenario_cell <- if (i == 1L) gsub("_", " ", sc, fixed = TRUE) else ""
      writeLines(sprintf(
        "%s & %s & %s & %s & %s \\\\",
        scenario_cell,
        z$method[i],
        fmt(z$rmse_gap_mean[i], z$rmse_gap_se[i], digits = 3),
        fmt(z$logscore_gap_mean[i], z$logscore_gap_se[i], digits = 3),
        fmt(z$selected_mean[i], z$selected_se[i], digits = 1)
      ), pred_con)
    }
    if (sc != tail(unique(predictive_summary$scenario), 1)) {
      writeLines("\\addlinespace", pred_con)
    }
  }
}
writeLines("\\bottomrule", pred_con)
writeLines("\\end{tabular}", pred_con)
close(pred_con)

ask <- detail[grepl("^ASK-PC", detail$method), , drop = FALSE]
pooled <- detail[detail$method == "ASK-PC pooled union", , drop = FALSE]
num <- function(x) as.numeric(x)
fmt_diag <- function(x, digits = 3) {
  if (!is.finite(x)) return("--")
  formatC(x, digits = digits, format = "f")
}
median_finite <- function(x) {
  x <- num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}
diagnostics <- data.frame(
  Diagnostic = c(
    "Benchmark cells",
    "Method rows",
    "Failed method rows",
    "ASK-PC rows checked",
    "ASK-PC objective-consistency failures",
    "Pooled-union fallback rows",
    "Median raw ASK-PC candidates",
    "Median near-deduped ASK-PC candidates",
    "Median pruned KKT residual",
    "Maximum pruned KKT residual"
  ),
  Value = c(
    length(unique(interaction(detail$scenario, detail$rho, detail$replication_id, drop = TRUE))),
    nrow(detail),
    sum(detail$status != "ok", na.rm = TRUE),
    nrow(ask),
    sum(ask$objective_consistent != "TRUE", na.rm = TRUE),
    sprintf("%d of %d", sum(pooled$optimizer_status == "embedded_pruned_feasible_solution", na.rm = TRUE), nrow(pooled)),
    fmt_diag(median_finite(ask$raw_candidates), 1),
    fmt_diag(median_finite(ask$near_deduped_candidates), 1),
    fmt_diag(median_finite(detail$kkt_residual[detail$method == "ASK-PC pooled-pruned 99%"]), 4),
    fmt_diag(max(num(detail$kkt_residual[detail$method == "ASK-PC pooled-pruned 99%"]), na.rm = TRUE), 4)
  ),
  stringsAsFactors = FALSE
)
write.csv(diagnostics, file.path(table_dir, "table_support_kernel_algorithm_diagnostics.csv"), row.names = FALSE)

diag_tex_path <- file.path(table_dir, "table_support_kernel_algorithm_diagnostics.tex")
diag_con <- file(diag_tex_path, open = "wt")
writeLines("\\begin{tabular}{lr}", diag_con)
writeLines("\\toprule", diag_con)
writeLines("Diagnostic & Value \\\\", diag_con)
writeLines("\\midrule", diag_con)
for (i in seq_len(nrow(diagnostics))) {
  writeLines(sprintf("%s & %s \\\\", diagnostics$Diagnostic[i], diagnostics$Value[i]), diag_con)
}
writeLines("\\bottomrule", diag_con)
writeLines("\\end{tabular}", diag_con)
close(diag_con)

fig_path <- file.path(fig_dir, "fig_support_kernel_competitor_frontier.pdf")
pdf(fig_path, width = 7.2, height = 7.4)
plot_detail <- detail[
  detail$method %in% main_table_methods &
    is.finite(detail$expected_code) &
    is.finite(detail$fkl),
  ,
  drop = FALSE
]
scenarios <- unique(plot_detail$scenario)
nr <- ceiling(length(scenarios) / 2)
op <- par(mfrow = c(nr, 2), mar = c(3.4, 3.7, 1.8, 0.6), oma = c(0, 0, 0, 0), mgp = c(2.2, 0.65, 0))
for (sc in scenarios) {
  z <- plot_detail[plot_detail$scenario == sc, , drop = FALSE]
  xlim <- range(z$expected_code[is.finite(z$expected_code)], na.rm = TRUE)
  ylim <- range(z$fkl[is.finite(z$fkl)], na.rm = TRUE)
  plot(
    NA,
    xlim = xlim,
    ylim = ylim,
    xlab = "Code or storage proxy",
    ylab = "Forward KL",
    main = gsub("_", " ", sc, fixed = TRUE),
    cex.main = 0.9
  )
  grid(col = "grey88", lty = "dotted")
  for (method in method_order) {
    zz <- z[z$method == method, , drop = FALSE]
    if (!nrow(zz)) next
    points(
      zz$expected_code,
      zz$fkl,
      pch = if (method == "ASK-PC pooled-pruned 99%") 19 else 21,
      bg = method_cols[[method]],
      col = method_cols[[method]],
      cex = 0.9
    )
  }
  if (identical(sc, scenarios[1])) {
    present <- main_table_methods[main_table_methods %in% z$method]
    legend(
      "topright",
      inset = c(0.02, 0.02),
      legend = unname(method_short[present]),
      pt.bg = method_cols[present],
      col = method_cols[present],
      pch = ifelse(present == "ASK-PC pooled-pruned 99%", 19, 21),
      bty = "o",
      bg = grDevices::adjustcolor("white", alpha.f = 0.92),
      box.col = "grey80",
      cex = 0.68
    )
  }
}
if (length(scenarios) %% 2 == 1) plot.new()
par(op)
dev.off()

cat("Wrote:\n")
cat("  ", file.path(table_dir, "table_support_kernel_competitor_best_fkl.csv"), "\n", sep = "")
cat("  ", tex_path, "\n", sep = "")
cat("  ", file.path(table_dir, "table_support_kernel_algorithm_diagnostics.csv"), "\n", sep = "")
cat("  ", file.path(table_dir, "table_support_kernel_predictive_baselines.csv"), "\n", sep = "")
cat("  ", fig_path, "\n", sep = "")
