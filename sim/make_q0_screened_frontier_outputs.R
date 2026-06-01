table_dir <- file.path("sim", "output", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

thresholds <- c(0.05, 0.10, 0.25, 0.50)
method_keep <- c(
  "ASK-PC pooled-pruned 99%",
  "Posterior clustering",
  "Top-M support atoms",
  "Credible support set",
  "Fixed hard dictionary"
)
pretty_method <- c(
  "ASK-PC pooled-pruned 99%" = "Pooled-pruned",
  "Posterior clustering" = "Cluster kernels",
  "Top-M support atoms" = "Top-M atoms",
  "Credible support set" = "Credible set",
  "Fixed hard dictionary" = "Fixed regions"
)

read_detail <- function(path, evidence, id_cols) {
  if (!file.exists(path)) return(NULL)
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$method %in% method_keep & is.finite(d$fkl) & is.finite(d$q0), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$evidence <- evidence
  d$cell <- do.call(paste, c(d[id_cols], sep = "::"))
  d
}

se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

summarize_frontier <- function(d) {
  rows <- list()
  for (thr in thresholds) {
    for (m in method_keep) {
      z <- d[d$method == m & d$q0 <= thr, , drop = FALSE]
      if (nrow(z)) {
        parts <- split(z, z$cell, drop = TRUE)
        best <- do.call(rbind, lapply(parts, function(x) x[which.min(x$fkl), , drop = FALSE]))
        row <- data.frame(
          evidence = d$evidence[1],
          q0_max = thr,
          method = m,
          n_cells = length(unique(best$cell)),
          tv_mean = mean(best$tv, na.rm = TRUE),
          tv_se = se(best$tv),
          fkl_mean = mean(best$fkl, na.rm = TRUE),
          fkl_se = se(best$fkl),
          q0_mean = mean(best$q0, na.rm = TRUE),
          q0_se = se(best$q0),
          expected_code_mean = mean(best$expected_code, na.rm = TRUE),
          expected_code_se = se(best$expected_code),
          active_kernels_mean = if ("active_kernels_001" %in% names(best)) mean(best$active_kernels_001, na.rm = TRUE) else NA_real_,
          stored_atoms_mean = if ("stored_atoms" %in% names(best)) mean(best$stored_atoms, na.rm = TRUE) else NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        row <- data.frame(
          evidence = d$evidence[1],
          q0_max = thr,
          method = m,
          n_cells = 0L,
          tv_mean = NA_real_,
          tv_se = NA_real_,
          fkl_mean = NA_real_,
          fkl_se = NA_real_,
          q0_mean = NA_real_,
          q0_se = NA_real_,
          expected_code_mean = NA_real_,
          expected_code_se = NA_real_,
          active_kernels_mean = NA_real_,
          stored_atoms_mean = NA_real_,
          stringsAsFactors = FALSE
        )
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  do.call(rbind, rows)
}

fmt <- function(x, digits = 3) {
  if (!is.finite(x)) return("--")
  formatC(x, format = "f", digits = digits)
}

fmt_mean <- function(mu, se_val, digits = 3) {
  if (!is.finite(mu)) return("--")
  if (!is.finite(se_val)) return(fmt(mu, digits))
  paste0(fmt(mu, digits), " (", fmt(se_val, digits), ")")
}

details <- Filter(Negate(is.null), list(
  read_detail(file.path(table_dir, "table_support_kernel_competitor_benchmark_detail.csv"), "Exact", c("scenario", "rho", "replication_id")),
  read_detail(file.path(table_dir, "table_large_end_to_end_askpc_detail_reliable_full_v2_augmented.csv"), "Large p=100", c("scenario", "rho", "replication_id")),
  read_detail(file.path(table_dir, "table_real_response_reduced_detail.csv"), "Reduced real response", c("dataset", "replication_id"))
))

if (!length(details)) stop("No detail files available for q0-screened frontier")
frontier <- do.call(rbind, lapply(details, summarize_frontier))
frontier$method_label <- unname(pretty_method[frontier$method])
write.csv(frontier, file.path(table_dir, "table_q0_screened_frontier.csv"), row.names = FALSE)

main <- frontier[frontier$evidence == "Large p=100" & frontier$q0_max %in% c(0.05, 0.10, 0.25, 0.50) & frontier$method %in% c("ASK-PC pooled-pruned 99%", "Posterior clustering", "Top-M support atoms", "Fixed hard dictionary"), , drop = FALSE]
main$method <- factor(main$method, levels = c("ASK-PC pooled-pruned 99%", "Posterior clustering", "Top-M support atoms", "Fixed hard dictionary"))
main <- main[order(main$q0_max, main$method), , drop = FALSE]

tex <- c(
  "\\begin{tabular}{llcccc}",
  "\\toprule",
  "$q_0$ screen & Method & cells & FKL & Storage & active/atoms \\\\",
  "\\midrule"
)
for (thr in unique(main$q0_max)) {
  z <- main[main$q0_max == thr, , drop = FALSE]
  for (i in seq_len(nrow(z))) {
    screen <- if (i == 1L) paste0("$q_0\\le", fmt(thr, 2), "$") else ""
    list_metric <- if (is.finite(z$active_kernels_mean[i])) z$active_kernels_mean[i] else z$stored_atoms_mean[i]
    tex <- c(
      tex,
      sprintf(
        "%s & %s & %d & %s & %s & %s \\\\",
        screen,
        z$method_label[i],
        z$n_cells[i],
        fmt_mean(z$fkl_mean[i], z$fkl_se[i]),
        fmt_mean(z$expected_code_mean[i], z$expected_code_se[i], 2),
        fmt(list_metric, 1)
      )
    )
  }
  if (thr != tail(unique(main$q0_max), 1)) tex <- c(tex, "\\addlinespace")
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(table_dir, "table_q0_screened_frontier.tex"))

cat("Wrote q0-screened frontier outputs\n")
