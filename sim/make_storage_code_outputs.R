table_dir <- file.path("sim", "output", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

escape_tex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x
}

storage_defs <- data.frame(
  family = c(
    "Top-M atoms and credible sets",
    "Fixed hard or group-representative regions",
    "Structured interval and graph kernels",
    "Posterior-cluster kernels",
    "Pooled-pruned hybrids",
    "Prior-changing BMA baselines"
  ),
  descriptor = c(
    "Sparse index sets for selected supports",
    "Active groups, capacities, and region type",
    "Intervals or graph communities, capacity, and bandwidth",
    "Cluster medoids, bandwidth, and cluster count",
    "Retained kernel descriptors after q-mass pruning",
    "Prior family, hyperparameters, and posterior approximation"
  ),
  primary_code = c(
    "Weighted sparse-support descriptor",
    "Group-pattern and capacity descriptor",
    "Region descriptor with capacity and bandwidth",
    "Medoid descriptor with bandwidth and cluster overhead",
    "Weighted descriptor cost of the refitted list",
    "Not a compression code for the unrestricted posterior"
  ),
  sensitivity_check = c(
    "Stored atom count and credible-set size",
    "Number of listed hard regions",
    "Active kernel count and total retained list size",
    "Cluster count, medoid size, and fallback weight",
    "Active kernel count, effective kernels, and retained q-mass threshold",
    "Target-shift TV/FKL and predictive diagnostics"
  ),
  stringsAsFactors = FALSE
)

write.csv(storage_defs, file.path(table_dir, "table_storage_code_definitions.csv"), row.names = FALSE)

tex_path <- file.path(table_dir, "table_storage_code_definitions.tex")
con <- file(tex_path, open = "w", encoding = "UTF-8")
on.exit(close(con), add = TRUE)
colspec <- "{@{}>{\\raggedright\\arraybackslash}m{0.19\\textwidth}>{\\raggedright\\arraybackslash}m{0.245\\textwidth}>{\\raggedright\\arraybackslash}m{0.22\\textwidth}>{\\raggedright\\arraybackslash}X@{}}"
table_width <- "0.985\\textwidth"
writeLines(paste0("\\setlength{\\tabcolsep}{2.6pt}"), con)
writeLines(paste0("\\renewcommand{\\tabularxcolumn}[1]{m{#1}}"), con)
writeLines(paste0("\\renewcommand{\\arraystretch}{1.0}"), con)
writeLines(paste0("\\begin{tabularx}{", table_width, "}", colspec), con)
writeLines("\\toprule", con)
writeLines("Family & Stored descriptor & Primary code & Sensitivity check \\\\", con)
writeLines("\\midrule", con)
writeLines("\\end{tabularx}", con)
writeLines("\\vspace{-0.25em}", con)
writeLines("\\renewcommand{\\arraystretch}{1.32}", con)
writeLines(paste0("\\begin{tabularx}{", table_width, "}", colspec), con)
for (i in seq_len(nrow(storage_defs))) {
  line <- paste(
    paste0("\\rule[-1.8ex]{0pt}{7.2ex}", escape_tex(storage_defs$family[i])),
    escape_tex(storage_defs$descriptor[i]),
    escape_tex(storage_defs$primary_code[i]),
    escape_tex(storage_defs$sensitivity_check[i]),
    sep = " & "
  )
  writeLines(paste0(line, " \\\\"), con)
}
writeLines("\\bottomrule", con)
writeLines("\\end{tabularx}", con)

fmt <- function(x) {
  if (!is.finite(x)) return(NA_real_)
  x
}

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

fallback_reading <- function(q0) {
  if (!is.finite(q0)) return("not applicable")
  if (q0 <= 0.10) return("low fallback")
  if (q0 <= 0.50) return("partial fallback")
  "large fallback"
}

method_keep <- c(
  "ASK-PC pooled-pruned 99%",
  "Posterior clustering",
  "Top-M support atoms",
  "Credible support set",
  "Dilution-prior BMA",
  "DPP-prior BMA",
  "Fixed hard dictionary"
)

summarize_storage <- function(path, study, code_col, fkl_col, q0_col,
                              active_col = NULL, total_col = NULL,
                              atoms_col = NULL, method_col = "method") {
  if (!file.exists(path)) return(NULL)
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[d[[method_col]] %in% method_keep, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  parts <- split(d, d[[method_col]], drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(z) {
    data.frame(
      study = study,
      method = z[[method_col]][1],
      fkl = mean_or_na(z[[fkl_col]]),
      q0 = mean_or_na(z[[q0_col]]),
      expected_code = mean_or_na(z[[code_col]]),
      active_list = if (!is.null(active_col) && active_col %in% names(z)) mean_or_na(z[[active_col]]) else NA_real_,
      total_list = if (!is.null(total_col) && total_col %in% names(z)) mean_or_na(z[[total_col]]) else NA_real_,
      stored_atoms = if (!is.null(atoms_col) && atoms_col %in% names(z)) mean_or_na(z[[atoms_col]]) else NA_real_,
      fallback_reading = fallback_reading(mean_or_na(z[[q0_col]])),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

storage_sensitivity <- do.call(rbind, Filter(Negate(is.null), list(
  summarize_storage(
    file.path(table_dir, "table_support_kernel_competitor_best_fkl.csv"),
    "Exact benchmark",
    code_col = "code_mean",
    fkl_col = "fkl_mean",
    q0_col = "q0_mean"
  ),
  summarize_storage(
    file.path(table_dir, "table_large_end_to_end_askpc_best_fkl.csv"),
    "Large p=100 benchmark",
    code_col = "expected_code_mean",
    fkl_col = "fkl_mean",
    q0_col = "q0_mean",
    active_col = "active_kernels_001_mean",
    total_col = "n_kernels_mean",
    atoms_col = "stored_atoms_mean"
  ),
  summarize_storage(
    file.path(table_dir, "table_real_response_reduced_best_fkl.csv"),
    "Reduced real-response benchmark",
    code_col = "expected_code_mean",
    fkl_col = "fkl_mean",
    q0_col = "q0_mean",
    active_col = "active_kernels_001_mean",
    total_col = "q_effective_kernels_mean",
    atoms_col = "stored_atoms_mean"
  )
)))

if (!is.null(storage_sensitivity) && nrow(storage_sensitivity)) {
  storage_sensitivity <- storage_sensitivity[order(storage_sensitivity$study, match(storage_sensitivity$method, method_keep)), ]
  write.csv(storage_sensitivity, file.path(table_dir, "table_storage_code_sensitivity.csv"), row.names = FALSE)
}

cat("Wrote storage code definition and sensitivity tables in ", table_dir, "\n", sep = "")
