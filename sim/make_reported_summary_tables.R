source(file.path("sim", "src", "support_kernel_benchmark.R"))

table_dir <- file.path("sim", "output", "tables")
result_dir <- file.path("sim", "output", "results")
object_dir <- file.path("sim", "output", "large_end_to_end")
tex_table_dir <- file.path("tex", "sim", "output", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tex_table_dir, recursive = TRUE, showWarnings = FALSE)

fmt <- function(x, digits = 3) {
  if (!is.finite(x)) return("--")
  formatC(x, format = "f", digits = digits)
}

tex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([_&#%$])", "\\\\\\1", x, perl = TRUE)
  x
}

cell_path <- function(scenario, rho, replication_id) {
  candidates <- list.files(
    object_dir,
    pattern = sprintf("^full_%s_rho.*_rep%03d_reference_reliable_full_v2\\.rds$", scenario, as.integer(replication_id)),
    full.names = TRUE
  )
  if (!length(candidates)) return(NA_character_)
  want <- as.numeric(rho)
  parsed <- sub(sprintf("^full_%s_rho", scenario), "", basename(candidates))
  parsed <- sub(sprintf("_rep%03d_reference_reliable_full_v2\\.rds$", as.integer(replication_id)), "", parsed)
  parsed <- as.numeric(gsub("p", ".", parsed, fixed = TRUE))
  hit <- which(abs(parsed - want) < 1e-8)
  if (length(hit)) candidates[hit[1]] else NA_character_
}

hard_ball_matrix_local <- function(supports, centers, radii, group_id = NULL, group_level = FALSE) {
  supports <- as.matrix(supports)
  if (group_level) {
    S <- support_group_counts(supports, group_id) > 0
    C <- support_group_counts(centers, group_id) > 0
  } else {
    S <- supports > 0
    C <- centers > 0
  }
  W <- matrix(0, nrow(S), nrow(C) * length(radii))
  k <- 0L
  for (i in seq_len(nrow(C))) {
    dist <- rowSums(sweep(S, 2, C[i, ], `!=`))
    for (r in radii) {
      k <- k + 1L
      W[, k] <- as.numeric(dist <= r)
    }
  }
  W
}

fit_group_hamming_example <- function(path, beta, alpha_floor = 1e-8) {
  obj <- readRDS(path)
  dat <- obj$data
  ref <- obj$reference_fit
  weights <- normalize_weights(ref$posterior)
  ord <- order(weights, decreasing = TRUE)
  centers <- unique(ref$supports[ord[seq_len(min(length(ord), 20L))], , drop = FALSE])
  centers <- centers[seq_len(min(nrow(centers), 5L)), , drop = FALSE]
  radii <- 0:3
  W_ball <- hard_ball_matrix_local(ref$supports, centers, radii, group_id = dat$group_id, group_level = TRUE)
  W <- cbind(1, W_ball)
  alpha <- estimate_family_alpha(W, weights, alpha_floor = alpha_floor)$alpha_truncated
  costs <- c(
    max(10, 4 * log(max(ncol(ref$supports), 2))),
    rep(log(ncol(ref$supports) + 1), ncol(W_ball)) + rep(radii, each = nrow(centers))
  )
  fit <- optimize_family_mixture(
    W,
    alpha,
    costs,
    weights,
    beta = beta,
    tau = 1e-3,
    distortion = "fkl",
    q0_min = 1e-3,
    q0_max = 1,
    max_iter = 120L,
    tol = 1e-7
  )
  list(obj = obj, centers = centers, radii = radii, W = W, alpha = alpha, costs = costs, fit = fit)
}

describe_region <- function(center_groups, radius) {
  if (!length(center_groups)) return("empty center")
  paste0("groups \\{", paste(center_groups, collapse = ","), "\\}")
}

reading_for_radius <- function(radius) {
  if (radius == 0L) return("exact group pattern")
  if (radius == 1L) return("one group change")
  paste0("up to ", radius, " group changes")
}

main_reading <- function(region, radius) {
  if (region == "fallback") return("fallback residual")
  if (radius == 0L) return("exact group pattern")
  paste0("radius ", radius, " group ball")
}

baseline <- read.csv(file.path(result_dir, "same_target_baselines_largep.csv"), stringsAsFactors = FALSE)
gh <- baseline[baseline$method == "Group-Hamming balls" & is.finite(baseline$fkl), , drop = FALSE]
if (!nrow(gh)) stop("No group-Hamming rows found in same_target_baselines_largep.csv")
best <- gh[which.min(gh$fkl), , drop = FALSE]
path <- cell_path(best$scenario, best$rho, best$replication_id)
if (!file.exists(path)) stop("Missing reference object: ", path)

fit <- fit_group_hamming_example(path, beta = best$beta)
obj <- fit$obj
group_active <- support_group_counts(fit$centers, obj$data$group_id) > 0

rows <- list(data.frame(
  region = "fallback",
  kernel_type = "unrestricted BMA",
  center = "all supports",
  radius = NA_integer_,
  q = fit$fit$q[1],
  alpha = fit$alpha[1],
  cost = fit$costs[1],
  reading = "residual unrestricted posterior",
  stringsAsFactors = FALSE
))

k <- 1L
for (i in seq_len(nrow(fit$centers))) {
  center_groups <- which(group_active[i, ])
  for (r in fit$radii) {
    idx <- 1L + k
    rows[[length(rows) + 1L]] <- data.frame(
      region = paste0("R", k),
      kernel_type = "group-Hamming ball",
      center = describe_region(center_groups, r),
      radius = r,
      q = fit$fit$q[idx],
      alpha = fit$alpha[idx],
      cost = fit$costs[idx],
      reading = reading_for_radius(r),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

report <- do.call(rbind, rows)
report <- report[order(report$q, decreasing = TRUE), , drop = FALSE]
report$scenario <- best$scenario
report$rho <- best$rho
report$replication_id <- best$replication_id
report$tv <- best$tv
report$fkl <- best$fkl
report$rkl <- best$rkl
report$q0 <- best$q0
report$cexp <- best$cexp
report$clist <- best$clist
report$klist <- best$klist
report$keff <- best$keff
write.csv(report, file.path(table_dir, "table_largep_group_hamming_report_columns.csv"), row.names = FALSE)

nonfallback <- report[report$region != "fallback", , drop = FALSE]
keys <- paste(nonfallback$center, nonfallback$radius, sep = "||")
display_rows <- lapply(split(nonfallback, keys), function(z) {
  z <- z[order(z$q, decreasing = TRUE), , drop = FALSE]
  z[1, c("kernel_type", "center", "radius", "alpha", "cost", "reading", "scenario", "rho",
         "replication_id", "tv", "fkl", "rkl", "q0", "cexp", "clist", "klist", "keff"), drop = FALSE] |>
    transform(
      q = sum(z$q),
      merged_columns = nrow(z),
      region = NA_character_
    )
})
display <- do.call(rbind, display_rows)
display <- display[order(display$q, decreasing = TRUE), , drop = FALSE]
display$region <- paste0("G", seq_len(nrow(display)))
display <- display[, c("region", "kernel_type", "center", "radius", "merged_columns", "q", "alpha",
                       "cost", "reading", "scenario", "rho", "replication_id", "tv", "fkl", "rkl",
                       "q0", "cexp", "clist", "klist", "keff")]
fallback <- report[report$region == "fallback", , drop = FALSE]
fallback$merged_columns <- 1L
display <- rbind(
  display,
  fallback[, c("region", "kernel_type", "center", "radius", "merged_columns", "q", "alpha",
               "cost", "reading", "scenario", "rho", "replication_id", "tv", "fkl", "rkl",
               "q0", "cexp", "clist", "klist", "keff")]
)
write.csv(display, file.path(table_dir, "table_largep_group_hamming_report.csv"), row.names = FALSE)

active <- display[display$q >= 1e-3, , drop = FALSE]
active_nonfallback <- active[active$region != "fallback", , drop = FALSE]
top_active <- head(active_nonfallback[order(active_nonfallback$q, decreasing = TRUE), , drop = FALSE], 8L)
remaining <- active_nonfallback[!(active_nonfallback$region %in% top_active$region), , drop = FALSE]
main_rows <- top_active
if (nrow(remaining)) {
  main_rows <- rbind(main_rows, data.frame(
    region = sprintf("remaining %d regions", nrow(remaining)),
    kernel_type = "group-Hamming balls",
    center = "see Appendix Table~\\ref{tab:largepghreportfull}",
    radius = NA_integer_,
    merged_columns = sum(remaining$merged_columns),
    q = sum(remaining$q),
    alpha = NA_real_,
    cost = sum(remaining$cost),
    reading = "remaining active report mass",
    scenario = best$scenario,
    rho = best$rho,
    replication_id = best$replication_id,
    tv = best$tv,
    fkl = best$fkl,
    rkl = best$rkl,
    q0 = best$q0,
    cexp = best$cexp,
    clist = best$clist,
    klist = best$klist,
    keff = best$keff,
    stringsAsFactors = FALSE
  ))
}
fallback <- display[display$region == "fallback", , drop = FALSE]
main_rows <- rbind(main_rows, fallback)
main_rows$main_reading <- vapply(seq_len(nrow(main_rows)), function(i) {
  main_reading(main_rows$region[i], main_rows$radius[i])
}, character(1))

main_tex <- c(
  "\\begin{tabularx}{\\textwidth}{@{}lXrrrrrX@{}}",
  "\\toprule",
  "Region & Center & Radius & Kernel cols. & $q_m$ & $\\alpha_m$ & Cost & Meaning \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(main_rows))) {
  main_tex <- c(main_tex, sprintf(
    "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
    tex_escape(main_rows$region[i]),
    main_rows$center[i],
    ifelse(is.na(main_rows$radius[i]), "--", as.character(main_rows$radius[i])),
    ifelse(is.na(main_rows$merged_columns[i]), "--", as.character(main_rows$merged_columns[i])),
    fmt(main_rows$q[i], 3),
    fmt(main_rows$alpha[i], 3),
    fmt(main_rows$cost[i], 2),
    tex_escape(main_rows$main_reading[i])
  ))
}
main_tex <- c(main_tex, "\\bottomrule", "\\end{tabularx}")
writeLines(main_tex, file.path(table_dir, "table_largep_group_hamming_report.tex"))

full_rows <- active[order(active$region == "fallback", -active$q), , drop = FALSE]
full_tex <- c(
  "\\begin{tabularx}{\\textwidth}{@{}llXrrrrr@{}}",
  "\\toprule",
  "Region & Type & Center & Radius & Kernel cols. & $q_m$ & $\\alpha_m$ & Cost \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(full_rows))) {
  full_tex <- c(full_tex, sprintf(
    "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
    tex_escape(full_rows$region[i]),
    tex_escape(full_rows$kernel_type[i]),
    full_rows$center[i],
    ifelse(is.na(full_rows$radius[i]), "--", as.character(full_rows$radius[i])),
    ifelse(is.na(full_rows$merged_columns[i]), "--", as.character(full_rows$merged_columns[i])),
    fmt(full_rows$q[i], 3),
    fmt(full_rows$alpha[i], 3),
    fmt(full_rows$cost[i], 2)
  ))
}
full_tex <- c(full_tex, "\\bottomrule", "\\end{tabularx}")
writeLines(full_tex, file.path(table_dir, "table_largep_group_hamming_report_full.tex"))

for (fn in c(
  "table_largep_group_hamming_report.tex",
  "table_largep_group_hamming_report_full.tex",
  "table_largep_group_hamming_report.csv",
  "table_largep_group_hamming_report_columns.csv"
)) {
  file.copy(file.path(table_dir, fn), file.path(tex_table_dir, fn), overwrite = TRUE)
}

message(sprintf(
  "Wrote group-Hamming report for %s rho %.2f rep %03d with TV %.3f, FKL %.3f, q0 %.3f.",
  best$scenario, best$rho, best$replication_id, best$tv, best$fkl, best$q0
))
