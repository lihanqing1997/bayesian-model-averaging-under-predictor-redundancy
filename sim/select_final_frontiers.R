table_dir <- file.path("sim", "output", "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}

rmse_tol <- 0.05
log_tol <- 0.03

choose_row <- function(d, mode, alpha_min = NA_real_) {
  if (nrow(d) == 0L) {
    return(NULL)
  }
  ok <- is.finite(d$RMSE_gap) & d$RMSE_gap <= rmse_tol &
    (is.na(d$logscore_gap) | d$logscore_gap >= -log_tol) &
    (!("diagnostic_pass" %in% names(d)) | is.na(d$diagnostic_pass) | d$diagnostic_pass)
  if (mode == "A") {
    ok <- ok & is.finite(d$alpha_bar) & d$alpha_bar >= alpha_min
  }
  cand <- d[ok, , drop = FALSE]
  if (nrow(cand) == 0L) {
    out <- d[1, , drop = FALSE]
    out$status <- "no pass"
    out$mode <- ifelse(mode == "A", "A distortion-controlled", "B prediction-compression")
    out$alpha_min <- alpha_min
    for (nm in c("capacity_b", "lambda", "family95", "support95", "RMSE_gap", "logscore_gap",
                 "alpha_bar", "distortion_bar", "retained_LCB95", "max_Rhat", "ESS_median")) {
      if (!nm %in% names(out)) out[[nm]] <- NA_real_
      out[[nm]] <- NA_real_
    }
    return(out)
  }
  cand <- cand[order(cand$family95, cand$support95, cand$RMSE_gap), , drop = FALSE]
  out <- cand[1, , drop = FALSE]
  out$status <- "pass"
  out$mode <- ifelse(mode == "A", "A distortion-controlled", "B prediction-compression")
  out$alpha_min <- alpha_min
  out
}

rows <- list()

direct <- read_if_exists(file.path(table_dir, "table_synthetic_direct_sampler.csv"))
retained <- read_if_exists(file.path(table_dir, "table_synthetic_recycling_independent.csv"))
if (!is.null(direct) && !is.null(retained)) {
  ret <- retained[retained$top_group_count == 12, ]
  if (nrow(ret) > 0L) {
    ret_mean <- aggregate(
      cbind(alpha_bar_mean, distortion_bar_mean) ~
        scenario + rho + capacity_b + lambda,
      data = ret,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    names(ret_mean)[names(ret_mean) == "alpha_bar_mean"] <- "alpha_bar"
    names(ret_mean)[names(ret_mean) == "distortion_bar_mean"] <- "distortion_bar"
    ret_mean$retained_LCB95 <- NA_real_
  } else {
    ret_mean <- direct[FALSE, c("scenario", "rho", "capacity_b", "lambda"), drop = FALSE]
    ret_mean$alpha_bar <- NA_real_
    ret_mean$distortion_bar <- NA_real_
    ret_mean$retained_LCB95 <- NA_real_
  }
  syn <- direct[direct$method == "repfam_direct", ]
  syn <- merge(syn, ret_mean[, c("scenario", "rho", "capacity_b", "lambda", "alpha_bar", "distortion_bar", "retained_LCB95")],
               by = c("scenario", "rho", "capacity_b", "lambda"), all.x = TRUE)
  syn$dataset <- "synthetic"
  syn$family95 <- syn$family_n95_mean
  syn$support95 <- syn$support_n95_mean
  syn$RMSE_gap <- syn$rmse_gap_mean
  syn$logscore_gap <- syn$log_score_gap_mean
  syn$max_Rhat <- syn$max_group_rhat_mean
  syn$ESS_median <- syn$median_group_ess_mean
  syn$diagnostic_pass <- with(syn, is.finite(max_Rhat) & max_Rhat <= 1.05 & ESS_median >= 400)
  for (key in unique(paste(syn$scenario, syn$rho))) {
    parts <- strsplit(key, " ", fixed = TRUE)[[1]]
    d <- syn[syn$scenario == parts[1] & as.character(syn$rho) == parts[2], , drop = FALSE]
    for (a in c(0.8, 0.9)) rows[[length(rows) + 1L]] <- choose_row(d, "A", a)
    rows[[length(rows) + 1L]] <- choose_row(d, "B", NA_real_)
  }
}

real_ret <- read_if_exists(file.path(table_dir, "table_real_spectroscopy_retained_validation_final.csv"))
real_sum <- NULL
for (dataset in c("tecator", "gasoline")) {
  path <- file.path("sim", "output", sprintf("%s_capacity_benchmark_summary.csv", dataset))
  if (file.exists(path)) {
    d <- read.csv(path, stringsAsFactors = FALSE)
    ub <- d[d$method == "unrestricted_bma", ]
    rf <- d[d$method == "representative_family", ]
    if (nrow(ub) > 0L && nrow(rf) > 0L) {
      rf$dataset <- dataset
      rf$RMSE_gap <- rf$rmse_mean - ub$rmse_mean[1]
      rf$logscore_gap <- rf$log_score_mean - ub$log_score_mean[1]
      real_sum <- rbind(real_sum, rf)
    }
  }
}
if (!is.null(real_sum)) {
  real_sum$capacity_b <- real_sum$capacity
  real_sum$family95 <- real_sum$family_n95_mean
  real_sum$support95 <- real_sum$support_n95_mean
  real_sum$alpha_bar <- real_sum$alpha_bar_mean
  real_sum$distortion_bar <- real_sum$distortion_bar_mean
  real_sum$retained_LCB95 <- NA_real_
  real_sum$max_Rhat <- real_sum$max_group_rhat_mean
  real_sum$ESS_median <- real_sum$median_group_ess_mean
  real_sum$diagnostic_pass <- with(real_sum, is.na(max_Rhat) | (max_Rhat <= 1.05 & ESS_median >= 400))
  if (!is.null(real_ret)) {
    rmini <- aggregate(
      candidate_cover_LCB95_mean ~ dataset + capacity_b + lambda + top_group_count,
      data = real_ret,
      FUN = min,
      na.rm = TRUE
    )
    real_sum <- merge(real_sum, rmini, by = c("dataset", "capacity_b", "lambda", "top_group_count"), all.x = TRUE)
    real_sum$retained_LCB95 <- real_sum$candidate_cover_LCB95_mean
  }
  for (dataset in unique(real_sum$dataset)) {
    d <- real_sum[real_sum$dataset == dataset, , drop = FALSE]
    d$scenario <- "spectroscopy"
    for (a in c(0.8, 0.9)) rows[[length(rows) + 1L]] <- choose_row(d, "A", a)
    rows[[length(rows) + 1L]] <- choose_row(d, "B", NA_real_)
  }
}

rib <- read_if_exists(file.path(table_dir, "table_riboflavin.csv"))
rib_ret <- read_if_exists(file.path(table_dir, "table_riboflavin_retained_validation.csv"))
if (!is.null(rib)) {
  d <- rib[rib$method == "repfam_direct" & rib$p0 == 300 & rib$grouping_K == 50, , drop = FALSE]
  if (!is.null(rib_ret)) {
    rr <- rib_ret[rib_ret$p0 == 300 & rib_ret$grouping_K == 50 & rib_ret$top_group_count == 12, ]
    if (nrow(rr) > 0L) {
      rr2 <- aggregate(
        cbind(alpha_bar_mean, distortion_bar_mean) ~ capacity_b + lambda,
        data = rr,
        FUN = function(x) mean(x, na.rm = TRUE)
      )
      names(rr2)[names(rr2) == "alpha_bar_mean"] <- "alpha_bar"
      names(rr2)[names(rr2) == "distortion_bar_mean"] <- "distortion_bar"
      d <- merge(d, rr2, by = c("capacity_b", "lambda"), all.x = TRUE)
    } else {
      d$alpha_bar <- NA_real_
      d$distortion_bar <- NA_real_
    }
  } else {
    d$alpha_bar <- NA_real_
    d$distortion_bar <- NA_real_
  }
  d$dataset <- "riboflavin"
  d$scenario <- "p0=300,K=50"
  d$family95 <- d$family_n95_mean
  d$support95 <- d$support_n95_mean
  d$RMSE_gap <- d$rmse_gap_mean
  d$logscore_gap <- d$log_score_gap_mean
  d$retained_LCB95 <- NA_real_
  d$max_Rhat <- d$max_group_rhat_mean
  d$ESS_median <- d$median_group_ess_mean
  d$diagnostic_pass <- with(d, is.finite(max_Rhat) & max_Rhat <= 1.05 & ESS_median >= 400)
  for (a in c(0.8, 0.9)) rows[[length(rows) + 1L]] <- choose_row(d, "A", a)
  rows[[length(rows) + 1L]] <- choose_row(d, "B", NA_real_)
}

all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
rows <- lapply(rows, function(x) {
  missing <- setdiff(all_names, names(x))
  for (nm in missing) x[[nm]] <- NA
  x[, all_names, drop = FALSE]
})
frontier <- do.call(rbind, rows)
keep <- c(
  "dataset", "scenario", "rho", "mode", "alpha_min", "capacity_b", "lambda",
  "family95", "support95", "RMSE_gap", "logscore_gap", "alpha_bar",
  "distortion_bar", "retained_LCB95", "max_Rhat", "ESS_median",
  "diagnostic_pass", "status"
)
for (nm in keep) if (!nm %in% names(frontier)) frontier[[nm]] <- NA
frontier <- frontier[, keep]
write.csv(frontier, file.path(table_dir, "table_modeA_modeB_frontier.csv"), row.names = FALSE)
write.csv(frontier, file.path(table_dir, "supp_table_full_frontier_grid.csv"), row.names = FALSE)
message("Mode A / Mode B frontier table generated")
