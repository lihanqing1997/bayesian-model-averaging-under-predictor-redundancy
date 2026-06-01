source(file.path("sim", "repfam_workflow.R"))

output_dir <- file.path("sim", "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

read_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

write_if_nonempty <- function(df, path) {
  if (!is.null(df) && nrow(df) > 0L) {
    write.csv(df, path, row.names = FALSE)
  }
}

full <- read_if_exists(file.path(output_dir, "full_simulation_results.csv"))
if (!is.null(full)) {
  full_tests <- compute_paired_tests(
    results = full,
    baseline_method = "unrestricted_bma",
    group_cols = c("scenario", "rho"),
    id_col = "replication",
    pair_cols = c("scenario", "rho", "replication"),
    n_boot = 1000L,
    seed = 20260522L
  )
  write_if_nonempty(full_tests, file.path(output_dir, "full_simulation_paired_tests.csv"))
}

direct <- read_if_exists(file.path(output_dir, "direct_sampler_validation_results.csv"))
if (!is.null(direct)) {
  direct_tests <- compute_paired_tests(
    results = direct,
    baseline_method = "unrestricted_bma",
    group_cols = c("experiment"),
    id_col = "replicate",
    n_boot = 1000L,
    seed = 20260523L
  )
  write_if_nonempty(direct_tests, file.path(output_dir, "direct_sampler_paired_tests.csv"))
}

direct_full <- read_if_exists(file.path(output_dir, "tables", "table_synthetic_direct_sampler_detail.csv"))
if (!is.null(direct_full)) {
  direct_full_tests <- compute_paired_tests(
    results = direct_full,
    baseline_method = "unrestricted_bma",
    group_cols = c("scenario", "rho", "capacity_b", "lambda"),
    id_col = "replication",
    pair_cols = c("scenario", "rho", "replication"),
    n_boot = 1000L,
    seed = 20260525L
  )
  write_if_nonempty(direct_full_tests, file.path(output_dir, "tables", "table_synthetic_direct_sampler_paired_tests.csv"))
}

riboflavin <- read_if_exists(file.path(output_dir, "tables", "table_riboflavin_detail.csv"))
if (!is.null(riboflavin)) {
  riboflavin_tests <- compute_paired_tests(
    results = riboflavin,
    baseline_method = "unrestricted_bma",
    group_cols = c("p0", "grouping_K", "capacity_b", "lambda"),
    id_col = "split",
    pair_cols = c("p0", "grouping_K", "split"),
    n_boot = 1000L,
    seed = 20260526L
  )
  write_if_nonempty(riboflavin_tests, file.path(output_dir, "tables", "table_riboflavin_paired_tests.csv"))
}

real_paths <- c(
  file.path(output_dir, "tecator_capacity_benchmark_results.csv"),
  file.path(output_dir, "gasoline_capacity_benchmark_results.csv")
)
real <- do.call(rbind, lapply(real_paths, read_if_exists))
if (!is.null(real) && nrow(real) > 0L) {
  real_tests <- compute_paired_tests(
    results = real,
    baseline_method = "unrestricted_bma",
    methods = setdiff(unique(real$method), "unrestricted_bma"),
    group_cols = c("dataset", "capacity", "lambda", "top_group_count"),
    id_col = "split",
    pair_cols = c("dataset", "split"),
    n_boot = 1000L,
    seed = 20260524L
  )
  write_if_nonempty(real_tests, file.path(output_dir, "real_benchmark_paired_tests_bootstrap.csv"))
}

message("Paired test computation complete. Outputs are written to sim/output when source files are present.")
