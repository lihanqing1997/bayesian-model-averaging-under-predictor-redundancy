source(file.path("sim", "src", "support_kernel_benchmark.R"))

args <- commandArgs(trailingOnly = TRUE)
mode_arg <- args[grepl("^--mode=", args)]
mode <- if (length(mode_arg)) sub("^--mode=", "", mode_arg[[1]]) else "smoke"
if (!mode %in% c("smoke", "pilot", "medium", "full")) {
  stop("--mode must be smoke, pilot, medium, or full")
}
checkpoint <- mode == "full" || any(args == "--checkpoint")
resume <- !any(args == "--no-resume")
workers_arg <- args[grepl("^--workers=", args)]
workers <- if (length(workers_arg)) as.integer(sub("^--workers=", "", workers_arg[[1]])) else if (mode == "full") 4L else 1L
chunk_arg <- args[grepl("^--chunk-size=", args)]
chunk_size <- if (length(chunk_arg)) as.integer(sub("^--chunk-size=", "", chunk_arg[[1]])) else if (mode == "full") 8L else 1L

set.seed(20260527)
out <- run_exact_competitor_benchmark(
  mode = mode,
  checkpoint = checkpoint,
  resume = resume,
  workers = workers,
  chunk_size = chunk_size
)

cat("\nWrote competitor benchmark outputs:\n")
cat("  sim/output/tables/table_support_kernel_competitor_benchmark_detail.csv\n")
cat("  sim/output/tables/table_support_kernel_competitor_benchmark_summary.csv\n")
if (checkpoint) {
  cat(sprintf("  sim/output/tables/table_support_kernel_competitor_benchmark_detail_%s_checkpoint.csv\n", mode))
}
cat(sprintf("Rows: detail=%d summary=%d\n", nrow(out$detail), nrow(out$summary)))
cat(sprintf("Workers=%d chunk_size=%d\n", workers, chunk_size))
