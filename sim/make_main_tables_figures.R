rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
scripts <- c(
  file.path("sim", "make_large_end_to_end_askpc_outputs.R"),
  file.path("sim", "make_support_kernel_competitor_outputs.R"),
  file.path("sim", "make_reported_summary_tables.R"),
  file.path("sim", "make_key_paired_comparisons.R"),
  file.path("sim", "make_adaptive_kernel_latex_tables.R"),
  file.path("sim", "redraw_main_figures.R"),
  file.path("sim", "make_semisynthetic_tecator_kernel_figure.R"),
  file.path("sim", "redraw_largep_frontiers_clean.R"),
  file.path("sim", "make_model_space_illustration.R"),
  file.path("sim", "make_practical_regime_examples.R")
)
for (script in scripts[file.exists(scripts)]) {
  status <- system2(rscript, script)
  if (!identical(status, 0L)) stop("main table/figure assembly failed for ", script)
}
cat("Main table and figure assembly scripts completed\n")
