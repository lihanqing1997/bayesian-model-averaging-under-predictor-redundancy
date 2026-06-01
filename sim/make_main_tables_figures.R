rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
scripts <- c(
  file.path("sim", "make_q0_screened_frontier_outputs.R"),
  file.path("sim", "make_storage_code_outputs.R"),
  file.path("sim", "make_adaptive_kernel_latex_tables.R"),
  file.path("sim", "redraw_main_figures.R"),
  file.path("sim", "make_family_mixture_latex_tables.R"),
  file.path("sim", "make_jmlr_tables.R"),
  file.path("sim", "make_jmlr_figures.R"),
  file.path("sim", "make_semisynthetic_tecator_kernel_figure.R"),
  file.path("sim", "make_realdata_adaptive_kernel_figures.R")
)
for (script in scripts[file.exists(scripts)]) {
  status <- system2(rscript, script)
  if (!identical(status, 0L)) stop("main table/figure assembly failed for ", script)
}
cat("Main table and figure assembly scripts completed\n")
