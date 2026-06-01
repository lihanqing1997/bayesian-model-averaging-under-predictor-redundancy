rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
scripts <- c(
  file.path("sim", "make_adaptive_kernel_latex_tables.R"),
  file.path("sim", "make_family_mixture_latex_tables.R")
)
for (script in scripts[file.exists(scripts)]) {
  status <- system2(rscript, script)
  if (!identical(status, 0L)) stop("table/figure assembly failed for ", script)
}
cat("Available table and figure assembly scripts completed\n")
