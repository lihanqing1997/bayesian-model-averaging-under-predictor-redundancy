rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
status <- system2(rscript, file.path("sim", "run_main_adaptive_kernel_smoke.R"))
if (!identical(status, 0L)) stop("main support-kernel smoke run failed")
