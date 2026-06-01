table_dir <- file.path("sim", "output", "tables")
fig_dir <- file.path("sim", "output", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

method_labels <- c(
  adaptive_support_kernel = "Adaptive kernels",
  fixed_hard_dictionary = "Fixed hard dictionary",
  topM_support_atoms = "Top-M atoms"
)
method_cols <- c(
  adaptive_support_kernel = "#1F77B4",
  fixed_hard_dictionary = "#D95F02",
  topM_support_atoms = "#2CA02C"
)

draw_rate_distortion <- function() {
  path <- file.path(table_dir, "table_adaptive_kernel_exact_detail.csv")
  if (!file.exists(path)) return(invisible(FALSE))
  exact <- read.csv(path, stringsAsFactors = FALSE)
  exact <- exact[exact$method %in% names(method_labels), , drop = FALSE]
  if (!nrow(exact)) return(invisible(FALSE))
  grDevices::pdf(file.path(fig_dir, "fig_adaptive_kernel_rate_distortion.pdf"), width = 6.7, height = 4.6)
  old <- par(mar = c(4.5, 4.8, 2.0, 1.0), mgp = c(2.7, 0.8, 0))
  on.exit({ par(old); grDevices::dev.off() }, add = TRUE)
  plot(
    exact$expected_code, exact$fkl,
    pch = 19,
    col = method_cols[exact$method],
    xlab = "Expected kernel code length",
    ylab = "Forward KL",
    main = "Fixed versus adaptive support-kernel compression",
    cex = 0.95,
    cex.main = 0.95,
    cex.lab = 0.95
  )
  grid(col = "grey88", lty = "dotted")
  present <- names(method_labels)[names(method_labels) %in% unique(exact$method)]
  legend(
    "topleft",
    inset = c(0.035, 0.035),
    legend = unname(method_labels[present]),
    col = method_cols[present],
    pch = 19,
    bty = "o",
    bg = grDevices::adjustcolor("white", alpha.f = 0.9),
    box.col = "grey85",
    cex = 0.82
  )
  invisible(TRUE)
}

draw_semisynthetic_realx <- function() {
  path <- file.path(table_dir, "table_adaptive_kernel_semisynthetic_realx_detail.csv")
  if (!file.exists(path)) return(invisible(FALSE))
  realx <- read.csv(path, stringsAsFactors = FALSE)
  realx <- realx[realx$method %in% names(method_labels), , drop = FALSE]
  if (!nrow(realx)) return(invisible(FALSE))
  grDevices::pdf(file.path(fig_dir, "fig_adaptive_kernel_semisynthetic_realx.pdf"), width = 6.7, height = 4.5)
  old <- par(mar = c(4.5, 4.8, 2.0, 1.0), mgp = c(2.7, 0.8, 0))
  on.exit({ par(old); grDevices::dev.off() }, add = TRUE)
  plot(
    realx$expected_code, realx$fkl,
    pch = 19,
    col = method_cols[realx$method],
    xlab = "Expected kernel code length",
    ylab = "Forward KL",
    main = "Semi-synthetic Tecator rate-distortion",
    cex = 1.05,
    cex.main = 0.95,
    cex.lab = 0.95
  )
  grid(col = "grey88", lty = "dotted")
  present <- names(method_labels)[names(method_labels) %in% unique(realx$method)]
  legend(
    "topleft",
    inset = c(0.035, 0.035),
    legend = unname(method_labels[present]),
    col = method_cols[present],
    pch = 19,
    bty = "o",
    bg = grDevices::adjustcolor("white", alpha.f = 0.9),
    box.col = "grey85",
    cex = 0.82
  )
  invisible(TRUE)
}

draw_rate_distortion()
draw_semisynthetic_realx()
message("Redrew main rate-distortion figures in ", fig_dir)
