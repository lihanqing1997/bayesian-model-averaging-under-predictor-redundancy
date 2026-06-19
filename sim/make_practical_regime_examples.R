output_dir <- file.path("sim", "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

blue <- "#2f6f9f"
orange <- "#ca6f1e"
green <- "#5b8e3e"
purple <- "#7b5aa6"
gray <- "#5c6370"

draw_spectroscopy_panel <- function() {
  x <- seq(0, 1, length.out = 160)
  y <- 0.18 + 0.32 * exp(-((x - 0.34) / 0.12)^2) +
    0.24 * exp(-((x - 0.68) / 0.16)^2)
  plot(
    x,
    y,
    type = "l",
    lwd = 2,
    col = blue,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "Spectroscopy"
  )
  box(col = "gray78")
  abline(v = seq(0.05, 0.95, length.out = 18), col = adjustcolor(gray, 0.18), lwd = 0.8)
  rect(0.27, 0.03, 0.43, 0.63, col = adjustcolor(orange, 0.12), border = NA)
  lines(x, y, lwd = 2, col = blue)
  arrows(0.27, 0.08, 0.43, 0.08, length = 0.06, angle = 25, code = 3, col = orange, lwd = 1.5)
  text(0.35, 0.13, "one absorption region", cex = 0.76, col = orange)
  points(c(0.31, 0.35, 0.39), approx(x, y, c(0.31, 0.35, 0.39))$y, pch = 19, col = orange, cex = 1.2)
  text(0.52, 0.54, "nearby wavelengths\nact as substitutes", cex = 0.78, col = gray)
}

draw_molecular_panel <- function() {
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "Molecular modules"
  )
  box(col = "gray78")
  module <- cbind(
    c(0.25, 0.35, 0.46, 0.38, 0.25, 0.51),
    c(0.62, 0.77, 0.67, 0.50, 0.44, 0.42)
  )
  other <- cbind(c(0.72, 0.82, 0.73), c(0.70, 0.58, 0.40))
  for (i in seq_len(nrow(module))) {
    for (j in seq_len(nrow(module))) {
      if (i < j && abs(i - j) <= 2) {
        segments(module[i, 1], module[i, 2], module[j, 1], module[j, 2],
                 col = adjustcolor(green, 0.28), lwd = 1.2)
      }
    }
  }
  symbols(0.36, 0.57, circles = 0.26, inches = FALSE, add = TRUE,
          bg = adjustcolor(green, 0.10), fg = adjustcolor(green, 0.45), lwd = 1.2)
  points(module, pch = 21, bg = green, col = "white", cex = 1.6, lwd = 1)
  points(other, pch = 21, bg = "#b4bac5", col = "white", cex = 1.3, lwd = 1)
  arrows(0.55, 0.57, 0.70, 0.57, length = 0.08, col = orange, lwd = 1.6)
  text(0.78, 0.57, "stable\nmodule signal", cex = 0.78, col = orange)
  text(0.36, 0.22, "several probes or genes\nshare one biological role", cex = 0.78, col = gray)
}

draw_sensor_panel <- function() {
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "Sensors, space, or lags"
  )
  box(col = "gray78")
  xs <- rep(seq(0.20, 0.72, length.out = 4), each = 3)
  ys <- rep(seq(0.32, 0.72, length.out = 3), times = 4)
  highlight <- xs > 0.34 & xs < 0.60 & ys > 0.44
  symbols(0.47, 0.58, circles = 0.23, inches = FALSE, add = TRUE,
          bg = adjustcolor(purple, 0.10), fg = adjustcolor(purple, 0.42), lwd = 1.2)
  points(xs, ys, pch = 21, bg = ifelse(highlight, purple, "#b4bac5"), col = "white", cex = 1.5, lwd = 1)
  axis_x <- seq(0.18, 0.78, length.out = 7)
  segments(axis_x[-length(axis_x)], 0.16, axis_x[-1], 0.16, col = adjustcolor(gray, 0.45), lwd = 1.2)
  points(axis_x, rep(0.16, length(axis_x)), pch = 21,
         bg = ifelse(seq_along(axis_x) %in% 3:5, purple, "#b4bac5"),
         col = "white", cex = 1.2, lwd = 1)
  text(0.50, 0.08, "neighboring locations or lags\ncan encode the same effect", cex = 0.78, col = gray)
}

draw_figure <- function(path, device) {
  device(path, width = 10.5, height = 3.2)
  op <- par(mfrow = c(1, 3), mar = c(2.2, 1.4, 2.4, 1.0), oma = c(0, 0, 0, 0))
  on.exit({
    par(op)
    dev.off()
  })
  draw_spectroscopy_panel()
  draw_molecular_panel()
  draw_sensor_panel()
}

draw_figure(file.path(output_dir, "practical_regime_examples.svg"), svg)
draw_figure(file.path(output_dir, "practical_regime_examples.pdf"), pdf)
