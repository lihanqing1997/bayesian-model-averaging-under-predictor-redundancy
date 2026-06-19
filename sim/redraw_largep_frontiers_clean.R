result_dir <- file.path("sim", "output", "results")
figure_dir <- file.path("sim", "output", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(result_dir, "largep_frontiers_clean_data.csv")
if (!file.exists(data_path)) {
  stop("Missing ", data_path, ". Run sim/run_largep_report_diagnostics.R first.")
}

plot_data <- read.csv(data_path, stringsAsFactors = FALSE)
required <- c("family", "fkl", "cexp", "clist", "q0")
missing <- setdiff(required, names(plot_data))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

draw_frontier <- function(device, file) {
  device(file, width = 10.2, height = 3.8)
  old <- par(mfrow = c(1, 3), mar = c(3.45, 4.15, 1.35, 0.65), oma = c(0, 0, 0, 0), las = 1, pty = "s")
  on.exit({ par(old); dev.off() }, add = TRUE)

  cols <- c("Pooled-pruned" = "#1b9e77",
            "Group-Hamming balls" = "#d95f02",
            "Hamming balls" = "#7570b3")
  pchs <- c("Pooled-pruned" = 17,
            "Group-Hamming balls" = 19,
            "Hamming balls" = 15)

  panel <- function(x, xlab, main, qline = FALSE, legend_panel = FALSE) {
    xpad <- diff(range(x, finite = TRUE)) * 0.08
    ypad <- diff(range(plot_data$fkl, finite = TRUE)) * 0.10
    if (!is.finite(xpad) || xpad == 0) xpad <- 1
    if (!is.finite(ypad) || ypad == 0) ypad <- 0.01
    plot(x, plot_data$fkl,
         pch = pchs[plot_data$family],
         col = cols[plot_data$family],
         xlab = xlab,
         ylab = "FKL",
         main = main,
         cex = 1.35,
         cex.lab = 1.04,
         cex.axis = 0.92,
         cex.main = 1.08,
         xlim = range(x, finite = TRUE) + c(-xpad, xpad),
         ylim = range(plot_data$fkl, finite = TRUE) + c(-ypad, ypad))
    grid(col = "gray90")
    if (qline) abline(v = 0.25, lty = 2, col = "gray45")
    if (legend_panel) {
      legend("topleft",
             inset = c(0.035, 0.025),
             legend = names(cols),
             col = cols,
             pch = pchs,
             bty = "n",
             cex = 0.64,
             pt.cex = 0.92)
    }
  }

  panel(plot_data$cexp, expression(C[exp]), "A", legend_panel = TRUE)
  panel(plot_data$clist, expression(C[list]), "B")
  panel(plot_data$q0, expression(q[0]), "C", qline = TRUE)
}

draw_frontier(grDevices::pdf, file.path(figure_dir, "largep_frontiers_clean.pdf"))
draw_frontier(function(file, width, height) {
  grDevices::png(file, width = width, height = height, units = "in", res = 300)
}, file.path(figure_dir, "largep_frontiers_clean.png"))

tex_figure_dir <- file.path("tex", "sim", "output", "figures")
if (dir.exists(tex_figure_dir)) {
  file.copy(file.path(figure_dir, "largep_frontiers_clean.pdf"),
            file.path(tex_figure_dir, "largep_frontiers_clean.pdf"),
            overwrite = TRUE)
  file.copy(file.path(figure_dir, "largep_frontiers_clean.png"),
            file.path(tex_figure_dir, "largep_frontiers_clean.png"),
            overwrite = TRUE)
}

cat("Redrew large-p frontier figure:\n")
cat("  ", file.path(figure_dir, "largep_frontiers_clean.pdf"), "\n", sep = "")
cat("  ", file.path(figure_dir, "largep_frontiers_clean.png"), "\n", sep = "")
if (dir.exists(tex_figure_dir)) {
  cat("Synced TeX-local copies:\n")
  cat("  ", file.path(tex_figure_dir, "largep_frontiers_clean.pdf"), "\n", sep = "")
  cat("  ", file.path(tex_figure_dir, "largep_frontiers_clean.png"), "\n", sep = "")
}
