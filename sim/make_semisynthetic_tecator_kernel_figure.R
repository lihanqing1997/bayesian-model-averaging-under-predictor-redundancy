source(file.path("sim", "design.R"))
source(file.path("sim", "exact_bma.R"))
source(file.path("sim", "src", "family_dictionary.R"))
source(file.path("sim", "src", "adaptive_support_kernel_compression.R"))

table_dir <- file.path("sim", "output", "tables")
fig_dir <- file.path("sim", "output", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

load_tecator_x <- function() {
  if (!requireNamespace("caret", quietly = TRUE)) return(NULL)
  data("tecator", package = "caret")
  as.matrix(absorp)
}

kernel_predictor_cover <- function(kernel, p, group_id) {
  cover <- integer(p)
  if (identical(kernel$type, "soft_interval")) {
    cover[kernel$params$interval] <- 1L
  } else if (identical(kernel$type, "soft_graph")) {
    cover[kernel$params$nodes] <- 1L
  } else if (identical(kernel$type, "soft_hamming") || identical(kernel$type, "hard_hamming_ball")) {
    cover[which(kernel$params$center_support > 0)] <- 1L
  } else if (identical(kernel$type, "soft_group_hamming") || identical(kernel$type, "hard_group_hamming_ball")) {
    cover[group_id %in% kernel$params$center_groups] <- 1L
  } else if (identical(kernel$type, "hard_active")) {
    cover[group_id %in% kernel$params$active_groups] <- 1L
  } else if (identical(kernel$type, "soft_capacity")) {
    cover[group_id %in% kernel$params$active_groups] <- 1L
  }
  cover
}

X_raw <- load_tecator_x()
if (is.null(X_raw)) {
  write.csv(data.frame(status = "caret_tecator_unavailable"),
            file.path(table_dir, "table_semisynthetic_tecator_band_overlap.csv"),
            row.names = FALSE)
  quit(save = "no")
}

set.seed(5001)
n <- nrow(X_raw)
p0 <- 12L
cols <- unique(round(seq(max(1, floor(ncol(X_raw) * 0.25)),
                         min(ncol(X_raw), ceiling(ncol(X_raw) * 0.75)),
                         length.out = p0)))
X0 <- X_raw[, cols, drop = FALSE]
test_id <- sample(seq_len(n), max(10L, floor(0.25 * n)))
train_id <- setdiff(seq_len(n), test_id)
z <- standardize_with_training(X0[train_id, , drop = FALSE], X0[test_id, , drop = FALSE])
true_band <- 5:7
beta <- rep(0, p0)
beta[true_band] <- seq(0.8, 1.0, length.out = length(true_band))
ytr <- as.numeric(z$X_train %*% beta + rnorm(nrow(z$X_train), sd = 1.25))
center <- mean(ytr)
group_id <- make_group_id(as.integer(p0 / 3L), 3L)
fit <- fit_exact_bma(z$X_train, ytr - center, theta = 0.10, tau2 = 4, a0 = 1, b0 = 1)
pool <- make_default_kernel_pool(fit$supports, fit$posterior, group_id, z$X_train, mode = "medium")
adapt <- learn_adaptive_kernel_dictionary(fit$supports, fit$posterior, group_id, z$X_train, pool,
                                          beta = 0.02, max_iter = 10L, mode = "medium")
d <- mixture_distortions(adapt$q, adapt$W, adapt$alpha, fit$posterior)
q <- adapt$q
ord <- order(q, decreasing = TRUE)
ord <- ord[adapt$q[ord] > 1e-4]
ord <- ord[!vapply(adapt$dictionary$kernels[ord], function(k) identical(k$type, "safety"), logical(1))]
ord <- head(ord, min(8L, length(ord)))

true_cover <- integer(p0)
true_cover[true_band] <- 1L
score <- rep(0, p0)
rows <- list()
for (i in seq_along(ord)) {
  k <- adapt$dictionary$kernels[[ord[i]]]
  cover <- kernel_predictor_cover(k, p0, group_id)
  score <- score + q[ord[i]] * cover
  inter <- sum(cover & true_cover)
  union <- sum(cover | true_cover)
  rows[[i]] <- data.frame(
    rank = i,
    kernel = k$name,
    type = k$type,
    q = q[ord[i]],
    covered_predictors = paste(which(cover > 0), collapse = ";"),
    true_band_covered = inter / length(true_band),
    jaccard_with_true_band = ifelse(union > 0, inter / union, 0),
    stringsAsFactors = FALSE
  )
}
detail <- if (length(rows)) do.call(rbind, rows) else data.frame()
selected <- as.integer(score > max(score, 1e-12) * 0.05)
overall <- data.frame(
  rank = 0L,
  kernel = "weighted_union_top_kernels",
  type = "summary",
  q = sum(q[ord]),
  covered_predictors = paste(which(selected > 0), collapse = ";"),
  true_band_covered = sum(selected & true_cover) / length(true_band),
  jaccard_with_true_band = sum(selected & true_cover) / max(1, sum(selected | true_cover)),
  stringsAsFactors = FALSE
)
out <- rbind(overall, detail)
out$tv <- d$tv
out$fkl <- d$kl_base_to_compressed
out$q0 <- q[1]
out$expected_code <- sum(q * kernel_costs(adapt$dictionary))
out$mode <- "medium"
write.csv(out, file.path(table_dir, "table_semisynthetic_tecator_band_overlap.csv"), row.names = FALSE)

pdf(file.path(fig_dir, "fig_semisynthetic_tecator_kernels.pdf"), width = 7.2, height = 4.1)
op <- par(mar = c(4.5, 4.5, 2.4, 1), mgp = c(2.8, 0.8, 0))
x <- seq_len(p0)
plot(x, score, type = "h", lwd = 7, lend = "butt", col = "#2A6FBB",
     xlab = "Reduced Tecator channel index",
     ylab = "Adaptive kernel q-weight cover",
     main = "Learned support kernels in semi-synthetic Tecator",
     ylim = c(0, max(score, 1e-6) * 1.45),
     cex.main = 0.95,
     cex.lab = 0.95)
rect(min(true_band) - 0.5, par("usr")[3], max(true_band) + 0.5, par("usr")[4],
     col = adjustcolor("#E15759", alpha.f = 0.15), border = NA)
lines(x, score, type = "h", lwd = 7, lend = "butt", col = "#2A6FBB")
text(mean(true_band), max(score, 1e-6) * 1.31, "Simulated active band",
     col = "#9B3A3A", cex = 0.72)
grid(nx = NA, ny = NULL)
par(op)
dev.off()

message("Semi-synthetic Tecator kernel figure written to ", file.path(fig_dir, "fig_semisynthetic_tecator_kernels.pdf"))
