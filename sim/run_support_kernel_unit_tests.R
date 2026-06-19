source(file.path("sim", "src", "family_mixture_compression.R"))
source(file.path("sim", "src", "support_kernels.R"))
source(file.path("sim", "src", "adaptive_support_kernel_compression.R"))

fail <- function(msg) stop(msg, call. = FALSE)
expect_true <- function(x, msg) if (!isTRUE(x)) fail(msg)
expect_lt <- function(x, y, msg) if (!(is.finite(x) && x < y)) fail(msg)

all_supports <- function(p) {
  as.matrix(expand.grid(rep(list(0:1), p)))
}

set.seed(20260527)
supports <- all_supports(4)
group_id <- c(1L, 1L, 2L, 2L)
raw_weights <- exp(-0.5 * rowSums(supports) + 0.8 * supports[, 1] + 0.4 * supports[, 3])
base_weights <- normalize_weights(raw_weights)

dict <- new_support_kernel_dictionary(list(
  safety_kernel(cost = 10),
  hard_active_set_kernel(1L, capacity = 1L, group_id = group_id),
  soft_hamming_kernel(c(1, 0, 1, 0), rho = 1),
  hard_hamming_ball_kernel(c(1, 0, 1, 0), radius = 1L),
  hard_group_hamming_ball_kernel(c(1L), group_id = group_id, radius = 1L)
))
W <- support_kernel_weight_matrix(supports, dict, group_id = group_id)
alpha <- estimate_kernel_alpha(W, base_weights, alpha_floor = 1e-8)$alpha_truncated
q <- c(0.15, 0.35, 0.20, 0.15, 0.15)

id <- verify_mixture_identities_exact(base_weights, W, alpha, q)
expect_lt(max(unlist(id)), 1e-9, "exact density-ratio identities failed")

fit_small <- optimize_family_mixture(
  W, alpha, kernel_costs(dict), base_weights,
  beta = 0.01, tau = 1e-4, distortion = "fkl",
  q0_min = 0.2, safety_index = 1L, max_iter = 120L, tol = 1e-8
)
expect_true(fit_small$q[1] >= 0.2 - 1e-10, "q0 lower bound was not enforced")
expect_lt(fit_small$kkt_residual, 10, "KKT residual is not finite or is implausibly large")

dict_one <- new_support_kernel_dictionary(dict$kernels[1:2])
W_one <- support_kernel_weight_matrix(supports, dict_one, group_id = group_id)
alpha_one <- estimate_kernel_alpha(W_one, base_weights, alpha_floor = 1e-8)$alpha_truncated
cost_one <- kernel_costs(dict_one)
fit_one <- optimize_family_mixture(
  W_one, alpha_one, cost_one, base_weights,
  beta = 0.01, tau = 0, distortion = "fkl",
  q0_min = 0.1, safety_index = 1L, max_iter = 160L, tol = 1e-8
)
fit_two <- optimize_family_mixture(
  W[, 1:3, drop = FALSE], alpha[1:3], kernel_costs(dict)[1:3], base_weights,
  beta = 0.01, tau = 0, distortion = "fkl",
  q_init = c(0.999 * fit_one$q, 0.001),
  q0_min = 0.1, safety_index = 1L, max_iter = 160L, tol = 1e-8
)
expect_true(
  fit_two$objective <= fit_one$objective + 1e-5,
  "adding a column and reoptimizing increased the empirical objective"
)

cover <- kernel_pool_cover_diagnostics(
  W_oracle = W[, 2:5, drop = FALSE],
  W_pool = W,
  base_weights = base_weights,
  alpha_floor = 1e-8,
  costs_oracle = kernel_costs(dict)[2:5],
  costs_pool = kernel_costs(dict)
)
expect_lt(attr(cover, "xi_max"), 1e-9, "oracle-cover diagnostic failed for a pool containing the oracle")

metric_pool <- generate_metric_ball_kernel_pool(supports, base_weights, group_id = group_id, top_centers = 3L, radii = 0:2)
metric_types <- unique(vapply(metric_pool, function(k) k$type, character(1)))
expect_true("hard_hamming_ball" %in% metric_types, "hard Hamming balls were not generated")
expect_true("hard_group_hamming_ball" %in% metric_types, "hard group-Hamming balls were not generated")

phi <- rowSums(supports[, group_id == 1L, drop = FALSE]) > 0
func <- posterior_functional_error(as.numeric(phi), q, W, alpha, base_weights)
expect_true(func$bound_holds, "TV functional bound failed")

blocks <- blocked_validation_summary(rnorm(100), block_lengths = c(5L, 10L))
expect_true(all(blocks$n_blocks >= 2L), "blocked MCSE summary failed")

ask <- fit_askpc_pooled_pruned(
  supports = supports,
  weights = base_weights,
  group_id = group_id,
  X = matrix(rnorm(nrow(supports) * ncol(supports)), nrow(supports), ncol(supports)),
  beta = 0.01,
  tau = 1e-4,
  q0_min = 0.05,
  prune_mass = 0.99,
  mode = "smoke",
  max_iter = 80L
)
expect_true(nrow(ask$summary) == 2L, "ASK-PC pooled-pruned summary is malformed")
expect_true(all(ask$summary$q0 >= 0.05 - 1e-10), "ASK-PC q0 lower bound failed")
expect_true(isTRUE(ask$pool_diagnostics$objective_consistent), "pooled union objective is worse than the pruned subset")
expect_true(ask$retained_q_mass >= ask$pooled$fit$q[1] + 0.99 * (1 - ask$pooled$fit$q[1]) - 1e-8, "pooled-pruned q-mass rule failed")
expect_true(nrow(ask$pool_counts) > 0, "candidate-pool diagnostics were not recorded")
expect_true(nrow(ask$pruned_top_kernels) >= 1L, "top selected kernel summary was not recorded")
expect_true(all(c("optimizer_status", "objective_change", "raw_candidates", "near_deduped_candidates") %in% names(ask$summary)), "ASK-PC summary lacks optimizer or pool diagnostics")

message("support-kernel unit tests passed")
