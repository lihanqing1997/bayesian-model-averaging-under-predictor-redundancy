make_group_id <- function(K, m) {
  rep(seq_len(K), each = m)
}

standardize_with_training <- function(X_train_raw, X_test_raw) {
  center <- colMeans(X_train_raw)
  scale <- apply(X_train_raw, 2, sd)
  scale[scale == 0] <- 1

  X_train <- sweep(X_train_raw, 2, center, "-")
  X_train <- sweep(X_train, 2, scale, "/")

  X_test <- sweep(X_test_raw, 2, center, "-")
  X_test <- sweep(X_test, 2, scale, "/")

  list(
    X_train = X_train,
    X_test = X_test,
    center = center,
    scale = scale
  )
}

simulate_block_design <- function(n, K, m, rho) {
  p <- K * m
  X <- matrix(0, nrow = n, ncol = p)

  for (k in seq_len(K)) {
    z <- rnorm(n)
    columns <- ((k - 1) * m + 1):(k * m)
    for (j in seq_along(columns)) {
      X[, columns[j]] <- sqrt(rho) * z + sqrt(1 - rho) * rnorm(n)
    }
  }

  X
}

simulate_block_regression <- function(n_train = 60,
                                      n_test = 200,
                                      K = 5,
                                      m = 3,
                                      rho = 0.8,
                                      active_groups = c(1, 3),
                                      active_reps = rep(1, length(active_groups)),
                                      beta_values = c(1.6, 1.2),
                                      sigma = 1,
                                      seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  p <- K * m
  group_id <- make_group_id(K, m)
  beta <- rep(0, p)

  for (i in seq_along(active_groups)) {
    group <- active_groups[i]
    rep_index <- active_reps[i]
    column <- (group - 1) * m + rep_index
    beta[column] <- beta_values[i]
  }

  X_train_raw <- simulate_block_design(n_train, K, m, rho)
  X_test_raw <- simulate_block_design(n_test, K, m, rho)
  standardized <- standardize_with_training(X_train_raw, X_test_raw)

  y_train_raw <- as.numeric(standardized$X_train %*% beta + rnorm(n_train, sd = sigma))
  y_test_raw <- as.numeric(standardized$X_test %*% beta + rnorm(n_test, sd = sigma))

  y_center <- mean(y_train_raw)
  y_train <- y_train_raw - y_center
  y_test <- y_test_raw - y_center

  list(
    X_train = standardized$X_train,
    X_test = standardized$X_test,
    y_train = y_train,
    y_test = y_test,
    y_train_raw = y_train_raw,
    y_test_raw = y_test_raw,
    y_center = y_center,
    beta = beta,
    true_support = as.integer(beta != 0),
    group_id = group_id,
    active_groups = active_groups,
    active_reps = active_reps,
    K = K,
    m = m,
    rho = rho,
    sigma = sigma
  )
}
