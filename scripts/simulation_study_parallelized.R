library(cli)
library(progressr)
library(future)
library(future.apply)
library(ggplot2)
library(Matrix)
library(pcd)

# resolves futures asynchronously (in parallel) in separate R sessions running
# in the background on the same machine
plan(multisession, workers = parallel::detectCores() - 2L)

# how is progress reported
handlers("cli")

# ─────────────────────────────────────────────────────────────
# Parameter grid
# ─────────────────────────────────────────────────────────────
data_settings <- data.frame(
  n = c(1000, 100),
  p = c(50, 300),
  stringsAsFactors = FALSE
)

other_factors <- expand.grid(
  design = c("independent", "correlated", "block", "ar1"),
  family = c("binomial", "gaussian"),
  k = c(5, 20),
  rho = c(0.5, 0.9),
  error_dist = c("normal", "t", "uniform"),
  seed = 1:10,
  stringsAsFactors = FALSE
)

sim_params <- merge(data_settings, other_factors, by = NULL)

pcd_params <- expand.grid(
  method = c("gradient", "newton"),
  prox_fun = c("soft", "mcp"),
  gamma = 1.7,
  intercept_update = c(
    "after_each_gradient",
    "after_each_newton",
    "after_each_exact",
    "after_sweep_gradient",
    "after_sweep_newton",
    "after_sweep_exact"
  ),
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────
# Global parameters
# ─────────────────────────────────────────────────────────────
lambda <- 0.1
tol <- 1e-8
max_iter <- 1e4
lr <- 1
# could be set to F if we were to test the effect, especially for gaussian fam
standardize <- TRUE

# ─────────────────────────────────────────────────────────────
# Metrics
# ─────────────────────────────────────────────────────────────
eval_metrics <- function(beta_true, beta_hat) {
  c(
    l1            = sum(abs(beta_hat - beta_true)),
    mse_full      = mean((beta_hat - beta_true)^2),
    intercept_err = abs(beta_hat[1] - beta_true[1])
  )
}

# ─────────────────────────────────────────────────────────────
# Simulate + fit once
# ─────────────────────────────────────────────────────────────
run_one <- function(sim_row, pcd_row) {
  sims <- switch(sim_row$design,
    independent = simulate_independent(
      n = sim_row$n, p = sim_row$p, k = sim_row$k,
      error_dist = sim_row$error_dist, family = sim_row$family,
      seed = sim_row$seed
    ),
    correlated = simulate_correlated(
      n = sim_row$n, p = sim_row$p, k = sim_row$k, rho = sim_row$rho,
      error_dist = sim_row$error_dist, family = sim_row$family,
      seed = sim_row$seed
    ),
    block = simulate_block_diagonal(
      n = sim_row$n, p = sim_row$p, k = sim_row$k, rho = sim_row$rho,
      block_sizes = floor(sim_row$p / 5),
      error_dist = sim_row$error_dist, family = sim_row$family,
      seed = sim_row$seed
    ),
    ar1 = simulate_ar1(
      n = sim_row$n, p = sim_row$p, k = sim_row$k, rho = sim_row$rho,
      burnin = 50, noise = "mnormal", sigma = diag(sim_row$p), df = 5,
      error_dist = sim_row$error_dist, family = sim_row$family,
      seed = sim_row$seed
    )
  )

  Xmat <- sims$X
  colnames(Xmat) <- paste0("X", seq_len(ncol(Xmat)))
  df <- data.frame(y = sims$y, Xmat)

  fit <- pcd(
    formula          = y ~ .,
    data             = df,
    lambda           = lambda,
    prox_fun         = pcd_row$prox_fun,
    intercept_update = pcd_row$intercept_update,
    family           = sim_row$family,
    method           = pcd_row$method,
    tol              = tol,
    max_iter         = max_iter,
    standardize      = standardize,
    lr               = lr,
    gamma            = pcd_row$gamma,
    verbose          = FALSE
  )

  # convergence diagnostics
  beta_path <- fit$path
  sweep_diff <- apply(abs(diff(beta_path)), 1, max)
  first_ok <- which(sweep_diff < tol)[1]
  conv_iter <- if (is.na(first_ok)) nrow(beta_path) else first_ok
  total_iter <- nrow(beta_path)
  n_updates <- fit$n_updates

  # metrics
  est_df <- fit$df
  beta_hat <- setNames(est_df$estimate, est_df$term)[
    c("(Intercept)", paste0("X", seq_len(sim_row$p)))
  ]
  metrics <- eval_metrics(sims$beta, beta_hat)

  # bundle results
  as.data.frame(
    c(
      as.list(sim_row), as.list(pcd_row),
      list(
        lambda = lambda,
        conv_iter = conv_iter,
        total_iter = total_iter,
        n_updates = n_updates
      ),
      metrics
    ),
    stringsAsFactors = FALSE
  )
}

# ─────────────────────────────────────────────────────────────
# Task grid
# ─────────────────────────────────────────────────────────────
runs_df <- expand.grid(
  sim_id = seq_len(nrow(sim_params)),
  pcd_id = seq_len(nrow(pcd_params)),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

total_runs <- nrow(runs_df)

# ─────────────────────────────────────────────────────────────
# Simulate + fit in parallel!
# ─────────────────────────────────────────────────────────────
run_one_safe <- function(idx, p) {
  library(pcd)
  library(Matrix)

  sp <- sim_params[runs_df$sim_id[idx], ]
  cfg <- pcd_params[runs_df$pcd_id[idx], ]

  res <- tryCatch(
    run_one(sp, cfg),
    error = function(e) {
      message(sprintf("run %d failed: %s", idx, e$message))
      NULL
    }
  )

  p(sprintf(
    "%s | %s | n=%d p=%d | %s | %s",
    sp$family, sp$design, sp$n, sp$p,
    cfg$method, cfg$intercept_update
  ))
  res
}

results_list <- with_progress({
  p <- progressor(steps = total_runs)

  future_lapply(
    seq_len(total_runs),
    run_one_safe,
    p               = p,
    future.seed     = TRUE,
    future.packages = c("pcd", "Matrix")
  )
})

results_df <- do.call(rbind, Filter(Negate(is.null), results_list))

# convert numeric columns
num_cols <- c(
  "n", "p", "k", "rho", "seed", "gamma", "lambda",
  "conv_iter", "total_iter",
  "l1", "mse_full", "intercept_err", "n_updates"
)
for (col in intersect(num_cols, names(results_df))) {
  results_df[[col]] <- as.numeric(results_df[[col]])
}

# ─────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────
write.csv(results_df, "results/full_sim_results.csv", row.names = FALSE)

# ─────────────────────────────────────────────────────────────
# Quick plotting, might delete later
# ─────────────────────────────────────────────────────────────

ggplot(
  results_df,
  aes(x = intercept_update, y = n_updates, fill = method)
) +
  geom_boxplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "n_updates until convergence",
    x = "Intercept-update strategy",
    y = "Numer of updates (log-scale)"
  )
