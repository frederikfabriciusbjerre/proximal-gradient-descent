library(pcd)
library(glmnet)
library(tidyverse)

# params
n <- 50
p <- 300
beta0_true <- 1.5
beta_true <- runif(p, -1, 1) * c(1, 0, 1, 0, 0, 0, 1, 0, 1, 1)
lambda <- 0.1
nsim <- 10
tol <- 1e-12
strategies <- c(
  "cyclic",
  "after_each_gradient",
  "after_each_newton",
  "after_each_exact",
  "after_sweep_gradient",
  "after_sweep_newton",
  "after_sweep_exact"
)

# single sim rep
simulate_once <- function(seed) {
  set.seed(seed)
  Xmat <- matrix(rnorm(n * p), n, p)
  eta <- beta0_true + Xmat %*% beta_true
  y <- rbinom(n, 1, plogis(eta))

  # prep
  dat <- as_tibble(Xmat)
  names(dat) <- paste0("x", seq_len(p))
  dat$y <- y
  form <- as.formula(paste("y ~", paste0("x", seq_len(p), collapse = " + ")))

  # glmnet baseline
  tm_glm <- system.time({
    fit_glm <- glmnet(
      x           = Xmat,
      y           = y,
      family      = "binomial",
      lambda      = lambda,
      standardize = TRUE,
      intercept   = TRUE
    )
  })
  coefs_glm <- as.vector(coef(fit_glm))
  errs_glm <- coefs_glm - c(beta0_true, beta_true)
  glm_row <- tibble(
    strategy       = "glmnet",
    abs_err_beta0  = abs(coefs_glm[1] - beta0_true),
    abs_err_slopes = mean(abs(coefs_glm[-1] - beta_true)),
    abs_err_full   = mean(abs(errs_glm)),
    sq_err_full    = mean(errs_glm^2),
    time_sec       = as.numeric(tm_glm["elapsed"])
  )

  # pcd
  pcd_rows <- purrr::map_df(strategies, function(s) {
    tm_pcd <- system.time({
      fit_pcd <- pcd(
        formula          = form,
        data             = dat,
        lambda           = lambda,
        prox_fun         = "soft",
        standardize      = TRUE,
        intercept_update = s,
        family           = "binomial",
        method           = "gradient",
        tol              = tol,
        max_iter         = Inf
      )
    })

    # extract estimates
    beta_hat <- fit_pcd$df$estimate
    errs <- beta_hat - c(beta0_true, beta_true)
    tibble(
      strategy       = s,
      abs_err_beta0  = abs(beta_hat[1] - beta0_true),
      abs_err_slopes = mean(abs(beta_hat[-1] - beta_true)),
      abs_err_full   = mean(abs(errs)),
      sq_err_full    = mean(errs^2),
      time_sec       = as.numeric(tm_pcd["elapsed"])
    )
  })

  bind_rows(glm_row, pcd_rows)
}

# simulate
sim_res <- map_df(seq_len(nsim), simulate_once)

# summarize
summary_res <- sim_res |>
  group_by(strategy) |>
  summarize(
    mean_abs_err_beta0 = mean(abs_err_beta0),
    mean_abs_err_slopes = mean(abs_err_slopes),
    mean_MAE_full = mean(abs_err_full),
    mean_MSE_full = mean(sq_err_full),
    mean_time_sec = mean(time_sec),
    .groups = "drop"
  ) |>
  arrange(mean_MAE_full)

print(summary_res)
