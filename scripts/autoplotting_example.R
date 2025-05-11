library(pcd)
library(glmnet)
library(dplyr)

set.seed(1)
n <- 200
p <- 5
β0_true <- 1
β_true <- c(-0.5, 1.0, 0, -1.2, 0.8)
Xmat <- matrix(rnorm(n * p), n, p)
eta <- β0_true + Xmat %*% β_true
prob <- 1 / (1 + exp(-eta))
y <- rbinom(n, 1, prob)

dat <- as.data.frame(Xmat)
colnames(dat) <- paste0("x", seq_len(p))
dat$y <- y

tol <- 1e-12

# fit with pcd
# using first order
fit_pcd <- pcd(
  formula = y ~ .,
  data = dat,
  lambda = 0.1,
  prox_fun = "soft",
  standardize = TRUE,
  intercept_update = "after_sweep_newton",
  family = "binomial",
  method = "gradient",
  tol = tol,
  max_iter = Inf
)

p <- autoplot(fit_pcd, tol = tol)

# using newton
fit_pcd_newton <- pcd(
  formula = y ~ .,
  data = dat,
  lambda = 0.1,
  prox_fun = "soft",
  standardize = TRUE,
  intercept_update = "after_sweep_newton",
  family = "binomial",
  method = "newton",
  tol = tol,
  max_iter = Inf
)

q <- autoplot(fit_pcd_newton, tol = tol)

# plot side by side
library(gridExtra)
library(ggplot2)
grid.arrange(p, q, ncol = 2)

# compare
tibble(
  term = fit_pcd$terms,
  estimate_pcd = fit_pcd$df$estimate,
  estimate_pcd_newton = fit_pcd_newton$df$estimate
) |>
  # mutate(abs_diff_pcd_grad_glm = abs(estimate_pcd - estimate_glmnet)) |>
  # mutate(abs_diff_pcd_newton_glm = abs(estimate_pcd_newton - estimate_glmnet)) |>
  mutate(abs_diff_pcd_grad_newton = abs(estimate_pcd - estimate_pcd_newton)) |>
  print()
