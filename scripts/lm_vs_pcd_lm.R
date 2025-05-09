library(pcd)
library(broom) # for tidy() on lm()

set.seed(1)
tol <- 1e-15

dat <- mtcars |>
  rename(y = mpg)

# pcd
fit_pcd <- pcd(
  formula          = y ~ .,
  data             = dat,
  lambda           = 0, # no penalty
  prox_fun         = "soft", # soft‐threshold with λ=0 is identity
  standardize      = TRUE,
  intercept_update = "cyclic", # any valid intercept scheme
  family           = "gaussian",
  tol              = 1e-12,
  max_iter         = Inf
)

# lm
fit_lm <- tidy(lm(mpg ~ ., data = mtcars))[, c("term", "estimate")]


# compare
comparison <- merge(
  fit_pcd, fit_lm,
  by = "term",
  suffixes = c("_pcd", "_lm")
)
comparison$abs_diff <- abs(comparison$estimate_pcd - comparison$estimate_lm)

print(comparison, row.names = FALSE, digits = 6)
