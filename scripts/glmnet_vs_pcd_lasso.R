library(pcd)
library(glmnet)
library(dplyr)

set.seed(1)
lambda <- 0.10
tol <- 1e-15

dat <- mtcars |>
  rename(y = mpg)

# pcd
fit_pcd <- pcd(
  formula          = y ~ .,
  data             = dat,
  lambda           = lambda,
  prox_fun         = "soft",
  standardize      = TRUE,
  intercept_update = "after_each_exact",
  family           = "gaussian",
  tol              = tol,
  max_iter         = Inf,
  lr               = 0.1
)

# glm
X_glm <- model.matrix(y ~ ., dat)[, -1]

fit_glm <- glmnet(
  x           = X_glm,
  y           = dat$y,
  alpha       = 1, # lasso
  lambda      = lambda,
  standardize = TRUE,
  intercept   = TRUE,
  thresh      = tol
)

coef_glm <- as.vector(coef(fit_glm))
names(coef_glm) <- rownames(coef(fit_glm))

# compare
tibble(
  term             = names(coef_glm),
  estimate_pcd     = fit_pcd$df$estimate,
  estimate_glmnet  = coef_glm
) |>
  mutate(
    abs_diff = abs(estimate_pcd - estimate_glmnet)
  )
