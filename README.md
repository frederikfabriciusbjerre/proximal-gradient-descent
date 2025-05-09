# pcd

This library was made for Project in Statistics @ UCPH. It is a lightweight R package providing **proximal coordinate-descent** solvers for penalized regression. Supports ordinary least squares (OLS), Lasso, and logistic regression (with L1 penalization), with flexible intercept-update schemes and choice of first-order (fixed-Lipschitz) or proximal-Newton methods (for logistic regression).

---

## Features

- **Generalized framework**: single `pcd()` wrapper for Gaussian and Binomial families  
- **Multiple penalties**: L1 (soft-threshold) or MCP (minimax concave)  
- **Intercept updates**:  
  - cyclic  
  - after each coordinate (gradient, Newton, or exact)  
  - after full sweep (gradient, Newton, or exact)  
- **Methods for logistic**:  
  - `"gradient"`: first-order proximal-gradient / CD with global Hessian bound  
  - `"newton"`: full proximal-Newton with working response and weights  
- **Tunable convergence**: user-specified `tol` and `max_iter`

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("frederikfabriciusbjerre/proximal-gradient-descent")
````

Or clone and build locally

---

## Usage

```r
library(pcd)
library(glmnet)
library(dplyr)

# Simulate a small logistic example
set.seed(1)
n <- 200 
p <- 5
beta0_true <- 1
beta_true  <- c(-0.5, 1.0, 0, -1.2, 0.8)
Xmat <- matrix(rnorm(n*p), n, p)
eta  <- beta0_true + Xmat %*% beta_true
y    <- rbinom(n, 1, 1/(1+exp(-eta)))
dat  <- as.data.frame(Xmat); names(dat) <- paste0("x",1:p)
dat$y <- y

# Lasso‐logistic via proximal‐gradient
fit_grad <- pcd(
  formula          = y ~ .,
  data             = dat,
  lambda           = 0.1,
  prox_fun         = "soft",
  standardize      = TRUE,
  intercept_update = "after_sweep_exact",
  family           = "binomial",
  method           = "gradient",
  tol              = 1e-12,
  max_iter         = Inf
)

# Lasso‐logistic via proximal‐Newton
fit_newt <- pcd(
  formula          = y ~ .,
  data             = dat,
  lambda           = 0.1,
  prox_fun         = "soft",
  standardize      = TRUE,
  intercept_update = "after_sweep_exact",
  family           = "binomial",
  method           = "newton",
  tol              = 1e-12,
  max_iter         = Inf
)

# Compare to glmnet
X_glm <- model.matrix(y ~ ., dat)[, -1]
glm_fit <- glmnet(X_glm, dat$y, family="binomial",
                  lambda=0.1, alpha=1, standardize=TRUE, thresh=1e-12)
glm_coefs <- as.vector(coef(glm_fit)); names(glm_coefs) <- rownames(coef(glm_fit))

tibble(
  term          = names(glm_coefs),
  grad_pcd      = fit_grad$estimate,
  newt_pcd      = fit_newt$estimate,
  glmnet        = glm_coefs
) |>
  mutate(
    diff_grad   = abs(grad_pcd - glmnet),
    diff_newton = abs(newt_pcd - glmnet)
  )
```
---

## License

This project is licensed under the **MIT License**. See [LICENSE.md](LICENSE) for details.

