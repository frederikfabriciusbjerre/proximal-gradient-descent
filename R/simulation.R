library(MASS) # for mvrnorm
library(freqdom) # for rar()
library(Matrix) # for nearPD, kronecker

#' Simulate independent design
#'
#' @param n Number of observations
#' @param p Number of features
#' @param k Number of non-zero coefficients
#' @param error_dist Error distribution: "normal", "t", "uniform"
#' @param beta Optional true coefficient vector (length p+1: intercept + slopes)
#' @param family "gaussian" or "binomial"
#' @param seed Random seed
#' @return List with X (n×p), y, beta (length p+1)
#' @export
simulate_independent <- function(n, p, k,
                                 error_dist = c("normal", "t", "uniform"),
                                 beta = NULL,
                                 family = c("gaussian", "binomial"),
                                 seed = 96) {
  set.seed(seed)

  # 1) raw features (n x p)
  Sigma <- diag(p)
  X <- MASS::mvrnorm(n, mu = runif(p, -2, 2), Sigma = Sigma)

  # 2) true beta = (intercept, slopes)
  if (is.null(beta)) {
    beta <- numeric(p + 1)
    nz <- sample.int(p, k)
    beta[nz + 1] <- rnorm(k)
    beta[1] <- runif(1, -2, 2)
  }

  # 3) error
  error_dist <- match.arg(error_dist)
  eps <- switch(error_dist,
    normal  = rnorm(n),
    t       = rt(n, df = 5),
    uniform = runif(n, -3, 3)
  )

  # 4) linear predictor & response
  eta <- as.numeric(beta[1] + X %*% beta[-1] + eps)
  if (family[1] == "binomial") {
    prob <- 1 / (1 + exp(-eta))
    y <- rbinom(n, 1, prob)
  } else if (family[1] == "gaussian") {
    y <- eta
  } else {
    stop("`family` must be 'gaussian' or 'binomial'")
  }

  list(X = X, y = y, beta = beta)
}


#' Simulate correlated design
#'
#' @inheritParams simulate_independent
#' @param rho Correlation coefficient (off‐diagonals uniform on [0, rho])
#' @export
simulate_correlated <- function(n, p, k,
                                rho = 0,
                                error_dist = c("normal", "t", "uniform"),
                                beta = NULL,
                                family = c("gaussian", "binomial"),
                                seed = 96) {
  set.seed(seed)

  # 1) build covariance
  R <- matrix(runif(p^2, 0, rho), p, p)
  diag(R) <- 1
  Sigma <- nearPD(R, keepDiag = TRUE, ensureSymmetry = TRUE)$mat

  # 2) raw features
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

  # 3) true beta
  if (is.null(beta)) {
    beta <- numeric(p + 1)
    nz <- sample.int(p, k)
    beta[nz + 1] <- rnorm(k)
    beta[1] <- runif(1, -2, 2)
  }

  # 4) error
  error_dist <- match.arg(error_dist)
  eps <- switch(error_dist,
    normal  = rnorm(n),
    t       = rt(n, df = 5),
    uniform = runif(n, -3, 3)
  )

  # 5) response
  eta <- as.numeric(beta[1] + X %*% beta[-1] + eps)
  if (family[1] == "binomial") {
    prob <- 1 / (1 + exp(-eta))
    y <- rbinom(n, 1, prob)
  } else if (family[1] == "gaussian") {
    y <- eta
  } else {
    stop("`family` must be 'gaussian' or 'binomial'")
  }

  list(X = X, y = y, beta = beta)
}


#' Simulate block‐diagonal design
#'
#' @inheritParams simulate_independent
#' @param rho Correlation coefficient for within‐block off‐diagonals
#' @param block_sizes Either
#'   * a single integer B, meaning “make B equally‐sized blocks” (so p must be divisible by B), or
#'   * a numeric vector of length B whose entries sum to p, giving each block’s size.
#' @export
simulate_block_diagonal <- function(n, p, k,
                                    rho = 0,
                                    block_sizes = NULL,
                                    error_dist = c("normal", "t", "uniform"),
                                    beta = NULL,
                                    family = c("gaussian", "binomial"),
                                    seed = 96) {
  set.seed(seed)
  if (is.null(block_sizes)) {
    stop("`block_sizes` must be provided")
  }

  # interpret block_sizes
  if (length(block_sizes) == 1L) {
    B <- block_sizes
    if (p %% B != 0L) {
      stop("When block_sizes is a single integer B, p must be divisible by B")
    }
    block_sizes_vec <- rep(p / B, B)
  } else {
    block_sizes_vec <- block_sizes
    if (sum(block_sizes_vec) != p) {
      stop("When block_sizes is a vector, its entries must sum to p")
    }
  }

  # build block covariance
  mats <- lapply(block_sizes_vec, function(bs) matrix(rho, nrow = bs, ncol = bs))
  Rblock <- Matrix::bdiag(mats)
  diag(Rblock) <- 1
  Sigma <- nearPD(Rblock, keepDiag = TRUE, ensureSymmetry = TRUE)$mat

  # raw features
  X <- MASS::mvrnorm(n, mu = runif(p, 0, rho), Sigma = Sigma)

  # true beta = (intercept, slopes)
  if (is.null(beta)) {
    beta <- numeric(p + 1)
    nz <- sample.int(p, k)
    beta[nz + 1] <- rnorm(k)
    beta[1] <- runif(1, -2, 2)
  }

  # errors
  error_dist <- match.arg(error_dist)
  eps <- switch(error_dist,
    normal  = rnorm(n),
    t       = rt(n, df = 5),
    uniform = runif(n, -3, 3)
  )

  # linear predictor + response
  eta <- as.numeric(beta[1] + X %*% beta[-1] + eps)
  if (family[1] == "binomial") {
    prob <- 1 / (1 + exp(-eta))
    y <- rbinom(n, 1, prob)
  } else if (family[1] == "gaussian") {
    y <- eta
  } else {
    stop("`family` must be 'gaussian' or 'binomial'")
  }

  list(X = X, y = y, beta = beta)
}


#' Simulate VAR(1) / multivariate AR(1) design
#'
#' @inheritParams simulate_independent
#' @param rho AR(1) coefficient
#' @param burnin Number of initial draws to discard
#' @param noise "mnormal" or "mt"
#' @param sigma Covariance for the noise (p×p)
#' @param df Degrees of freedom for t‐noise
#' @param Psi AR coefficient array (p×p×lag)
#' @export
simulate_ar1 <- function(n, p, k,
                         rho = 0,
                         burnin = 50,
                         noise = c("mnormal", "mt"),
                         sigma = NULL,
                         df = 5,
                         Psi = NULL,
                         error_dist = c("normal", "t", "uniform"),
                         beta = NULL,
                         family = c("gaussian", "binomial"),
                         seed = 96) {
  set.seed(seed)
  noise <- match.arg(noise)
  if (is.null(sigma)) sigma <- diag(p)
  if (is.null(Psi)) Psi <- array(diag(p) * rho, dim = c(p, p, 1))

  # 1) raw AR(1) features
  X <- freqdom::rar(
    n = n, d = p,
    Psi = Psi, burnin = burnin,
    noise = noise, sigma = sigma, df = df
  )

  # 2) true beta
  if (is.null(beta)) {
    beta <- numeric(p + 1)
    nz <- sample.int(p, k)
    beta[nz + 1] <- rnorm(k)
    beta[1] <- runif(1, -2, 2)
  }

  # 3) error
  error_dist <- match.arg(error_dist)
  eps <- switch(error_dist,
    normal  = rnorm(n),
    t       = rt(n, df = 5),
    uniform = runif(n, -3, 3)
  )

  # 4) response
  eta <- as.numeric(beta[1] + X %*% beta[-1] + eps)
  if (family[1] == "binomial") {
    prob <- 1 / (1 + exp(-eta))
    y <- rbinom(n, 1, prob)
  } else if (family[1] == "gaussian") {
    y <- eta
  } else {
    stop("`family` must be 'gaussian' or 'binomial'")
  }

  list(X = X, y = y, beta = beta)
}
