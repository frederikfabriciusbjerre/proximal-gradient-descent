library(MASS) # for mvrnorm
library(freqdom) # for rar()

#' Simulate independent design
#'
#' @param n Number of observations
#' @param p Number of features
#' @param k Number of non-zero coefficients
#' @param error_dist Error distribution: "normal", "t", "uniform"
#' @param beta Optional true coefficient vector
#' @param seed Random seed
#' @return List with X, y, beta
#' @export
simulate_independent <- function(n, p, k,
                                 error_dist = c("normal", "t", "uniform"),
                                 beta = NULL,
                                 seed = 96) {
  set.seed(seed)
  # covariance matrix
  Sigma <- diag(p)

  # model matrix
  X <- MASS::mvrnorm(n, mu = runif(p, -2, 2), Sigma = Sigma)

  # init beta if not provided
  if (is.null(beta)) {
    beta <- numeric(p)
    nz <- sample.int(p, k)
    beta[nz] <- rnorm(k)
    beta[1] <- runif(1, -2, 2) # make sure we have non-zero intercept
  }

  # generate errors
  error_dist <- match.arg(error_dist)
  eps <- switch(error_dist,
    "normal" = rnorm(n, 0, 1),
    "t" = rt(n, df = 5),
    "uniform" = runif(n, -3, 3)
  )

  # response
  y <- X %*% beta + eps

  # return
  list(X = X, y = as.numeric(y), beta = beta)
}

#' Simulate correlated design
#'
#' @inheritParams simulate_independent
#' @param rho Correlation coefficient used to simulate from uniform 0 to rho
#' @export
simulate_correlated <- function(n, p, k,
                                rho = 0,
                                error_dist = c("normal", "t", "uniform"),
                                beta = NULL,
                                seed = 96) {
  set.seed(seed)

  # correlation matrix
  R <- runif(p^2, 0, rho) |> matrix(p, p)
  diag(R) <- 1

  # covariance matrix
  S <- rep(1, p) # sd is set to 1 for all variables
  Sigma <- outer(S, S) * R

  # model matrix
  X <- MASS::mvrnorm(n, mu = runif(p, -2, 2), Sigma = Sigma)

  # init beta if not provided
  if (is.null(beta)) {
    beta <- numeric(p)
    nz <- sample.int(p, k)
    beta[nz] <- rnorm(k)
    beta[1] <- runif(1, -2, 2) # make sure we have non-zero intercept
  }

  # generate errors
  error_dist <- match.arg(error_dist)
  eps <- switch(error_dist,
    "normal" = rnorm(n, 0, 1),
    "t" = rt(n, df = 5),
    "uniform" = runif(n, -3, 3)
  )

  # response
  y <- X %*% beta + eps

  # returns
  list(X = X, y = y, beta = beta)
}

#' Simulate block-diagonal design
#'
#' @inheritParams simulate_independent
#' @param rho Correlation coefficient for block matrix (non-diagonals)
#' @param block_sizes Sizes of blocks in the block-diagonal matrix
#' @export
simulate_block_diagonal <- function(n, p, k,
                                    rho = 0,
                                    block_sizes = NULL,
                                    error_dist = c("normal", "t", "uniform"),
                                    beta = NULL,
                                    seed = 96) {
  set.seed(seed)
  if (is.null(block_sizes)) stop("block_sizes must be provided")
  if (rho < 0 || rho > 1) stop("rho must be in [0,1]")

  # block correlation matrix
  R <- Matrix::kronecker(
    diag(k),
    matrix(rho, nrow = block_sizes, ncol = block_sizes)
  )
  diag(R) <- 1

  # covariance matrix
  S <- rep(1, p) # sd is set to 1 for all variables
  Sigma <- outer(S, S) * R

  # model matrix
  X <- MASS::mvrnorm(n, mu = runif(p, 0, rho), Sigma = Sigma)

  # initalize beta if not provided
  if (is.null(beta)) {
    beta <- numeric(p)
    nz <- sample.int(p, k)
    beta[nz] <- rnorm(k)
  }

  # generate errors
  error_dist <- match.arg(error_dist)
  eps <- switch(error_dist,
    "normal" = rnorm(n, 0, 1),
    "t" = rt(n, df = 5),
    "uniform" = runif(n, -3, 3)
  )

  # response
  y <- X %*% beta + eps

  # returns
  list(X = X, y = y, beta = beta)
}

#' Simulate VAR(1) / multivariate AR(1) design
#'
#' @inheritParams simulate_independent
#' @param rho Correlation coefficient for AR(1) process
#' @param burnin Number of initial observations to discard
#' @param noise Type of noise to use for AR process: "mnormal" or "mt"
#' @param sigma Covariance matrix for the noise
#' @param df Degrees of freedom for t-distribution for noise
#' @param Psi AR coefficient matrix
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
                         seed = 96) {
  set.seed(seed)
  noise <- match.arg(noise)
  if (is.null(sigma)) sigma <- diag(p)
  if (is.null(Psi)) {
    Psi <- array(0, dim = c(p, p, 1))
    Psi[, , 1] <- diag(p) * rho
  }

  # generate model matrix using autoregressive process
  X <- freqdom::rar(
    n = n,
    d = p,
    Psi = Psi,
    burnin = burnin,
    noise = noise,
    sigma = sigma,
    df = df
  )

  # initialize beta if not provided
  if (is.null(beta)) {
    beta <- numeric(p)
    nz <- sample.int(p, k)
    beta[nz] <- rnorm(k)
  }
  # generate errors
  error_dist <- match.arg(error_dist)
  eps <- switch(error_dist,
    "normal" = rnorm(n, 0, 1),
    "t" = rt(n, df = 5),
    "uniform" = runif(n, -3, 3)
  )

  # response
  y <- X %*% beta + eps

  # returns
  list(X = X, y = y, beta = beta)
}
