# ────────────────────────────────────────────────────────────────────────
#  PUBLIC WRAPPERS
# ────────────────────────────────────────────────────────────────────────

#' Proximal coordinate descent for penalized regression
#'
#' This function implements proximal coordinate descent for penalized regression
#' using the specified proximal operator. It supports both linear regression,
#' lasso regression, and logistic regression with the option to use
#' either soft-thresholding or minimax concave penalty (MCP) as the proximal
#' operator.
#'
#' @param formula an object of class "formula"
#' (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment
#' containing the variables in the model.
#' @param lambda the penalty parameter (scalar). lambda >= 0.
#' @param prox_fun the proximal operator to use. Can be either "identity",
#' "soft", or "mcp". For "mcp", the `gamma` parameter must be specified.
#' @param intercept_update the strategy for updating the intercept term.
#' * cyclic: update the intercept term when the first predictor is visited.
#' * after_each_\*: update the intercept term after each slope update.
#' * after_sweep_gradient: update the intercept term after each sweep a single
#'   gradient descent step
#' * after_sweep_newton: update the intercept term after each sweep a single
#'   Newton's method step
#' * after_sweep_exact: update the intercept term after using Newtons method
#'   until convergence
#' @param family If "gaussian", the model is fitted using ordinary least squares.
#' If "binomial", the model is fitted using logistic regression.
#' @param method If family is "binomial", this argument specifies if
#' the method used for coordinate descent should be a first-order method
#' ("gradient), using an upper bound of the Hessian, or a proximal Newton method.
#' @param tol the convergence tolerance. The algorithm stops when the maximum
#' change in coefficients is less than this value.
#' @param max_iter the maximum number of iterations. If -1, the algorithm
#' will run until convergence.
#' @param standardize if TRUE, standardize the predictors before fitting the model.
#' @param lr Learning rate for the intercept update. This is only relevant if
#' `intercept_update` is set to `after_each_exact` or `after_sweep_exact`.
#' @param gamma the concavity parameter for MCP. This parameter is only
#' relevant if `prox_fun` is "mcp". It should be a positive number or `Inf`,
#' which indicates that the MCP is equivalent to soft-thresholding.
#' @param verbose if TRUE, print progress messages.
#' @returns data.frame(term, estimate)
#' @export
pcd <- function(
    formula, data,
    lambda,
    prox_fun = c("identity", "soft", "mcp"),
    intercept_update = c(
      "cyclic",
      "after_each_gradient",
      "after_each_newton",
      "after_each_exact",
      "after_sweep_gradient",
      "after_sweep_newton",
      "after_sweep_exact"
    ),
    family = c("gaussian", "binomial"),
    method = c("gradient", "newton"),
    tol = 1e-8,
    max_iter = 1e5,
    standardize = TRUE,
    lr = 1,
    gamma = -1, # MCP concavity parameter
    verbose = FALSE) {
  # checks
  if (missing(formula)) stop("formula is required")
  if (missing(data)) stop("data is required")
  if (missing(lambda)) stop("lambda is required")
  if (missing(prox_fun)) stop("prox_fun is required")
  stopifnot(is.numeric(lambda), length(lambda) == 1, lambda >= 0)
  if (max_iter < 1 && max_iter != -1) {
    stop("max_iter must be a positive integer or -1")
  }
  if (max_iter == -1L || is.infinite(max_iter)) {
    max_iter <- Inf
  } else {
    max_iter <- as.integer(max_iter)
  }

  # match arguments
  intercept_update <- match.arg(intercept_update)
  family <- match.arg(family)
  method <- match.arg(method)
  prox_fun <- switch(match.arg(prox_fun),
    identity = prox_identity,
    soft = prox_soft,
    mcp = if (gamma > 0) prox_mcp(gamma = gamma) else stop("gamma must be positive")
  )

  # get model frame from formula
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(attr(mf, "terms"), mf)

  # standardize data using mean and sd (population sd)
  if (standardize && ncol(X) > 1) {
    means <- colMeans(X[, -1, drop = FALSE])
    sds <- apply(X[, -1, drop = FALSE], 2, sd_population)
    zero_sd <- sds == 0 | !is.finite(sds)
    if (any(zero_sd)) {
      warning(sum(zero_sd), " zero-variance predictors dropped.")
      X <- X[, c(TRUE, !zero_sd), drop = FALSE]
      means <- means[!zero_sd]
      sds <- sds[!zero_sd]
    }
    X[, -1] <- sweep(X[, -1], 2, means, "-")
    X[, -1] <- sweep(X[, -1], 2, sds, "/")
  }

  # logistic regression
  if (family == "binomial") {
    if (missing(method)) stop("method is required")
    # call internal helper
    engine_out <- cd_engine_logistic(X, y,
      lambda = lambda,
      prox_fun = prox_fun,
      tol = tol,
      max_iter = max_iter,
      intercept_update = intercept_update,
      method = method,
      lr = lr,
      verbose = verbose
    )

    # linear regression
  } else if (family == "gaussian") {
    # call internal helper
    engine_out <- cd_engine(X, y,
      lambda = lambda,
      prox_fun = prox_fun,
      tol = tol,
      max_iter = max_iter,
      intercept_update = intercept_update,
      alpha = 1,
      lr = 1,
      verbose = verbose
    )
  } else {
    stop("Unsupported family for logistic regression.")
  }

  # extract beta and path
  beta_std <- engine_out$beta
  beta_path <- engine_out$beta_path
  n_updates <- engine_out$n_updates

  colnames(beta_path) <- colnames(X)

  # save return variables
  intercept <- beta_std[1]
  slopes <- beta_std[-1]

  # transform back to original scale
  if (standardize && length(slopes)) {
    intercept <- intercept - sum(slopes * means / sds)
    slopes <- slopes / sds
  }

  # return data frame
  out_df <- data.frame(
    term = colnames(X),
    estimate = c(intercept, slopes),
    row.names = NULL
  )

  structure(
    list(
      df = out_df,
      path = beta_path,
      terms = colnames(X),
      n_updates = n_updates
    ),
    class = "pcd_out"
  )
}


# ────────────────────────────────────────────────────────────────────────
#  INTERNAL HELPERS
# ────────────────────────────────────────────────────────────────────────


# ───────────────────────── coordinate descent helpers ─────────────────────────

#' Proximal coordinate descent for penalized regression
#'
#' @inheritParams pcd
#' @param X design matrix (n x p)
#' @param y response vector (n x 1)
#' @returns numeric vector of coefficients (p x 1), beta
cd_engine <- function(
    X, y,
    lambda = 0,
    prox_fun = prox_soft,
    tol = 1e-8,
    max_iter = 1e5,
    intercept_update = c(
      "cyclic",
      "after_each_gradient",
      "after_each_newton",
      "after_each_exact",
      "after_sweep_gradient",
      "after_sweep_newton",
      "after_sweep_exact"
    ),
    alpha = 1, # elastic‐net mixing parameter (α=1 ⇒ lasso)
    lr = 1,
    verbose = FALSE) {
  intercept_update <- match.arg(intercept_update)
  stopifnot(
    is.matrix(X), is.numeric(y), nrow(X) == length(y),
    is.function(prox_fun),
    alpha >= 0, alpha <= 1
  )

  # initialize
  n <- nrow(X)
  p <- ncol(X)
  has_intercept <- all(X[, 1] == 1)
  col_norm_squared <- colSums(X^2)
  beta <- numeric(p)
  r <- y # residual = y − Xβ
  iter <- 0L

  # initialize beta path for convergence plotting
  beta_path <- list()
  n_updates <- 0
  repeat {
    iter <- iter + 1L
    delta_max <- 0 # max change in coefficients

    # which coords to cycle over
    j_seq <- if (has_intercept && intercept_update != "cyclic") 2:p else 1:p

    # update sloaps
    for (j in j_seq) {
      x_j <- X[, j]
      # skip if colnorm squared is 0 or infinite
      if (!is.finite(col_norm_squared[j]) || col_norm_squared[j] == 0) next

      # coordinate update
      # denominater
      denom <- (col_norm_squared[j] / n) + (1 - alpha) * lambda

      # z is the unpenalized least‐squares coordinate update
      z <- (sum(x_j * r) / n +
        (col_norm_squared[j] / n) * beta[j]
      )

      # update beta using the proximal operator
      beta_new <- prox_fun(z, alpha * lambda * (j > 1)) / denom
      n_updates <- n_updates + 1L

      # save the change
      change <- beta_new - beta[j]
      if (!is.finite(change) || change == 0) next

      # Update the residual: since
      #   r = z - (beta0 + X %*% beta),
      # changing beta[j] by 'change' gives
      #   r_i_new = z_i - (beta0 + sum_k X[i,k] * beta_k_new)
      #           = [z_i - (beta0 + sum_k X[i,k] * beta_k_old)] - X[i,j] * change
      #           = r_i_old - X[i,j] * change
      r <- r - change * x_j
      beta[j] <- beta_new

      # track max change for convergence
      delta_max <- max(delta_max, abs(change))

      # intercept after each
      if (has_intercept && startsWith(intercept_update, "after_each") && j > 1) {
        # gradient and hessian for intercept update
        grad_0 <- -sum(r)
        H_0 <- n
        delta0 <- switch(intercept_update,
          after_each_gradient = -grad_0 / n,
          after_each_newton = -grad_0 / H_0,
          # calculate exact intercept update with Newton's method
          after_each_exact = {
            newton_tol <- 1e-12
            newton_max <- 50
            acc <- 0
            for (k in seq_len(newton_max)) {
              n_updates <- n_updates + 1L
              g <- -sum(r)
              if (abs(g) < newton_tol) break
              step0 <- -g / H_0
              r <- r - step0
              acc <- acc + step0
            }
            acc
          }
        )
        # update intercept
        if (delta0 != 0) {
          n_updates <- n_updates + 1L
          beta[1] <- beta[1] + delta0
          r <- r - delta0
          delta_max <- max(delta_max, abs(delta0))
        }
      }
    }

    # intercept after sweep
    if (has_intercept && startsWith(intercept_update, "after_sweep")) {
      # recompute eta and pi
      eta <- as.numeric(X %*% beta)
      pi <- 1 / (1 + exp(-eta))

      # gradient and hessian for intercept update
      grad_0 <- -sum(r)
      H_0 <- n
      delta0 <- switch(intercept_update,
        after_sweep_gradient = -grad_0 / n,
        after_sweep_newton = -grad_0 / H_0,
        # calculate exact intercept update with Newton's method
        after_sweep_exact = {
          newton_tol <- 1e-12
          newton_max <- 50
          acc <- 0
          for (k in seq_len(newton_max)) {
            n_updates <- n_updates + 1L
            g <- -sum(r)
            if (abs(g) < newton_tol) break
            step0 <- -g / H_0
            r <- r - step0
            acc <- acc + step0
          }
          acc
        }
      )
      # update intercept
      if (delta0 != 0) {
        n_updates <- n_updates + 1L
        beta[1] <- beta[1] + delta0
        if (intercept_update != "after_sweep_exact") {
          r <- r - delta0
        }
        delta_max <- max(delta_max, abs(delta0))
      }
    }

    if (verbose && iter %% 1000L == 0L) {
      message(
        "iter=", iter,
        "  max|Δβ|=", formatC(delta_max, digits = 3, format = "e")
      )
    }

    # record full beta vector for convergence plotting
    beta_path[[iter]] <- beta

    if (delta_max < tol || iter >= max_iter) break
  }

  # convert beta path to matrix
  beta_path <- do.call(rbind, beta_path)

  return(list(beta_path = beta_path, beta = beta, n_updates))
}

#' Logistic regression via coordinate descent, logistic loss
#' @inheritParams pcd
#' @param X design matrix (n x p)
#' @param y response vector (n x 1)
#' @returns numeric vector of coefficients (p x 1), beta
cd_engine_logistic <- function(X, y,
                               lambda,
                               standardize = TRUE,
                               intercept_update = c(
                                 "cyclic",
                                 "after_each_gradient",
                                 "after_each_newton",
                                 "after_each_exact",
                                 "after_sweep_gradient",
                                 "after_sweep_newton",
                                 "after_sweep_exact"
                               ),
                               method = c("gradient", "newton"),
                               prox_fun = prox_soft,
                               tol = 1e-8,
                               max_iter = Inf,
                               lr = 1,
                               verbose = FALSE) {
  # checks
  intercept_update <- match.arg(intercept_update)
  method <- match.arg(method)
  stopifnot(is.numeric(lambda), length(lambda) == 1, lambda >= 0)
  stopifnot(is.numeric(tol), length(tol) == 1, tol > 0)
  stopifnot(
    is.numeric(max_iter), length(max_iter) == 1,
    (max_iter >= 1 && max_iter == floor(max_iter)) ||
      is.infinite(max_iter)
  )

  # initialize
  n <- nrow(X)
  p <- ncol(X)
  has_intercept <- all(X[, 1] == 1)
  beta <- numeric(p)
  iter <- 0L

  # initialize beta path for convergence plotting
  beta_path <- list()
  n_updates <- 0L

  # first‐order proximal gradient method
  if (method == "gradient") {
    # Lipschitz constants for each coordinate
    L_j <- colSums(X^2) / (4 * n)

    repeat {
      iter <- iter + 1L
      delta_max <- 0 # max change in coefficients

      # which coords to cycle over
      j_seq <- if (has_intercept && intercept_update != "cyclic") 2:p else 1:p

      # do one sweep of prox‐grad updates
      eta <- X %*% beta
      pi <- 1 / (1 + exp(-eta))

      for (j in j_seq) {
        x_j <- X[, j]
        # gradient
        grad_j <- sum((pi - y) * x_j) / n
        t_j <- 1 / L_j[j]

        # un‐penalized step + soft‐threshold
        z <- beta[j] - t_j * grad_j
        beta_new <- prox_fun(z, lambda * t_j * (j > 1))

        n_updates <- n_updates + 1L

        change <- beta_new - beta[j]
        if (!is.finite(change) || change == 0) next

        beta[j] <- beta_new
        delta_max <- max(delta_max, abs(change))

        # update eta and pi incrementally
        eta <- eta + change * x_j
        pi <- 1 / (1 + exp(-eta))

        # intercept after each slope update
        if (has_intercept && startsWith(intercept_update, "after_each") && j > 1) {
          # gradient and hessian for intercept update
          grad_0 <- sum(pi - y)
          H_0 <- sum(pi * (1 - pi))
          delta0 <- switch(intercept_update,
            after_each_gradient = -4 * grad_0 / n,
            after_each_newton = -grad_0 / H_0,
            after_each_exact = {
              grad_tol <- 1e-12
              max_inner <- 1000 # safety cap
              acc <- 0

              for (k in seq_len(max_inner)) {
                n_updates <- n_updates + 1L
                grad_0 <- 4 * sum(pi - y) / n # average gradient wrt β0
                if (abs(grad_0) < grad_tol) break

                step0 <- -lr * grad_0
                eta <- eta + step0
                pi <- 1 / (1 + exp(-eta))

                acc <- acc + step0
              }
              acc
            }
          )
          # update intercept
          if (delta0 != 0) {
            n_updates <- n_updates + 1L
            beta[1] <- beta[1] + delta0
            if (intercept_update != "after_each_exact") {
              eta <- eta + delta0
              pi <- 1 / (1 + exp(-eta))
            }
            delta_max <- max(delta_max, abs(delta0))
            w <- pi * (1 - pi)
            z <- eta + (y - pi) / w
            r <- z - eta
          }
        }
      }

      # intercept after sweep
      if (has_intercept && startsWith(intercept_update, "after_sweep")) {
        # recompute eta and pi
        eta <- as.numeric(X %*% beta)
        pi <- 1 / (1 + exp(-eta))

        # gradient and hessian for intercept update
        grad_0 <- sum(pi - y)
        H_0 <- sum(pi * (1 - pi))
        delta0 <- switch(intercept_update,
          after_sweep_gradient = -4 * grad_0 / n,
          after_sweep_newton = -grad_0 / H_0,
          after_sweep_exact = {
            grad_tol <- 1e-12
            max_inner <- 1000 # safety cap
            acc <- 0

            for (k in seq_len(max_inner)) {
              n_updates <- n_updates + 1L
              grad_0 <- 4 * sum(pi - y) / n # average gradient wrt β0
              if (abs(grad_0) < grad_tol) break

              step0 <- -lr * grad_0
              eta <- eta + step0
              pi <- 1 / (1 + exp(-eta))

              acc <- acc + step0
            }
            acc
          }
        )
        # update intercept
        if (delta0 != 0) {
          n_updates <- n_updates + 1L
          beta[1] <- beta[1] + delta0
          if (intercept_update != "after_sweep_exact") {
            eta <- eta + delta0
            pi <- 1 / (1 + exp(-eta))
          }
          delta_max <- max(delta_max, abs(delta0))
        }
      }

      if (verbose && iter %% 1e3 == 0L) {
        message("iter=", iter, "  max|Δβ|=", formatC(delta_max, digits = 3, format = "e"))
      }

      # record full beta vector for convergence plotting
      beta_path[[iter]] <- beta
      if (delta_max < tol || iter >= max_iter) break
    }
  }

  # proximal Newton coordinate descent
  else if (method == "newton") {
    repeat {
      iter <- iter + 1L
      delta_max <- 0

      # do one sweep of prox‐Newton updates
      eta <- as.numeric(X %*% beta)
      pi <- 1 / (1 + exp(-eta))
      w <- pi * (1 - pi)
      z <- eta + (y - pi) / w

      # incremental residual for weighted LS
      r <- z - eta

      # same j_seq logic
      j_seq <- if (has_intercept && intercept_update != "cyclic") 2:p else 1:p

      # main loop
      for (j in j_seq) {
        x_j <- X[, j]
        # partial residual
        r_j <- r + x_j * beta[j]

        # weighted‐LS grad and hessian
        grad_j <- sum(w * x_j * r_j) / n
        H_j <- sum(w * x_j^2) / n

        # update
        beta_new <- prox_fun(grad_j, lambda * (j > 1)) / H_j
        change <- beta_new - beta[j]
        if (!is.finite(change) || change == 0) next

        beta[j] <- beta_new

        n_updates <- n_updates + 1L

        delta_max <- max(delta_max, abs(change))

        # update residual
        r <- r - x_j * change

        # update eta and pi incrementally
        eta <- eta + change * x_j
        pi <- 1 / (1 + exp(-eta))

        # intercept after each update
        if (has_intercept && startsWith(intercept_update, "after_each") && j > 1) {
          # gradient and hessian for intercept update
          grad_0 <- sum(pi - y)
          H_0 <- sum(pi * (1 - pi))
          delta0 <- switch(intercept_update,
            after_each_gradient = -4 * grad_0 / n,
            after_each_newton = -grad_0 / H_0,
            after_each_exact = {
              grad_tol <- 1e-12
              max_inner <- 1000 # safety cap
              acc <- 0

              for (k in seq_len(max_inner)) {
                n_updates <- n_updates + 1L
                grad_0 <- 4 * sum(pi - y) / n
                if (abs(grad_0) < grad_tol) break

                step0 <- -lr * grad_0
                eta <- eta + step0
                pi <- 1 / (1 + exp(-eta))

                acc <- acc + step0
              }
              acc
            }
          )
          # update intercept
          if (delta0 != 0) {
            n_updates <- n_updates + 1L
            beta[1] <- beta[1] + delta0
            if (intercept_update != "after_each_exact") {
              eta <- eta + delta0
              pi <- 1 / (1 + exp(-eta))
            }
            delta_max <- max(delta_max, abs(delta0))
            w <- pi * (1 - pi)
            z <- eta + (y - pi) / w
            r <- z - eta
          }
        }
      }

      # intercept after sweep
      if (has_intercept && startsWith(intercept_update, "after_sweep")) {
        # recompute eta and pi
        eta <- as.numeric(X %*% beta)
        pi <- 1 / (1 + exp(-eta))

        # gradient and hessian for intercept update
        grad_0 <- sum(pi - y)
        H_0 <- sum(pi * (1 - pi))
        delta0 <- switch(intercept_update,
          after_sweep_gradient = -4 * grad_0 / n,
          after_sweep_newton = -sum(pi - y) / H_0,
          after_sweep_exact = {
            grad_tol <- 1e-12
            max_inner <- 1000 # safety cap
            acc <- 0

            for (k in seq_len(max_inner)) {
              n_updates <- n_updates + 1L
              grad_0 <- 4 * sum(pi - y) / n
              if (abs(grad_0) < grad_tol) break

              step0 <- -lr * grad_0
              eta <- eta + step0
              pi <- 1 / (1 + exp(-eta))

              acc <- acc + step0
            }
            acc
          }
        )
        if (delta0 != 0) {
          n_updates <- n_updates + 1L
          beta[1] <- beta[1] + delta0
          delta_max <- max(delta_max, abs(delta0))
          if (intercept_update != "after_sweep_exact") {
            eta <- eta + delta0
            pi <- 1 / (1 + exp(-eta))
          }
        }
      }

      if (verbose && iter %% 10L == 0L) {
        message("iter=", iter, "  max|Δβ|=", formatC(delta_max, format = "e", digits = 2))
      }

      # record full beta vector for convergence plotting
      beta_path[[iter]] <- beta

      if (delta_max < tol || iter >= max_iter) break
    }
  }

  # convert beta path to matrix
  beta_path <- do.call(rbind, beta_path)

  return(list(beta_path = beta_path, beta = beta, n_updates = n_updates))
}


# ───────────────────────── proximal function helpers ──────────────────────────

#' Soft-threshold proximal operator (ℓ₁ penalty)
#' @param z vector of coefficients
#' @param lambda penalty parameter. If lambda = 0, then this is the identity
#' operator.
#' @keywords internal
prox_soft <- function(z, lambda) sign(z) * pmax(abs(z) - lambda, 0)

#' Minimax concave penalty proximal operator (MCP)
#' @param z vector of coefficients
#' @param lambda penalty parameter
#' @param gamma concavity parameter (default = Inf)
#' If gamma = Inf then MCP = soft-thresholding
#' @keywords internal
prox_mcp <- function(gamma) {
  function(z, lambda) {
    thresh1 <- lambda
    thresh2 <- gamma * lambda
    absz <- abs(z)

    # initialize zero vector
    out <- numeric(length(z))

    # initialize index masks
    idx1 <- absz <= thresh1
    idx2 <- absz > thresh1 & absz <= thresh2
    idx3 <- absz > thresh2

    # assign non-zero values to index masks in a vectorized way
    out[idx2] <- sign(z[idx2]) * ((absz[idx2] - thresh1) / (1 - 1 / gamma))
    out[idx3] <- z[idx3]
    out
  }
}

# ─────────────────────────────────── MISC ─────────────────────────────────────

#' Population standard deviation
#'
#' @keywords internal
sd_population <- function(x) {
  n <- length(x)
  sqrt((n - 1) / n) * sd(x)
}
