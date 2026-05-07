#' @importFrom survival Surv
#' @export
survival::Surv

#' Generic function for basehaz
#'
#' @param object Any object.
#' @param ... Additional arguments.
#'
#' @return A numeric vector of baseline hazard.
#' @export
basehaz <- function(object, ...) {
  UseMethod("basehaz")
}

#' Generic function for diagnose
#'
#' @param object Any object.
#' @param ... Additional arguments.
#' @return Called for its side effect of producing diagnostic plots.
#' Returns \code{NULL} invisibly.
#' @export
diagnose <- function(object, ...) {
  UseMethod("diagnose")
}

# Internal: standard grouped Cox training data.
cox_data <- function(formula, data, group, offset) {
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  time <- y[, 1]
  status <- y[, 2]
  x <- stats::model.matrix(formula, data)[, -1]

  # Properties of the data
  n_samples <- nrow(x)

  # Standardize the covariates
  x <- scale(x)
  x_center <- attr(x, "scaled:center")
  x_scale <- attr(x, "scaled:scale")

  # Check the offset and group arguments
  if (missing(offset) || is.null(offset)) offset <- rep(0.0, n_samples)
  if (missing(group) || is.null(group)) group <- rep(0, n_samples)
  if (!is.factor(group)) group <- factor(group)

  # Sort the data by time
  sorted <- order(time, decreasing = TRUE)
  time <- time[sorted]
  status <- status[sorted]
  x <- x[sorted, , drop = FALSE]
  attr(x, "scaled:center") <- x_center
  attr(x, "scaled:scale") <- x_scale

  offset <- offset[sorted]
  group <- group[sorted]

  list(
    x = x, time = time, status = status, group = group, offset = offset
  )
}

# Internal: capture the rhs design specification used at fit time.
cox_design <- function(formula, data) {
  mf <- stats::model.frame(formula, data)
  terms_rhs <- stats::delete.response(stats::terms(mf))
  x_template <- stats::model.matrix(terms_rhs, mf)
  list(
    terms_rhs = terms_rhs,
    xlevels = stats::.getXlevels(stats::terms(mf), mf),
    x_contrasts = attr(x_template, "contrasts")
  )
}

# Internal: replay the training design-matrix construction on new data.
cox_model_matrix <- function(design_info, newdata, feature_names) {
  x <- stats::model.matrix(design_info$terms_rhs,
    newdata,
    xlev = design_info$xlevels,
    contrasts.arg = design_info$x_contrasts
  )
  if ("(Intercept)" %in% colnames(x)) {
    x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  }
  if (!setequal(colnames(x), feature_names)) {
    stop("newdata design matrix columns do not match the fitted model")
  }
  x <- x[, feature_names, drop = FALSE]
  if (!is.null(design_info$x_center)) {
    x <- sweep(x, 2, design_info$x_center[feature_names], "-")
  }
  x
}

# Internal: validate and align group labels at prediction time.
cox_prediction_group <- function(newgroup, group_levels) {
  group <- factor(newgroup, levels = group_levels)
  if (anyNA(group)) {
    invalid_groups <- unique(as.character(newgroup[is.na(group)]))
    stop(
      "newgroup contains unseen group labels: ",
      paste(invalid_groups, collapse = ", ")
    )
  }
  group
}

# Internal: score observations with group-specific coefficient columns.
cox_linear_predictor <- function(x, group, coefficients) {
  group_levels <- colnames(coefficients)
  lp <- numeric(nrow(x))
  for (k in seq_along(group_levels)) {
    idx <- which(group == group_levels[k])
    if (length(idx) > 0L) {
      lp[idx] <- x[idx, , drop = FALSE] %*% coefficients[, k]
    }
  }
  lp
}

# Internal: compute risk set from hazard, time, and group indices
cox_risk_set <- function(hazard, time, group_idxs) {
  risk_set <- numeric(length(hazard))
  for (k in seq_along(group_idxs)) {
    idx <- group_idxs[[k]]
    risk_set[idx] <- ave_max(cumsum(hazard[idx]), time[idx])
  }
  risk_set
}

# Internal: compute per-group linear predictor offset from theta matrix.
cox_theta_offset <- function(theta, n_features, n_groups, x_by_group,
                             stacked_group_idxs) {
  theta_mat <- matrix(as.numeric(theta), nrow = n_features)
  offset <- numeric(sum(lengths(stacked_group_idxs)))
  for (k in seq_len(n_groups)) {
    offset[stacked_group_idxs[[k]]] <-
      x_by_group[[k]] %*% theta_mat[, k]
  }
  offset
}

# Internal: construct a dense block-diagonal matrix from a list of blocks.
dense_block_diag <- function(blocks) {
  total_rows <- sum(vapply(blocks, nrow, 1L))
  total_cols <- sum(vapply(blocks, ncol, 1L))
  result <- matrix(0, total_rows, total_cols)
  r_off <- 0L
  c_off <- 0L
  for (blk in blocks) {
    nr <- nrow(blk)
    nc <- ncol(blk)
    result[r_off + seq_len(nr), c_off + seq_len(nc)] <- blk
    r_off <- r_off + nr
    c_off <- c_off + nc
  }
  result
}

null_basis <- function(constraints, n_parameters, tol = 1e-8) {
  if (is.null(constraints) || nrow(constraints) == 0L) {
    basis <- diag(n_parameters)
  } else {
    qr_obj <- qr(t(constraints), tol = tol, LAPACK = FALSE)
    rank <- qr_obj$rank
    if (rank >= n_parameters) {
      basis <- matrix(0, nrow = n_parameters, ncol = 0L)
    } else {
      q_full <- qr.Q(qr_obj, complete = TRUE)
      basis <- q_full[, seq.int(rank + 1L, n_parameters), drop = FALSE]
    }
  }
  colnames(basis) <- if (ncol(basis) > 0L) paste0("phi", seq_len(ncol(basis))) else character(0)
  basis
}

safe_solve <- function(x) {
  tryCatch(
    solve(x),
    error = function(e) MASS::ginv(x)
  )
}

#' Compute penalty for a numeric vector using specified regularization type.
#'
#' This function computes the penalty for a given parameter vector `x` based on
#' the specified penalty approach.
#'
#' @param x A numeric vector for which the penalty is computed.
#' @param penalty A character string specifying the penalty type. Valid options
#' are "lasso", "MCP", or "SCAD".
#' @param lambda A numeric value representing the regularization parameter.
#' @param gamma A numeric value used in the penalty for MCP and SCAD.
#'
#' @return A numeric value representing the computed penalty.
#'
#' @details The computation differs according to the penalty type:
#'   \describe{
#'     \item{lasso}{Returns the L1 norm of `x` scaled by `lambda`.}
#'     \item{MCP}{Applies a minimax concave penalty where the penalty function
#'                behaves linearly when the absolute values are below
#'                `lambda * gamma` and quadratically otherwise.}
#'     \item{SCAD}{Applies the smoothly clipped absolute deviation penalty;
#'                 different expressions are used when `abs(x)` is below
#'                 `lambda`, between `lambda` and `gamma * lambda`, or above
#'                 these thresholds.}
#'   }
penalty_value <- function(x, penalty, lambda, gamma) {
  if (lambda == 0) {
    return(0)
  }
  x_abs <- abs(x)
  switch(penalty,
    lasso = lambda * sum(x_abs),
    MCP = {
      condition <- x_abs <= lambda * gamma
      sum(lambda * x_abs[condition] - 0.5 * x_abs[condition]^2 / gamma) +
        sum(0.5 * lambda * gamma^2 * (!condition))
    },
    SCAD = {
      condition1 <- x_abs <= lambda
      condition2 <- x_abs > lambda & x_abs <= gamma * lambda
      sum(lambda * x_abs[condition1]) +
        sum(
          (
            2 * gamma * lambda * x_abs[condition2] -
              x_abs[condition2]^2 -
              lambda^2
          ) / (2 * (gamma - 1))
        ) +
        sum(0.5 * (gamma + 1) * lambda^2 * (!condition1 & !condition2))
    },
    stop("Invalid penalty type. Please choose 'lasso', 'MCP', or 'SCAD'.")
  )
}

coxtrans_null_scores <- function(formula, data, group, offset = NULL) {
  data <- cox_data(formula, data, group, offset)
  x <- data$x
  time <- data$time
  status <- data$status
  group <- data$group
  offset <- data$offset

  group_levels <- levels(group)
  n_features <- ncol(x)
  scores <- matrix(0, nrow = n_features, ncol = length(group_levels))
  colnames(scores) <- group_levels
  rownames(scores) <- colnames(x)

  for (i in seq_along(group_levels)) {
    idx <- which(group == group_levels[i])
    wls <- approx_likelihood(offset[idx], time[idx], status[idx])
    if (length(idx) > 1) {
      scores[, i] <- colMeans(
        sweep(x[idx, , drop = FALSE], 1, wls$residuals * wls$weights, `*`)
      )
    }
  }

  scores
}

coxtrans_penalty_bounds <- function(
  formula, data, group, target, prior_matrix = NULL, offset = NULL
) {
  if (!is.factor(group)) group <- factor(group)
  target_level <- as.character(target)
  scores <- coxtrans_null_scores(formula, data, group, offset)
  if (!target_level %in% colnames(scores)) {
    stop("target '", target_level, "' not found in group levels")
  }
  lambda1 <- max(abs(scores[, target_level]), na.rm = TRUE)
  source_levels <- setdiff(colnames(scores), target_level)
  if (length(source_levels) == 0) {
    return(list(lambda1 = lambda1, lambda2 = 0, lambda3 = 0))
  }
  target_scores <- scores[, target_level]
  local_scores <- sweep(scores[, source_levels, drop = FALSE], 1, target_scores)
  lambda2 <- max(abs(local_scores), na.rm = TRUE)
  if (is.null(prior_matrix)) {
    source_sizes <- tabulate(group)[match(source_levels, levels(group))]
    prior_matrix <- matrix(source_sizes / sum(source_sizes), nrow = 1)
    colnames(prior_matrix) <- source_levels
  }
  if (ncol(prior_matrix) != length(source_levels)) {
    stop(
      "prior_matrix must have ", length(source_levels),
      " columns (one per source)"
    )
  }

  source_scores <- scores[, source_levels, drop = FALSE]
  lambda3 <- vapply(seq_len(nrow(prior_matrix)), function(g) {
    prior_scores <- as.numeric(
      prior_matrix[g, , drop = FALSE] %*% t(source_scores)
    )
    max(abs(target_scores - prior_scores), na.rm = TRUE)
  }, numeric(1))
  names(lambda3) <- rownames(prior_matrix)
  list(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3)
}

#' Calculate the maximum value of the penalty parameter lambda
#'
#' @param formula A formula expression for regression models, in the form
#' \code{response ~ predictors}. The response must be a survival object as
#' returned by the \link[survival]{Surv} function.
#' @param data A data frame containing the variables in the model.
#' @param group A factor specifying the group of each sample.
#' @param offset A numeric vector specifying the offset.
#'
#' @return The maximum value of the penalty parameter lambda, which shrinks all
#' the coefficients to zero.
#' @export
calc_lambda_max <- function(formula, data, group, offset = NULL) {
  if (!is.factor(group)) group <- factor(group)
  target_levels <- levels(group)
  max(vapply(target_levels, function(target) {
    bounds <- coxtrans_penalty_bounds(
      formula, data, group, target, offset = offset
    )
    max(bounds$lambda1, bounds$lambda2, bounds$lambda3)
  }, numeric(1)))
}

#' Simulate survival data for a multi-source Cox model
#'
#' @param beta A vector of length p representing the common coefficients,
#' where p is the number of features.
#' @param eta A matrix of size p x K representing the group-specific
#' coefficients, where K is the number of groups.
#' @param lambda A vector of length K representing the baseline hazard's scale
#' parameters.
#' @param gamma A vector of length K representing the baseline hazard's shape
#' parameters.
#' @param dist A string specifying the distribution of the baseline hazard,
#' either "exponential", "weibull", or "gompertz".
#' @param maxt A positive number specifying the maximum time to simulate.
#' @param n_samples A vector of length K specifying the number of samples per
#' group, or a single number specifying the total number of samples.
#' @param seed An integer specifying the random seed, with a default value of 0.
#' @param sigma An optional p x p covariance matrix for the covariates. If
#' \code{NULL} (default), an identity matrix is used.
#'
#' @return A data frame with columns "id", "group", "X1", "X2", ..., "Xp",
#' "time", and "status".
#'
#' @export
#'
#' @examples
#' beta <- c(1, 1)
#' eta <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
#' lambda <- c(1, 2)
#' gamma <- c(2, 1)
#' dist <- c("gompertz", "weibull")
#' maxt <- 3
#' n_samples <- 100
#' df <- simsurv_tl(beta, eta, lambda, gamma, dist, maxt, n_samples)
#' df
simsurv_tl <- function(
  beta, eta, lambda, gamma, dist, maxt, n_samples, seed = 0, sigma = NULL
) {
  set.seed(seed)
  n_groups <- ncol(eta)
  n_features <- nrow(eta)
  names(beta) <- stringr::str_c("X", seq_len(n_features))
  if (length(n_samples) == 1) n_samples <- rep(n_samples, n_groups)
  mu <- rep(0, n_features)
  if (is.null(sigma)) sigma <- diag(n_features)
  covs <- c()
  times <- c()
  for (k in 1:n_groups) {
    cov <- data.frame(
      id = seq_len(n_samples[k]),
      MASS::mvrnorm(n_samples[k], mu = mu, Sigma = sigma)
    )
    time <- simsurv::simsurv(
      lambdas = lambda[k], gammas = gamma[k], dist = dist[k],
      x = cov, betas = beta + eta[, k], maxt = maxt
    )
    cov$group <- k
    time$group <- k
    covs <- rbind(covs, cov)
    times <- rbind(times, time)
  }
  df <- merge(covs, times, by = c("id", "group"), all.x = TRUE)
  names(df)[names(df) == "eventtime"] <- "time"
  df$id <- seq_len(nrow(df))
  df
}
