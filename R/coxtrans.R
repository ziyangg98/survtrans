#' Transfer Learning Cox Model with Prior Constraints
#'
#' Fits a Cox proportional hazards model that transfers survival information
#' from source domains to a target domain using three-layer penalization:
#' sparse (lambda1), local sharing (lambda2), and prior transfer (lambda3).
#'
#' @param formula A formula with a \code{\link[survival]{Surv}} response.
#' @param data A data frame containing the variables in the model.
#' @param group A factor indicating the group of each observation.
#' @param target The target group level.
#' @param lambda1 Sparse penalty (non-negative scalar). Default 0.
#' @param lambda2 Local penalty for source-target coefficient differences
#'   (non-negative scalar). Default 0.
#' @param lambda3 Prior transfer penalty (non-negative scalar or vector of
#'   length G). When \code{prior_matrix} has G rows, \code{lambda3} should
#'   have length G; a scalar is recycled. Default 0.
#' @param prior_matrix Optional G x (K-1) weight matrix defining transfer
#'   constraints. Each row specifies a weighted combination of source
#'   coefficients that the target should be close to. If \code{NULL},
#'   a single prior using sample-weighted source means is used.
#' @param penalty Penalty type: \code{"lasso"}, \code{"MCP"}, or
#'   \code{"SCAD"}.
#' @param gamma Concavity parameter for MCP/SCAD. Default 3.7 (SCAD) or
#'   3.0 (MCP).
#' @param vartheta Initial augmented Lagrangian parameter. If \code{NULL},
#'   lasso uses 1.0 and non-convex penalties use a curvature-safe value.
#' @param control A \link{survtrans_control} object.
#' @param ... Additional arguments passed to \code{survtrans_control}.
#'
#' @return An object of class \code{coxtrans} with components:
#'   \item{coefficients}{Matrix (p x K) of group-specific coefficients.
#'     Each column is the full beta for that group (target first).}
#'   \item{prior_matrix}{The prior weight matrix used.}
#'   \item{prior_effects}{Matrix of prior constraint residuals.}
#'   \item{active_local}{Logical matrix (p x K-1) of active local
#'     constraints.}
#'   \item{active_prior}{Logical matrix (p x G) of active prior
#'     constraints.}
#'   \item{iter}{Number of ADMM iterations.}
#'   \item{message}{Convergence message.}
#'   \item{history}{Matrix of per-iteration diagnostics.}
#'
#' @export
#'
#' @examples
#' formula <- Surv(time, status) ~ . - group - id
#'
#' # Basic usage with default prior (sample-weighted source mean)
#' fit <- coxtrans(
#'   formula, sim2, sim2$group, 1,
#'   lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
#' )
#' summary(fit)
#'
#' # Custom prior matrix (two tissue-based priors)
#' pm <- rbind(
#'   tissue_A = c(0.5, 0.5, 0, 0),
#'   tissue_B = c(0, 0, 0.5, 0.5)
#' )
#' colnames(pm) <- c("2", "3", "4", "5")
#' fit2 <- coxtrans(
#'   formula, sim2, sim2$group, 1,
#'   lambda1 = 0.075, lambda2 = 0.04, lambda3 = c(0.04, 0.04),
#'   prior_matrix = pm, penalty = "SCAD"
#' )
#' summary(fit2)
coxtrans <- function(
  formula, data, group, target,
  lambda1 = 0.0, lambda2 = 0.0, lambda3 = 0.0,
  prior_matrix = NULL,
  penalty = c("lasso", "MCP", "SCAD"),
  gamma = switch(penalty,
    SCAD = 3.7,
    MCP = 3,
    1
  ), vartheta = NULL,
  control, ...
) {
  penalty <- match.arg(penalty, choices = c("lasso", "MCP", "SCAD"))
  if (missing(control)) control <- survtrans_control(...)
  vartheta_min <- switch(penalty,
    MCP = 1.5 / gamma,
    SCAD = 1.5 / (gamma - 1),
    0
  )
  if (is.null(vartheta)) vartheta <- if (penalty == "lasso") 1 else vartheta_min
  if (!is.numeric(vartheta) || length(vartheta) != 1L || vartheta <= 0) {
    stop("Invalid value for vartheta")
  }
  vartheta <- max(vartheta, vartheta_min)
  set_omp_threads(control$nthreads)
  if (lambda1 < 0 || lambda2 < 0 || any(lambda3 < 0)) {
    stop("Lambda parameters must be non-negative")
  }

  design_info <- cox_design(formula, data)

  # Preprocess data
  data <- cox_data(formula, data, group = group)
  x <- data$x
  x_scale <- attr(x, "scaled:scale")
  x_center <- attr(x, "scaled:center")
  time <- data$time
  status <- data$status
  group <- droplevels(data$group)

  group_levels <- levels(group)
  target_level <- as.character(target)
  if (!target_level %in% group_levels) {
    stop(
      "target '", target_level, "' not found in group levels: ",
      paste(group_levels, collapse = ", ")
    )
  }
  group_levels <- c(target_level, group_levels[group_levels != target_level])

  n_groups <- length(group_levels)
  group_idxs <- lapply(group_levels, function(g) which(group == g))
  n_samples_group <- lengths(group_idxs)
  n_samples_total <- nrow(x)
  n_features <- ncol(x)
  n_parameters <- n_features * n_groups

  # Setup prior_matrix
  if (is.null(prior_matrix)) {
    source_sizes <- n_samples_group[-1]
    prior_matrix <- matrix(source_sizes / sum(source_sizes), nrow = 1)
    colnames(prior_matrix) <- group_levels[-1]
  }
  n_priors <- nrow(prior_matrix)
  if (ncol(prior_matrix) != n_groups - 1) {
    stop("prior_matrix must have ", n_groups - 1, " columns (one per source)")
  }
  if (length(lambda3) == 1L) lambda3 <- rep(lambda3, n_priors)
  if (length(lambda3) != n_priors) {
    stop("lambda3 must be length 1 or ", n_priors, " (one per prior)")
  }

  constraints <- coxtrans_build_constraints(n_features, n_groups, prior_matrix)
  contr_pen <- constraints$contr_pen
  contr_cross <- constraints$contr_cross
  blocks <- coxtrans_build_penalty_blocks(
    constraints$sparse_idx, constraints$local_idx, constraints$prior_idx,
    lambda1, lambda2, lambda3
  )

  # Per-group data (stacked by group order)
  x_by_group <- lapply(group_idxs, function(idx) x[idx, , drop = FALSE])
  time_stacked <- time[unlist(group_idxs)]
  status_stacked <- status[unlist(group_idxs)]

  sizes <- lengths(group_idxs)
  ends <- cumsum(sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  stacked_group_idxs <- Map(seq, starts, ends)

  # Initialize ADMM variables
  theta <- numeric(n_parameters)
  eta <- numeric(nrow(contr_pen))
  nu <- numeric(nrow(contr_pen))

  # Initialize the training process
  n_iterations <- 0
  msg <- ""
  converged <- FALSE
  history <- matrix(NA_real_, nrow = control$maxit, ncol = 9)
  colnames(history) <- c(
    "Iteration", "Primal.Residual", "Dual.Residual", "Primal.Epsilon",
    "Dual.Epsilon", "Augmented.Parameter", "Log.Likelihood",
    "Penalty.Term", "Total.Loss"
  )

  w <- numeric(n_samples_total)
  z <- numeric(n_samples_total)

  # Block index ranges for assembling xwx/xwz
  blk_idx <- lapply(seq_len(n_groups), function(k) {
    ((k - 1) * n_features + 1):(k * n_features)
  })

  repeat {
    n_iterations <- n_iterations + 1

    # Compute per-group offset from theta
    offset <- cox_theta_offset(
      theta, n_features, n_groups, x_by_group, stacked_group_idxs
    )

    # Calculate the weights and working response (IRLS)
    for (k in seq_len(n_groups)) {
      idx <- stacked_group_idxs[[k]]
      wls <- approx_likelihood(
        offset = offset[idx],
        time = time_stacked[idx],
        status = status_stacked[idx]
      )
      w[idx] <- wls$weights
      z[idx] <- wls$residuals + offset[idx]
    }

    # Block-structure assembly for X'WX and X'Wz
    xwx <- matrix(0, n_parameters, n_parameters)
    xwz <- numeric(n_parameters)
    for (k in seq_len(n_groups)) {
      idx <- stacked_group_idxs[[k]]
      w_k <- w[idx]
      z_k <- z[idx]
      a_k <- crossprod(x_by_group[[k]], w_k * x_by_group[[k]]) / n_samples_total
      b_k <- crossprod(x_by_group[[k]], w_k * z_k) / n_samples_total
      xwx[blk_idx[[k]], blk_idx[[k]]] <- a_k
      xwz[blk_idx[[k]]] <- b_k
    }

    # Solve the linear system
    lhs <- xwx + vartheta * contr_cross
    rhs <- xwz + vartheta * crossprod(contr_pen, eta - nu / vartheta)
    theta <- solve(lhs, rhs)

    # Constraint products
    c_theta <- as.numeric(contr_pen %*% theta)

    # Update auxiliary variables
    eta_old <- eta
    eta <- coxtrans_apply_penalty_blocks(
      c_theta + nu / vartheta, blocks, vartheta, penalty, gamma
    )
    nu <- nu + vartheta * (c_theta - eta)

    # Primal and dual residuals
    r_norm <- sqrt(sum((c_theta - eta)^2))
    s_norm <- vartheta * sqrt(sum(crossprod(contr_pen, eta - eta_old)^2))

    eps_pri <- sqrt(nrow(contr_pen)) * control$abstol +
      control$reltol * max(sqrt(sum(c_theta^2)), sqrt(sum(eta^2)))
    dual_vec <- crossprod(contr_pen, nu)
    eps_dual <- sqrt(n_parameters) * control$abstol +
      control$reltol * sqrt(sum(dual_vec^2))
    vartheta_next <- vartheta
    if (r_norm / eps_pri > 10 * s_norm / eps_dual) {
      vartheta_next <- 2 * vartheta
    } else if (s_norm / eps_dual > 10 * r_norm / eps_pri) {
      vartheta_next <- max(vartheta / 2, vartheta_min)
    }

    # Compute loss (recompute offset from updated theta)
    offset <- cox_theta_offset(
      theta, n_features, n_groups, x_by_group, stacked_group_idxs
    )
    hazard <- exp(offset)
    risk_set <- cox_risk_set(hazard, time_stacked, stacked_group_idxs)
    loss <- sum(status_stacked * (offset - log(risk_set)))

    loss_penalty <- coxtrans_penalty_value(
      c_theta, blocks, penalty, gamma
    ) * n_samples_total
    loss_total <- loss - loss_penalty

    # Check convergence
    if (r_norm < eps_pri && s_norm < eps_dual) {
      converged <- TRUE
      msg <- sprintf("Convergence reached at iteration %d.", n_iterations)
    } else if (is.infinite(loss) || is.nan(loss)) {
      converged <- TRUE
      msg <- "Log-likelihood is not finite. Stopping."
    } else if (n_iterations >= control$maxit) {
      converged <- TRUE
      msg <- sprintf(
        "Maximum iterations reached (%d).", control$maxit
      )
    }

    if (control$verbose) {
      cat(sprintf(
        "Iter %d | primal %.4e (tol %.4e) | dual %.4e (tol %.4e) | loss %.4f\n",
        n_iterations, r_norm, eps_pri, s_norm, eps_dual, loss_total
      ))
    }
    history[n_iterations, ] <- c(
      n_iterations, r_norm, s_norm, eps_pri, eps_dual, vartheta,
      loss, loss_penalty, loss_total
    )

    if (converged) break
    vartheta <- vartheta_next
  }

  # Post-processing: recover coefficients from ADMM solution
  coefficients <- matrix(theta, nrow = n_features, ncol = n_groups)

  # Extract active constraint flags from eta (tolerance-based comparison)
  eps <- .Machine$double.eps^0.5
  eta_sparse <- eta[constraints$sparse_idx]
  eta_local <- matrix(
    eta[constraints$local_idx], nrow = n_features, ncol = n_groups - 1
  )
  eta_prior <- vapply(constraints$prior_idx, function(idx) eta[idx], numeric(n_features))

  active_local <- abs(eta_local) < eps # p x (K-1)
  active_prior <- abs(eta_prior) < eps # p x G
  active_sparse <- abs(eta_sparse) < eps # p

  # Post-processing: jointly resolve constraints per feature
  for (j in seq_len(n_features)) {
    merged_srcs <- which(active_local[j, ]) # source indices fused to target

    if (active_sparse[j]) {
      # Sparse: β₁ = 0; fused sources are also zero
      coefficients[j, 1] <- 0
      if (length(merged_srcs) > 0) coefficients[j, merged_srcs + 1] <- 0
      next
    }

    active_g <- which(active_prior[j, ])
    if (length(active_g) > 0) {
      # Prior active: reconstruct target from source weighted average
      g <- active_g[1]
      prior_w <- prior_matrix[g, ]
      is_merged <- rep(FALSE, n_groups - 1)
      is_merged[merged_srcs] <- TRUE
      w_local <- sum(prior_w[is_merged])
      w_free <- sum(prior_w[!is_merged])
      if (w_free > 0) {
        # β₁ = Σ_{free sources} w_i β_{i+1} / (1 - w_local)
        free_srcs <- which(!is_merged)
        coefficients[j, 1] <- sum(
          prior_w[free_srcs] * coefficients[j, free_srcs + 1]
        ) / (1 - w_local)
      } else {
        # All sources fused to target: β₁ = mean(all groups)
        coefficients[j, 1] <- mean(coefficients[j, ])
      }
      # Prior-induced sparsity: reconstruction below ADMM primal tolerance
      # is indistinguishable from zero (symmetric to active_sparse).
      if (abs(coefficients[j, 1]) < control$abstol) {
        coefficients[j, 1] <- 0
      }
      # Sync fused sources to target
      if (length(merged_srcs) > 0) {
        coefficients[j, merged_srcs + 1] <- coefficients[j, 1]
      }
    } else if (length(merged_srcs) > 0) {
      # Local fusion only, no prior: average fused groups
      fused_groups <- c(1, merged_srcs + 1)
      coefficients[j, fused_groups] <- mean(coefficients[j, fused_groups])
    }
  }

  colnames(coefficients) <- group_levels
  rownames(coefficients) <- colnames(x)

  coefficients <- sweep(coefficients, 1, x_scale, "/")
  x <- sweep(x, 2, x_scale, "*")

  # Prior effects (constraint residuals)
  prior_effects <- t(eta_prior)
  rownames(prior_effects) <- rownames(prior_matrix)
  colnames(prior_effects) <- colnames(x)

  # Return the fitted model
  fit <- list(
    coefficients = coefficients,
    prior_matrix = prior_matrix,
    prior_effects = prior_effects,
    active_local = active_local,
    active_prior = active_prior,
    iter = n_iterations,
    message = msg,
    history = history[seq_len(n_iterations), , drop = FALSE],
    penalty = penalty,
    lambda1 = lambda1,
    lambda2 = lambda2,
    lambda3 = lambda3,
    gamma = gamma,
    formula = formula,
    design_info = utils::modifyList(design_info, list(
      x_center = x_center,
      x_scale = x_scale
    )),
    call = match.call(),
    time = time,
    status = status,
    group = group,
    x = x
  )
  class(fit) <- "coxtrans"
  fit
}

coxtrans_build_constraints <- function(n_features, n_groups, prior_matrix) {
  n_priors <- nrow(prior_matrix)
  n_parameters <- n_features * n_groups
  n_constraints <- n_features * (1L + (n_groups - 1L) + n_priors)
  contr_pen <- matrix(0, n_constraints, n_parameters)

  sparse_idx <- seq_len(n_features)
  for (j in seq_len(n_features)) {
    contr_pen[j, j] <- 1
  }

  local_start <- n_features
  local_idx <- local_start + seq_len(n_features * (n_groups - 1L))
  if (n_groups > 1L) {
    for (k in seq.int(2L, n_groups)) {
      for (j in seq_len(n_features)) {
        row <- local_start + (k - 2L) * n_features + j
        contr_pen[row, (k - 1L) * n_features + j] <- 1
        contr_pen[row, j] <- -1
      }
    }
  }

  prior_start <- n_features + n_features * (n_groups - 1L)
  prior_idx <- lapply(seq_len(n_priors), function(g) {
    prior_start + (g - 1L) * n_features + seq_len(n_features)
  })
  for (g in seq_len(n_priors)) {
    for (j in seq_len(n_features)) {
      row <- prior_start + (g - 1L) * n_features + j
      contr_pen[row, j] <- 1
      for (i in seq_len(n_groups - 1L)) {
        contr_pen[row, i * n_features + j] <- -prior_matrix[g, i]
      }
    }
  }

  list(
    contr_pen = contr_pen,
    contr_cross = crossprod(contr_pen),
    sparse_idx = sparse_idx,
    local_idx = local_idx,
    prior_idx = prior_idx
  )
}

coxtrans_build_penalty_blocks <- function(sparse_idx, local_idx, prior_idx,
                                          lambda1, lambda2, lambda3) {
  c(
    list(list(idx = sparse_idx, lambda = lambda1)),
    if (length(local_idx) > 0L) list(list(idx = local_idx, lambda = lambda2)),
    lapply(seq_along(prior_idx), function(g) {
      list(idx = prior_idx[[g]], lambda = lambda3[g])
    })
  )
}

coxtrans_apply_penalty_blocks <- function(
  values, blocks, vartheta, penalty, gamma
) {
  for (block in blocks) {
    values[block$idx] <- threshold_prox(
      values[block$idx], vartheta, penalty, block$lambda, gamma
    )
  }
  values
}

coxtrans_penalty_value <- function(values, blocks, penalty, gamma) {
  sum(vapply(blocks, function(block) {
    penalty_value(values[block$idx], penalty, block$lambda, gamma)
  }, numeric(1)))
}

#' @importFrom utils head
coxtrans_coefficient_map <- function(coefficients, active_prior = NULL,
                                     prior_matrix = NULL) {
  stopifnot(is.matrix(coefficients), ncol(coefficients) >= 1)
  n_groups <- ncol(coefficients)
  n_features <- nrow(coefficients)
  has_prior <- !is.null(active_prior) && !is.null(prior_matrix)

  target_constrained <- if (has_prior) {
    apply(active_prior, 1, any)
  } else {
    rep(FALSE, n_features)
  }

  phi_list <- lapply(seq_len(n_features), function(j) {
    if (target_constrained[j]) {
      unique(coefficients[j, -1])
    } else {
      unique(coefficients[j, ])
    }
  })

  n_unique <- lengths(phi_list)
  offsets <- c(0L, head(cumsum(n_unique), -1L))
  link <- matrix(0, n_groups * n_features, sum(n_unique))

  for (j in seq_len(n_features)) {
    phi_j <- phi_list[[j]]
    off_j <- offsets[j]

    for (k in if (target_constrained[j]) 2:n_groups else seq_len(n_groups)) {
      row <- (k - 1L) * n_features + j
      link[row, off_j + which.min(abs(coefficients[j, k] - phi_j))] <- 1
    }

    if (target_constrained[j]) {
      g_active <- which(active_prior[j, ])[1]
      w <- prior_matrix[g_active, ]
      for (i in seq_len(n_groups - 1)) {
        pos <- which.min(abs(coefficients[j, i + 1] - phi_j))
        link[j, off_j + pos] <- link[j, off_j + pos] + w[i]
      }
    }
  }
  link
}

#' Print method for a \code{coxtrans} object
#'
#' @param x An object of class \code{coxtrans}.
#' @param ... Additional arguments passed to \code{print.summary.coxtrans}.
#' @return Invisibly returns \code{x}.
#' @export
print.coxtrans <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Extract the coefficients from a \code{coxtrans} object
#' @param object An object of class \code{coxtrans}.
#' @param ... Additional arguments (unused).
#' @return A named numeric vector containing the coefficients of the fitted
#' \code{coxtrans} object. The names indicate the group(s) to which the
#' coefficients belong. Zero coefficients are removed.
#' @export
coef.coxtrans <- function(object, ...) {
  coefficients <- object$coefficients
  n_features <- nrow(coefficients)
  active_prior <- object$active_prior

  phi_list <- lapply(seq_len(n_features), function(j) {
    if (!is.null(active_prior) && any(active_prior[j, ])) {
      # Target constrained by prior: free params are source unique values
      unique(coefficients[j, -1])
    } else {
      unique(coefficients[j, ])
    }
  })

  psi <- unlist(phi_list)
  names(psi) <- unlist(lapply(seq_along(phi_list), function(j) {
    vals <- phi_list[[j]]
    stringr::str_c(rownames(coefficients)[j], seq_along(vals), sep = ".")
  }))
  psi
}

#' Log-likelihood for a \code{coxtrans} object
#'
#' @param object An object of class \code{coxtrans}.
#' @param ... Additional arguments (unused).
#' @return A numeric value representing the log-likelihood of the fitted
#' \code{coxtrans} object.
#' @export
logLik.coxtrans <- function(object, ...) {
  group_levels <- colnames(object$coefficients)
  group_idxs <- lapply(group_levels, function(g) which(object$group == g))
  lp <- cox_linear_predictor(object$x, object$group, object$coefficients)
  risk_set <- cox_risk_set(exp(lp), object$time, group_idxs)
  sum(object$status * (lp - log(risk_set)))
}

#' Prediction method for \code{coxtrans} objects.
#' @param object An object of class \code{coxtrans}.
#' @param newdata Optional new data for making predictions. If omitted,
#'   predictions are made using the data used for fitting the model.
#' @param newgroup Optional new group for making predictions. If omitted,
#'   predictions are made using the groups from the original data.
#' @param type The type of prediction to perform. Options include:
#'   \describe{
#'     \item{\code{"lp"}}{The linear predictor.}
#'     \item{\code{"risk"}}{The risk score \eqn{\exp(\text{lp})}.}
#'   }
#' @param ... Additional arguments (unused).
#' @return A numeric vector of predictions.
#' @export
predict.coxtrans <- function(
  object, newdata = NULL, newgroup = NULL,
  type = c("lp", "risk"), ...
) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    x <- object$x
    group <- object$group
  } else {
    x <- cox_model_matrix(
      object$design_info,
      newdata,
      rownames(object$coefficients)
    )
    if (is.null(newgroup)) {
      if ("group" %in% names(newdata)) {
        newgroup <- newdata$group
      } else {
        stop("newgroup must be supplied when newdata has no group column")
      }
    }
    group <- cox_prediction_group(newgroup, levels(object$group))
  }

  lp <- cox_linear_predictor(x, group, object$coefficients)
  if (type == "risk") lp <- exp(lp)
  lp
}

#' Predict the cumulative baseline hazard function for \code{coxtrans} objects
#'
#' @param object An object of class \code{coxtrans}.
#' @param ... Additional arguments (unused).
#'
#' @return A \code{data.frame} with one row for each time point, and columns
#' containing the event time, the cumulative baseline hazard function, and the
#' strata.
#' @export
basehaz.coxtrans <- function(object, ...) {
  group_levels <- colnames(object$coefficients)
  group_idxs <- lapply(group_levels, function(g) which(object$group == g))
  lp <- cox_linear_predictor(object$x, object$group, object$coefficients)
  risk_set <- cox_risk_set(exp(lp), object$time, group_idxs)

  basehaz_list <- lapply(seq_along(group_idxs), function(k) {
    idx <- group_idxs[[k]]
    time_rev <- rev(object$time[idx])
    status_rev <- rev(object$status[idx])
    bh <- cumsum(status_rev / rev(risk_set[idx]))
    data.frame(
      time    = time_rev[status_rev == 1],
      basehaz = bh[status_rev == 1],
      strata  = group_levels[k]
    )
  })
  do.call(rbind, basehaz_list)
}

#' Variance-covariance matrix for a \code{coxtrans} object.
#' @param object An object of class \code{coxtrans}.
#' @param ... Additional arguments (unused).
#' @return A matrix representing the variance-covariance matrix of the
#' coefficients.
#' @export
vcov.coxtrans <- function(object, ...) {
  time <- object$time
  status <- object$status
  group <- object$group
  x <- object$x

  n_samples <- nrow(x)
  coefficients <- object$coefficients
  n_groups <- ncol(coefficients)
  coef_groups <- colnames(coefficients)
  group_idxs <- lapply(coef_groups, function(g) which(group == g))

  psi <- coef(object)
  link_matrix <- coxtrans_coefficient_map(
    coefficients, object$active_prior, object$prior_matrix
  )

  z <- dense_block_diag(lapply(group_idxs, function(idx) x[idx, ])) %*%
    link_matrix
  time <- unlist(lapply(group_idxs, function(idx) time[idx]))
  status <- unlist(lapply(group_idxs, function(idx) status[idx]))

  is_nonzero <- as.vector(psi) != 0
  n_nonzero <- sum(is_nonzero)
  z1 <- as.matrix(z[, is_nonzero, drop = FALSE])
  psi1 <- psi[is_nonzero]
  lp <- z1 %*% psi1

  sizes_v <- lengths(group_idxs)
  ends_v <- cumsum(sizes_v)
  starts_v <- c(1L, head(ends_v, -1L) + 1L)
  stacked_idxs <- Map(seq, starts_v, ends_v)

  gradients <- matrix(0, nrow = n_samples, ncol = n_nonzero)
  hessians <- matrix(0, nrow = n_samples, ncol = n_nonzero^2)
  for (k in seq_len(n_groups)) {
    idx <- stacked_idxs[[k]]
    ghs <- calc_grad_hess(lp[idx], z1[idx, ], time[idx], status[idx])
    gradients[idx, ] <- ghs$grad
    hessians[idx, ] <- ghs$hess
  }
  hess <- matrix(colSums(hessians), n_nonzero, n_nonzero)
  hess_inv <- safe_solve(hess)
  grad_cov <- crossprod(gradients)
  vcov <- hess_inv %*% grad_cov %*% hess_inv
  dimnames(vcov) <- list(
    names(psi1), names(psi1)
  )
  vcov
}

#' Summary method for a \code{coxtrans} object
#'
#' @param object An object of class \code{coxtrans}.
#' @param conf.int A numeric value between 0 and 1 indicating the confidence
#' level of the confidence interval. Default is 0.95.
#' @param target_only Logical; if \code{TRUE}, only the coefficients for the
#' target group are shown in the summary. Default is \code{TRUE}.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class \code{summary.coxtrans}, with the following
#' components:
#' \item{\code{n}, \code{nevent}}{Number of observations and number of events,
#' respectively, in the fit.}
#' \item{\code{logLik}}{The log partial likelihood at the final value.}
#' \item{\code{coefficients}}{A matrix with one row for each coefficient, and
#' columns containing the coefficient, the hazard ratio exp(coef), standard
#' error, Wald statistic, and P value.}
#' \item{\code{conf.int}}{A matrix with one row for each coefficient, containing
#' the confidence limits for exp(coef).}
#'
#' @export
summary.coxtrans <- function(object, conf.int = 0.95, target_only = TRUE, ...) {
  # Extract necessary components from the object
  n_samples <- nrow(object$x)
  n_events <- sum(object$status)
  loglik <- logLik(object)

  # Standard errors
  vcov_matrix <- vcov(object)
  if (is.null(vcov_matrix)) {
    stop("Variance-covariance matrix is not available.")
  }

  coefficients <- coef(object)
  is_nonzero <- coefficients != 0
  coefficients <- coefficients[is_nonzero]
  if (target_only) {
    n_features <- nrow(object$coefficients)
    link_matrix <- as.matrix(coxtrans_coefficient_map(
      object$coefficients, object$active_prior, object$prior_matrix
    )[
      seq_len(n_features), is_nonzero,
      drop = FALSE
    ])
    coefficients <- as.vector(link_matrix %*% coefficients)
    vcov_matrix <- link_matrix %*% vcov_matrix %*% t(link_matrix)
    names(coefficients) <- rownames(object$coefficients)
  }

  se <- sqrt(diag(vcov_matrix))
  z_scores <- coefficients / se
  p_values <- stats::pchisq(z_scores^2, 1, lower.tail = FALSE)
  coef_matrix <- cbind(
    coefficients, exp(coefficients), se, z_scores, p_values
  )
  dimnames(coef_matrix) <- list(
    names(coefficients), c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
  )

  z <- stats::qnorm((1 + conf.int) / 2)
  conf_int_matrix <- cbind(
    exp(coefficients), exp(-coefficients),
    exp(coefficients - z * se), exp(coefficients + z * se)
  )
  dimnames(conf_int_matrix) <- list(
    names(coefficients), c(
      "exp(coef)", "exp(-coef)",
      stringr::str_c("lower .", round(100 * conf.int, 2)),
      stringr::str_c("upper .", round(100 * conf.int, 2), sep = "")
    )
  )

  # Classify feature structure
  coefs_full <- object$coefficients
  n_feat <- nrow(coefs_full)
  feat_names <- rownames(coefs_full)
  a_sparse <- coefs_full[, 1] == 0
  a_local <- object$active_local
  a_prior <- object$active_prior
  prior_names <- rownames(object$prior_matrix)

  feat_type <- character(n_feat)
  feat_prior_label <- character(n_feat) # which priors are active
  for (j in seq_len(n_feat)) {
    if (a_sparse[j] && all(a_local[j, ])) {
      feat_type[j] <- "zero"
    } else if (all(a_local[j, ])) {
      feat_type[j] <- "global"
    } else if (any(a_prior[j, ]) && !any(a_local[j, ])) {
      feat_type[j] <- "prior"
    } else if (any(a_local[j, ])) {
      feat_type[j] <- "partial"
    } else {
      feat_type[j] <- "local"
    }
    # Record which priors are active
    active_g <- which(a_prior[j, ])
    if (length(active_g) > 0) {
      if (!is.null(prior_names)) {
        feat_prior_label[j] <- paste(prior_names[active_g], collapse = ",")
      } else {
        feat_prior_label[j] <- paste(active_g, collapse = ",")
      }
    }
  }
  names(feat_type) <- feat_names
  names(feat_prior_label) <- feat_names

  # Create a summary list
  summary_list <- list(
    n = n_samples,
    nevent = n_events,
    logLik = loglik,
    call = object$call,
    coefficients = coef_matrix,
    conf.int = conf_int_matrix,
    feature_type = feat_type,
    feature_prior_label = feat_prior_label
  )

  class(summary_list) <- "summary.coxtrans"
  summary_list
}

#' Print method for a \code{summary.coxtrans} object
#' @description This function prints a summary of the results of a
#' \code{coxtrans} model in a formatted and user-friendly manner, including the
#' model call, number of samples, number of events, coefficients, and confidence
#' intervals. It also includes
#' significance stars for p-values.
#' @param x A summary object produced from a fitted \code{coxtrans} model. This
#' object contains information such as model coefficients and confidence
#' intervals.
#' @param digits An integer controlling the number of significant digits to
#' print for numeric values.
#' @param signif.stars Logical; if \code{TRUE}, significance stars are printed
#' along with the p-values.
#' @param ... Additional arguments (unused).
#' @return The function prints the summary of the \code{coxtrans} model and
#' returns the object \code{x} invisibly.
#' @details The function provides a formatted output that includes:
#' \itemize{
#'   \item \strong{Call:} The original function call that produced the model.
#'   \item \strong{n and events:} The total number of samples and the number of
#'         events (e.g., deaths).
#'   \item \strong{Coefficients:} The regression coefficients, their standard
#'         errors, z-values, and p-values, formatted in a table. Significance
#'         stars are shown next to p-values if \code{signif.stars} is
#'         \code{TRUE}.
#'   \item \strong{Confidence intervals:} The exponentiated coefficients along
#'         with their confidence intervals.
#' }
#' @export
print.summary.coxtrans <- function(
  x, digits = max(getOption("digits") - 3, 3),
  signif.stars = getOption("show.signif.stars"), ...
) {
  # Print call
  cat("Call:\n")
  print(x$call)
  cat("\n")

  # Print number of samples and events
  cat("  n=", x$n, ", number of events=", x$nevent, "\n\n", sep = "")

  is_nonzero <- x$coefficients[, "coef"] != 0

  # Print coefficients with formatted output
  stats::printCoefmat(x$coefficients[is_nonzero, , drop = FALSE],
    digits = digits, signif.stars = signif.stars,
    cs.ind = 1:3, tst.ind = 4, P.values = TRUE, has.Pvalue = TRUE
  )

  # Print confidence intervals
  print(
    format(x$conf.int[is_nonzero, , drop = FALSE], digits = digits),
    quote = FALSE
  )

  # Print feature structure
  if (!is.null(x$feature_type)) {
    ft <- x$feature_type
    fpl <- x$feature_prior_label
    cat("\nFeature structure:\n")

    prior_f <- names(ft[ft == "prior"])
    global_f <- names(ft[ft == "global"])
    partial_f <- names(ft[ft == "partial"])
    local_f <- names(ft[ft == "local"])
    zero_f <- names(ft[ft == "zero"])

    # Build label-value pairs, then align output
    entries <- list()
    if (length(prior_f) > 0) {
      prior_groups <- split(prior_f, fpl[prior_f])
      entries <- c(entries, lapply(names(prior_groups), function(pg) {
        label <- if (nchar(pg) > 0 && pg != "1") {
          paste0("Prior [", pg, "]")
        } else {
          "Prior transfer"
        }
        list(label = label, value = paste(prior_groups[[pg]], collapse = ", "))
      }))
    }
    if (length(global_f) > 0) {
      entries <- c(entries, list(list(
        label = "Shared (local)",
        value = paste(global_f, collapse = ", ")
      )))
    }
    if (length(partial_f) > 0) {
      entries <- c(entries, list(list(
        label = "Partial shared",
        value = paste(partial_f, collapse = ", ")
      )))
    }
    if (length(local_f) > 0) {
      entries <- c(entries, list(list(
        label = "Group-specific",
        value = paste(local_f, collapse = ", ")
      )))
    }
    if (length(zero_f) > 0) {
      preview <- paste(
        zero_f[seq_len(min(3L, length(zero_f)))],
        collapse = ", "
      )
      if (length(zero_f) > 3L) preview <- paste0(preview, ", ...")
      entries <- c(entries, list(list(
        label = "Sparse (zero)",
        value = paste0(length(zero_f), " features (", preview, ")")
      )))
    }

    if (length(entries) > 0) {
      max_width <- max(nchar(sapply(entries, `[[`, "label")))
      for (e in entries) {
        cat("  ", formatC(e$label, width = -max_width), ": ",
          e$value, "\n",
          sep = ""
        )
      }
    }
  }

  invisible(x)
}

#' Refit a coxtrans model with hard constraints
#'
#' Extracts the constraint structure (active set, sharing pattern, prior
#' constraints) from a penalized coxtrans fit, encodes them as hard constraints
#' via a link matrix, and solves an unpenalized stratified Cox model in the
#' reparametrized space. Returns debiased MLE coefficients and valid standard
#' errors for target-group features.
#'
#' @param object A fitted \code{coxtrans} object.
#' @param ... Additional arguments (unused).
#'
#' @return An object of class \code{refit.coxtrans} with components:
#'   \item{coefficients}{Matrix (p x 5) with columns \code{coef},
#'     \code{exp(coef)}, \code{se(coef)}, \code{z}, \code{Pr(>|z|)}.
#'     Rows correspond to original features; inactive features have all zeros
#'     except \code{Pr(>|z|)} = 1.}
#'   \item{coxph_fit}{The underlying \code{coxph} object, or \code{NULL}.}
#'   \item{n}{Total number of observations used in refit.}
#'   \item{nevent}{Number of events.}
#'   \item{active_set}{Integer vector of active feature indices.}
#'
#' @details
#' The refit procedure:
#' \enumerate{
#'   \item Build the link matrix \eqn{L} from the penalized fit, encoding
#'     shared coefficients, prior constraints, and group-specific parameters.
#'   \item Construct reparametrized design \eqn{Z = \mathrm{bdiag}(X_1, \ldots,
#'     X_K) L_{[:, \text{nonzero}]}}.
#'   \item Fit \code{coxph(Surv ~ Z + strata(group))} — unpenalized
#'     stratified Cox with hard constraints built into the design.
#'   \item Map estimates back to target coefficients via \eqn{L}, compute
#'     standard errors via delta method.
#' }
#'
#' @export
refit <- function(object, ...) {
  UseMethod("refit")
}

#' @rdname refit
#' @export
refit.coxtrans <- function(object, ...) {
  p <- nrow(object$coefficients)
  feature_names <- rownames(object$coefficients)
  coef_target <- object$coefficients[, 1]
  active_set <- which(abs(coef_target) > 1e-6)

  # Empty result template
  empty_coefmat <- matrix(0, nrow = p, ncol = 5)
  empty_coefmat[, 5] <- 1
  colnames(empty_coefmat) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
  rownames(empty_coefmat) <- feature_names

  # Feature constraint type from penalized fit
  feat_type <- character(p)
  names(feat_type) <- feature_names
  a_local <- object$active_local
  a_prior <- object$active_prior
  for (j in seq_len(p)) {
    if (abs(coef_target[j]) < 1e-6) {
      feat_type[j] <- "zero"
    } else if (any(a_prior[j, ])) {
      feat_type[j] <- "prior"
    } else if (all(a_local[j, ])) {
      feat_type[j] <- "shared"
    } else if (any(a_local[j, ])) {
      feat_type[j] <- "partial"
    } else {
      feat_type[j] <- "specific"
    }
  }

  make_result <- function(coefmat, coxph_fit = NULL, n = 0L, nevent = 0L) {
    out <- list(
      coefficients = coefmat, coxph_fit = coxph_fit,
      n = n, nevent = nevent, active_set = active_set,
      feature_type = feat_type
    )
    class(out) <- "refit.coxtrans"
    out
  }

  if (length(active_set) == 0) {
    return(make_result(empty_coefmat))
  }

  # 1. Link matrix encoding constraint structure
  link <- coxtrans_coefficient_map(
    object$coefficients, object$active_prior, object$prior_matrix
  )
  psi <- coef(object)
  is_nonzero <- psi != 0
  n_free <- sum(is_nonzero)
  if (n_free == 0) {
    return(make_result(empty_coefmat))
  }

  # 2. Reparametrized design matrix Z
  group_levels <- colnames(object$coefficients)
  group_idxs <- lapply(group_levels, function(g) which(object$group == g))

  z_mat <- (dense_block_diag(
    lapply(group_idxs, function(idx) object$x[idx, , drop = FALSE])
  ) %*% link)[, is_nonzero, drop = FALSE]
  colnames(z_mat) <- paste0("psi", seq_len(n_free))

  # 3. Stack response and strata
  time_s <- unlist(lapply(group_idxs, function(idx) object$time[idx]))
  status_s <- unlist(lapply(group_idxs, function(idx) object$status[idx]))
  strata_s <- factor(rep(seq_along(group_idxs), lengths(group_idxs)))

  # 4. Unpenalized stratified Cox
  df_refit <- data.frame(
    z_mat,
    time = time_s,
    status = status_s,
    strata_var = strata_s
  )
  fml <- stats::as.formula(
    paste(
      "survival::Surv(time, status) ~",
      paste(colnames(z_mat), collapse = " + "),
      "+ strata(strata_var)"
    ),
    env = list2env(list(strata = survival::strata), parent = baseenv())
  )
  coxph_fit <- tryCatch(
    survival::coxph(fml, data = df_refit),
    error = function(e) NULL
  )
  if (is.null(coxph_fit)) {
    return(make_result(empty_coefmat))
  }

  # 5. Map back to target coefficients via delta method
  psi_hat <- rep(0, length(psi))
  psi_hat[is_nonzero] <- stats::coef(coxph_fit)

  target_link <- link[seq_len(p), , drop = FALSE]
  beta <- as.vector(target_link %*% psi_hat)

  l_active <- target_link[, is_nonzero, drop = FALSE]
  v_mat <- l_active %*% stats::vcov(coxph_fit) %*% t(l_active)
  se <- sqrt(pmax(diag(v_mat), 0))

  # 6. Build coefficient matrix
  coefmat <- empty_coefmat
  for (j in active_set) {
    z_val <- if (se[j] > 0) beta[j] / se[j] else 0
    p_val <- if (se[j] > 0) stats::pchisq(z_val^2, 1, lower.tail = FALSE) else 1
    coefmat[j, ] <- c(beta[j], exp(beta[j]), se[j], z_val, p_val)
  }

  make_result(coefmat, coxph_fit, n = nrow(df_refit), nevent = sum(status_s))
}

#' @export
print.refit.coxtrans <- function(x, digits = max(getOption("digits") - 3, 3),
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {
  cat("Refit coxtrans (unpenalized MLE, hard constraints)\n")
  cat("  n=", x$n, ", events=", x$nevent, "\n\n", sep = "")

  active <- x$active_set
  if (length(active) == 0) {
    cat("  No active features.\n")
    return(invisible(x))
  }

  # Only print features with p < 0.2 (notable effects)
  pvals <- x$coefficients[active, "Pr(>|z|)"]
  show <- active[pvals < 0.2]
  n_omit <- length(active) - length(show)

  if (length(show) > 0) {
    stats::printCoefmat(
      x$coefficients[show, , drop = FALSE],
      digits = digits, signif.stars = signif.stars,
      cs.ind = 1:3, tst.ind = 4, P.values = TRUE, has.Pvalue = TRUE
    )
  }
  if (n_omit > 0) {
    cat("  (", n_omit, " non-significant features omitted)\n", sep = "")
  }

  # Constraint structure
  ft <- x$feature_type[active]
  types <- list(
    "Prior transfer" = names(ft[ft == "prior"]),
    "Shared (local)" = names(ft[ft == "shared"]),
    "Partial shared" = names(ft[ft == "partial"]),
    "Group-specific" = names(ft[ft == "specific"])
  )
  types <- types[lengths(types) > 0]
  if (length(types) > 0) {
    cat("\nConstraint structure:\n")
    w <- max(nchar(names(types)))
    for (nm in names(types)) {
      cat("  ", formatC(nm, width = -w), ": ",
        paste(types[[nm]], collapse = ", "), "\n",
        sep = ""
      )
    }
  }

  invisible(x)
}

#' Diagnose Cox Transfer Model's Optimization Process
#'
#' @param object An object of class \code{coxtrans}.
#' @param ... Additional arguments (currently unused).
#' @details This function produces two plots:
#' - Residuals Convergence: Plots the evolution of primal and dual residuals
#' along with their tolerance levels.
#' - Loss Decomposition: Plots the negative log-likelihood, total loss, and
#' penalty term.
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#' producing diagnostic plots.
#' @export
diagnose.coxtrans <- function(object, ...) {
  history <- as.list(as.data.frame(object$history))
  iter_range <- range(history$Iteration)
  colors <- c(
    Augmented = "#1B9E77", Primal = "#D95F02", Dual = "#7570B3",
    NLL = "#E7298A", Penalty = "#66A61E", Total = "#E6AB02"
  )
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(
    mfrow = c(1, 2), mar = c(5, 4.5, 4, 4) + 0.1, mgp = c(2.5, 0.8, 0),
    cex.axis = 0.9, cex.lab = 1.1, cex.main = 1.2
  )

  plot_residuals <- function() {
    y_lim <- range(
      0, history$Primal.Residual, history$Dual.Residual,
      history$Primal.Epsilon, history$Dual.Epsilon
    ) * 1.05
    plot(NA,
      xlim = iter_range, ylim = y_lim, xlab = "Iteration", ylab = "Residuals",
      main = "Residuals Convergence"
    )
    grid(col = "gray90", lty = 2)
    lines(
      history$Iteration, history$Primal.Residual,
      col = colors["Primal"], lwd = 2
    )
    lines(
      history$Iteration, history$Dual.Residual,
      col = colors["Dual"], lwd = 2
    )
    lines(
      history$Iteration, history$Primal.Epsilon,
      col = colors["Primal"], lty = 2, lwd = 1.5
    )
    lines(
      history$Iteration, history$Dual.Epsilon,
      col = colors["Dual"], lty = 2, lwd = 1.5
    )
    par(new = TRUE)
    plot(
      history$Iteration, history$Augmented.Parameter,
      type = "l", col = colors["Augmented"], lwd = 2, axes = FALSE,
      xlab = "", ylab = "",
      ylim = range(history$Augmented.Parameter) + c(-0.05, 0.05)
    )
    axis(4, col.axis = colors["Augmented"], col = colors["Augmented"])
    mtext("Augmented Parameter",
      side = 4, line = 2.5, col = colors["Augmented"], cex = 0.9
    )
    legend(
      "top",
      legend = c("Primal", "Dual", "Residual", "Tolerance", expression(rho)),
      col = c(
        colors["Primal"], colors["Dual"], "black", "black",
        colors["Augmented"]
      ),
      lty = c(1, 1, 1, 2, 1),
      lwd = c(2, 2, 1.5, 1.5, 2),
      horiz = FALSE,
      ncol = 3,
      xpd = NA,
      inset = c(0, 0.01),
      cex = 0.75
    )
  }

  plot_objective <- function() {
    y_lim_obj <- range(
      -history$Log.Likelihood, -history$Total.Loss
    ) + c(-0.05, 0.05)
    plot(NA,
      xlim = iter_range, ylim = y_lim_obj, xlab = "Iteration", ylab = "Loss",
      main = "Loss Decomposition"
    )
    grid(col = "gray90", lty = 2)
    lines(
      history$Iteration, -history$Log.Likelihood,
      col = colors["NLL"], lwd = 2
    )
    lines(
      history$Iteration, -history$Total.Loss,
      col = colors["Total"], lwd = 2
    )
    par(new = TRUE)
    plot(
      history$Iteration, history$Penalty.Term,
      type = "l", col = colors["Penalty"], lwd = 2, axes = FALSE,
      xlab = "", ylab = "", ylim = range(history$Penalty.Term) + c(-0.05, 0.05)
    )
    axis(4, col.axis = colors["Penalty"], col = colors["Penalty"])
    mtext("Penalty Term",
      side = 4, line = 2.5, col = colors["Penalty"], cex = 0.9
    )
    legend(
      "top",
      legend = c("Total", "NLL", "Penalty"),
      col = colors[c("Total", "NLL", "Penalty")],
      lty = 1,
      lwd = 2,
      horiz = FALSE,
      ncol = 3,
      xpd = NA,
      inset = c(0, 0.01),
      cex = 0.75
    )
  }

  plot_residuals()
  plot_objective()
}
