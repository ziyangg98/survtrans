#' Symmetric Multi-Task Cox Model with Global Centers
#'
#' Fits a target-free multi-task Cox model across all groups using three
#' symmetric penalties: group-wise sparsity, pairwise fusion, and shrinkage
#' toward user-defined global centers.
#'
#' @param formula A formula with a \code{\link[survival]{Surv}} response.
#' @param data A data frame containing the variables in the model.
#' @param group A factor indicating the group of each observation.
#' @param w A non-negative matrix with one row per global center and one column
#'   per group. Each row must sum to one and defines a convex combination of the
#'   group-specific coefficients.
#' @param lambda1 Sparse penalty applied to every group-specific coefficient.
#' @param lambda2 Fusion penalty applied to every pairwise group difference.
#' @param lambda3 Center penalty (scalar or vector of length \code{nrow(w)}).
#' @param penalty Penalty type: \code{"lasso"}, \code{"MCP"}, or
#'   \code{"SCAD"}.
#' @param gamma Concavity parameter for MCP/SCAD. Default 3.7 (SCAD) or
#'   3.0 (MCP).
#' @param vartheta Fixed augmented Lagrangian parameter. Default 1.0.
#' @param control A \link{survtrans_control} object.
#' @param ... Additional arguments passed to \code{survtrans_control}.
#'
#' @return An object of class \code{coxmtl}.
#' @importFrom stats BIC
#' @export
coxmtl <- function(
  formula, data, group, w,
  lambda1 = 0.0, lambda2 = 0.0, lambda3 = 0.0,
  penalty = c("lasso", "MCP", "SCAD"),
  gamma = switch(penalty,
    SCAD = 3.7,
    MCP = 3,
    1
  ), vartheta = 1.0,
  control, ...
) {
  penalty <- match.arg(penalty, choices = c("lasso", "MCP", "SCAD"))
  if (missing(control)) control <- survtrans_control(...)
  set_omp_threads(control$nthreads)
  if (lambda1 < 0 || lambda2 < 0 || any(lambda3 < 0)) {
    stop("Lambda parameters must be non-negative")
  }

  design_info <- cox_design(formula, data)

  data <- cox_data(formula, data, group = group)
  x <- data$x
  x_center <- attr(x, "scaled:center")
  x_scale <- attr(x, "scaled:scale")
  time <- data$time
  status <- data$status
  group <- droplevels(data$group)

  group_levels <- levels(group)
  n_groups <- length(group_levels)
  n_features <- ncol(x)
  n_parameters <- n_features * n_groups
  feature_names <- colnames(x)
  param_names <- as.vector(outer(feature_names, group_levels, paste, sep = ":"))
  w <- coxmtl_validate_weights(w, group_levels)

  n_centers <- nrow(w)
  if (length(lambda3) == 1L) lambda3 <- rep(lambda3, n_centers)
  if (length(lambda3) != n_centers) {
    stop("lambda3 must be length 1 or ", n_centers, " (one per row of w)")
  }

  constraints <- coxmtl_build_constraints(feature_names, group_levels, w)
  contr_pen <- constraints$contr_pen
  contr_cross <- constraints$contr_cross
  group_idxs <- lapply(group_levels, function(g) which(group == g))
  n_samples_total <- nrow(x)
  x_by_group <- lapply(group_idxs, function(idx) x[idx, , drop = FALSE])
  time_stacked <- unlist(lapply(group_idxs, function(idx) time[idx]))
  status_stacked <- unlist(lapply(group_idxs, function(idx) status[idx]))

  sizes <- lengths(group_idxs)
  ends <- cumsum(sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  stacked_group_idxs <- Map(seq, starts, ends)
  blk_idx <- lapply(seq_len(n_groups), function(k) {
    ((k - 1L) * n_features + 1L):(k * n_features)
  })
  blocks <- coxmtl_build_penalty_blocks(
    constraints$sparse_idx, constraints$pair_rows, constraints$center_idx,
    lambda1, lambda2, lambda3
  )

  theta <- numeric(n_parameters)
  eta <- numeric(nrow(contr_pen))
  nu <- numeric(nrow(contr_pen))

  n_iterations <- 0L
  msg <- ""
  converged <- FALSE
  history <- matrix(NA_real_, nrow = control$maxit, ncol = 9)
  colnames(history) <- c(
    "Iteration", "Primal.Residual", "Dual.Residual", "Primal.Epsilon",
    "Dual.Epsilon", "Augmented.Parameter", "Log.Likelihood",
    "Penalty.Term", "Total.Loss"
  )

  weights <- numeric(n_samples_total)
  z <- numeric(n_samples_total)
  lag_aug_prev <- Inf
  loss_total_prev <- Inf

  repeat {
    n_iterations <- n_iterations + 1L

    offset <- cox_theta_offset(
      theta, n_features, n_groups, x_by_group, stacked_group_idxs
    )
    for (k in seq_len(n_groups)) {
      idx <- stacked_group_idxs[[k]]
      wls <- approx_likelihood(
        offset = offset[idx],
        time = time_stacked[idx],
        status = status_stacked[idx]
      )
      weights[idx] <- wls$weights
      z[idx] <- wls$residuals + offset[idx]
    }

    xwx <- matrix(0, n_parameters, n_parameters)
    xwz <- numeric(n_parameters)
    for (k in seq_len(n_groups)) {
      idx <- stacked_group_idxs[[k]]
      w_k <- weights[idx]
      z_k <- z[idx]
      a_k <- crossprod(x_by_group[[k]], w_k * x_by_group[[k]]) / n_samples_total
      b_k <- crossprod(x_by_group[[k]], w_k * z_k) / n_samples_total
      xwx[blk_idx[[k]], blk_idx[[k]]] <- a_k
      xwz[blk_idx[[k]]] <- b_k
    }

    lhs <- xwx + vartheta * contr_cross
    rhs <- xwz + vartheta * crossprod(contr_pen, eta - nu / vartheta)
    theta <- solve(lhs, rhs)

    c_theta <- as.numeric(contr_pen %*% theta)
    eta_old <- eta
    eta <- coxmtl_apply_penalty_blocks(
      c_theta + nu / vartheta, blocks, vartheta, penalty, gamma
    )
    nu <- nu + vartheta * (c_theta - eta)

    r_norm <- sqrt(sum((c_theta - eta)^2))
    s_norm <- vartheta * sqrt(sum(crossprod(contr_pen, eta - eta_old)^2))

    eps_pri <- sqrt(nrow(contr_pen)) * control$abstol +
      control$reltol * max(sqrt(sum(c_theta^2)), sqrt(sum(eta^2)))
    dual_vec <- crossprod(contr_pen, nu)
    eps_dual <- sqrt(n_parameters) * control$abstol +
      control$reltol * sqrt(sum(dual_vec^2))

    # re-evaluate at updated theta for convergence diagnostics
    offset <- cox_theta_offset(
      theta, n_features, n_groups, x_by_group, stacked_group_idxs
    )
    hazard <- exp(offset)
    risk_set <- cox_risk_set(hazard, time_stacked, stacked_group_idxs)
    loss <- sum(status_stacked * (offset - log(risk_set)))

    loss_penalty <- coxmtl_penalty_value(
      c_theta, blocks, penalty, gamma
    ) * n_samples_total
    loss_total <- loss - loss_penalty
    lag_aug <- loss_total +
      sum(nu * (c_theta - eta)) +
      0.5 * vartheta * sum((c_theta - eta)^2)

    if (r_norm < eps_pri && s_norm < eps_dual) {
      converged <- TRUE
      msg <- sprintf("Convergence reached at iteration %d.", n_iterations)
    } else if (is.infinite(loss) || is.nan(loss)) {
      converged <- TRUE
      msg <- "Log-likelihood is not finite. Stopping."
    } else if (n_iterations > 1L && (
      abs(lag_aug - lag_aug_prev) / (abs(lag_aug_prev) + 1) < control$fdev ||
        abs(loss_total - loss_total_prev) /
          (abs(loss_total_prev) + 1) < control$fdev
    )) {
      converged <- TRUE
      msg <- sprintf("Objective stabilized at iteration %d.", n_iterations)
    } else if (n_iterations >= control$maxit) {
      converged <- TRUE
      msg <- sprintf("Maximum iterations reached (%d).", control$maxit)
    }
    lag_aug_prev <- lag_aug
    loss_total_prev <- loss_total

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
  }

  eps <- .Machine$double.eps^0.5
  binding <- abs(eta) < eps # TRUE where constraint is tight (eta ≈ 0)
  active_sparse <- matrix(
    binding[constraints$sparse_idx],
    nrow = n_features,
    ncol = n_groups
  )
  colnames(active_sparse) <- group_levels
  rownames(active_sparse) <- feature_names

  if (constraints$n_pairs > 0L) {
    active_pair <- matrix(
      binding[constraints$pair_rows],
      nrow = n_features,
      ncol = constraints$n_pairs
    )
    colnames(active_pair) <- constraints$pair_labels
    rownames(active_pair) <- feature_names
  } else {
    active_pair <- matrix(FALSE, nrow = n_features, ncol = 0L)
    rownames(active_pair) <- feature_names
  }

  active_w <- matrix(
    binding[constraints$center_rows],
    nrow = n_features,
    ncol = n_centers * n_groups
  )
  colnames(active_w) <- constraints$center_labels
  rownames(active_w) <- feature_names

  # Project unstandardized theta onto the null space of binding constraints
  binding_constr <- contr_pen[binding, , drop = FALSE]
  basis <- null_basis(binding_constr, n_parameters)
  rownames(basis) <- param_names
  theta_unscaled <- theta / rep(x_scale, times = n_groups)
  phi <- if (ncol(basis) > 0L) as.numeric(crossprod(basis, theta_unscaled)) else numeric(0)
  names(phi) <- colnames(basis)
  beta_exact <- if (ncol(basis) > 0L) as.numeric(basis %*% phi) else numeric(n_parameters)
  coefficients <- matrix(
    beta_exact, nrow = n_features,
    dimnames = list(feature_names, group_levels)
  )
  x <- sweep(x, 2, x_scale, "*")

  fit <- list(
    coefficients = coefficients,
    phi = phi,
    basis = basis,
    w = w,
    active_sparse = active_sparse,
    active_pair = active_pair,
    active_w = active_w,
    pair_index = constraints$pair_index,
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
    x_scale = x_scale,
    call = match.call(),
    time = time,
    status = status,
    group = group,
    x = x
  )
  class(fit) <- "coxmtl"
  fit
}

coxmtl_validate_weights <- function(w, group_levels, tol = 1e-8) {
  if (missing(w) || is.null(w)) {
    stop("w must be provided")
  }
  w <- as.matrix(w)
  if (!is.numeric(w)) {
    stop("w must be a numeric matrix")
  }
  if (nrow(w) < 1L) {
    stop("w must have at least one row")
  }
  if (ncol(w) != length(group_levels)) {
    stop("w must have ", length(group_levels), " columns (one per group)")
  }
  if (!is.null(colnames(w))) {
    if (!setequal(colnames(w), group_levels)) {
      stop("colnames(w) must match group levels")
    }
    w <- w[, group_levels, drop = FALSE]
  } else {
    colnames(w) <- group_levels
  }
  if (is.null(rownames(w))) {
    rownames(w) <- paste0("w", seq_len(nrow(w)))
  }
  if (any(!is.finite(w))) {
    stop("w must contain only finite values")
  }
  if (any(w < 0)) {
    stop("w must be non-negative")
  }
  rs <- rowSums(w)
  if (any(abs(rs - 1) > tol)) {
    stop("Each row of w must sum to 1")
  }
  w
}

coxmtl_build_constraints <- function(feature_names, group_levels, w) {
  n_features <- length(feature_names)
  n_groups <- length(group_levels)
  param_names <- as.vector(outer(feature_names, group_levels, paste, sep = ":"))

  pair_index <- if (n_groups < 2L) {
    matrix(character(0), nrow = 0L, ncol = 2L)
  } else {
    t(utils::combn(group_levels, 2))
  }
  n_pairs <- nrow(pair_index)
  n_centers <- nrow(w)
  pair_labels <- if (n_pairs) {
    paste(pair_index[, 1], pair_index[, 2], sep = "~")
  } else {
    character(0)
  }
  center_labels <- unlist(lapply(rownames(w), function(lbl) {
    paste(lbl, group_levels, sep = "->")
  }))

  sparse_contr <- diag(n_groups)
  pair_contr <- matrix(0, nrow = n_pairs, ncol = n_groups)
  if (n_pairs) {
    pair_contr[
      cbind(seq_len(n_pairs), match(pair_index[, 1], group_levels))
    ] <- 1
    pair_contr[
      cbind(seq_len(n_pairs), match(pair_index[, 2], group_levels))
    ] <- -1
  }
  center_contr <- do.call(rbind, lapply(seq_len(n_centers), function(g) {
    diag(n_groups) - matrix(
      w[g, ],
      nrow = n_groups,
      ncol = n_groups,
      byrow = TRUE
    )
  }))

  group_contr <- rbind(sparse_contr, pair_contr, center_contr)
  contr_pen <- kronecker(group_contr, diag(n_features))
  colnames(contr_pen) <- param_names

  n_sparse <- n_groups * n_features
  n_pair_rows <- n_pairs * n_features
  n_center_rows <- n_centers * n_groups * n_features
  sparse_idx <- seq_len(n_sparse)
  offset_pair <- n_sparse
  offset_center <- n_sparse + n_pair_rows
  center_idx <- lapply(seq_len(n_centers), function(g) {
    offset_center + (g - 1L) * n_groups * n_features +
      seq_len(n_groups * n_features)
  })

  rownames(contr_pen) <- c(
    paste0("sparse:", param_names),
    if (n_pairs) {
      paste0(
        "pair:", rep(pair_labels, each = n_features),
        ":", feature_names
      )
    },
    paste0(
      "center:", rep(rownames(w), each = n_groups * n_features),
      "->", rep(group_levels, times = n_centers, each = n_features),
      ":", feature_names
    )
  )

  list(
    contr_pen = contr_pen,
    contr_cross = crossprod(contr_pen),
    sparse_idx = sparse_idx,
    pair_rows = if (n_pair_rows) {
      offset_pair + seq_len(n_pair_rows)
    } else {
      integer(0)
    },
    center_rows = offset_center + seq_len(n_center_rows),
    center_idx = center_idx,
    pair_index = pair_index,
    pair_labels = pair_labels,
    center_labels = center_labels,
    n_pairs = n_pairs
  )
}

coxmtl_build_penalty_blocks <- function(sparse_idx, pair_rows, center_idx,
                                        lambda1, lambda2, lambda3) {
  c(
    list(list(idx = sparse_idx, lambda = lambda1)),
    if (length(pair_rows) > 0L) list(list(idx = pair_rows, lambda = lambda2)),
    lapply(seq_along(center_idx), function(g) {
      list(idx = center_idx[[g]], lambda = lambda3[g])
    })
  )
}

coxmtl_apply_penalty_blocks <- function(
  values, blocks, vartheta, penalty, gamma
) {
  for (block in blocks) {
    values[block$idx] <- threshold_prox(
      values[block$idx], vartheta, penalty, block$lambda, gamma
    )
  }
  values
}

coxmtl_penalty_value <- function(values, blocks, penalty, gamma) {
  sum(vapply(blocks, function(block) {
    penalty_value(values[block$idx], penalty, block$lambda, gamma)
  }, numeric(1)))
}

#' Print method for a \code{coxmtl} object
#'
#' @param x An object of class \code{coxmtl}.
#' @param ... Additional arguments passed to \code{print.summary.coxmtl}.
#' @return Invisibly returns \code{x}.
#' @export
print.coxmtl <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' Extract the coefficients from a \code{coxmtl} object
#'
#' @param object An object of class \code{coxmtl}.
#' @param ... Additional arguments (unused).
#' @return A named numeric vector containing the free parameters of the fitted
#'   \code{coxmtl} object.
#' @export
coef.coxmtl <- function(object, ...) {
  object$phi
}

#' Log-likelihood for a \code{coxmtl} object
#'
#' @param object An object of class \code{coxmtl}.
#' @param ... Additional arguments (unused).
#' @return A numeric value representing the log-likelihood of the fitted
#'   \code{coxmtl} object.
#' @export
logLik.coxmtl <- function(object, ...) {
  group_levels <- colnames(object$coefficients)
  group_idxs <- lapply(group_levels, function(g) which(object$group == g))
  lp <- cox_linear_predictor(object$x, object$group, object$coefficients)
  risk_set <- cox_risk_set(exp(lp), object$time, group_idxs)
  sum(object$status * (lp - log(risk_set)))
}

#' Prediction method for \code{coxmtl} objects
#'
#' @param object An object of class \code{coxmtl}.
#' @param newdata Optional new data for making predictions. If omitted,
#'   predictions are made using the data used for fitting the model.
#' @param newgroup Optional new group labels for making predictions. If omitted,
#'   predictions use the groups from the original data, or the \code{group}
#'   column in \code{newdata} when available.
#' @param type The type of prediction to perform: \code{"lp"} for the linear
#'   predictor or \code{"risk"} for \eqn{\exp(\text{lp})}.
#' @param ... Additional arguments (unused).
#' @return A numeric vector of predictions.
#' @export
predict.coxmtl <- function(
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

#' Predict the cumulative baseline hazard function for \code{coxmtl} objects
#'
#' @param object An object of class \code{coxmtl}.
#' @param ... Additional arguments (unused).
#' @return A \code{data.frame} with event time, cumulative baseline hazard, and
#'   strata columns.
#' @export
basehaz.coxmtl <- function(object, ...) {
  group_levels <- colnames(object$coefficients)
  group_idxs <- lapply(group_levels, function(g) which(object$group == g))
  lp <- cox_linear_predictor(object$x, object$group, object$coefficients)
  risk_set <- cox_risk_set(exp(lp), object$time, group_idxs)

  out <- lapply(seq_along(group_idxs), function(k) {
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
  do.call(rbind, out)
}

#' Variance-covariance matrix for a \code{coxmtl} object
#'
#' @param object An object of class \code{coxmtl}.
#' @param ... Additional arguments (unused).
#' @return A matrix representing the variance-covariance matrix of the free
#'   parameters.
#' @export
vcov.coxmtl <- function(object, ...) {
  phi <- coef(object)
  n_phi <- length(phi)
  if (n_phi == 0L) {
    empty <- matrix(0, nrow = 0L, ncol = 0L)
    dimnames(empty) <- list(character(0), character(0))
    return(empty)
  }

  group_levels <- colnames(object$coefficients)
  group_idxs <- lapply(group_levels, function(g) which(object$group == g))
  z <- dense_block_diag(lapply(group_idxs, function(idx) {
    object$x[idx, , drop = FALSE]
  })) %*% object$basis
  time <- unlist(lapply(group_idxs, function(idx) object$time[idx]))
  status <- unlist(lapply(group_idxs, function(idx) object$status[idx]))
  lp <- as.numeric(z %*% phi)

  sizes_v <- lengths(group_idxs)
  ends_v <- cumsum(sizes_v)
  starts_v <- c(1L, head(ends_v, -1L) + 1L)
  stacked_idxs <- Map(seq, starts_v, ends_v)

  gradients <- matrix(0, nrow = nrow(z), ncol = n_phi)
  hessians <- matrix(0, nrow = nrow(z), ncol = n_phi^2)
  for (k in seq_along(group_idxs)) {
    idx <- stacked_idxs[[k]]
    ghs <- calc_grad_hess(
      lp[idx], z[idx, , drop = FALSE], time[idx], status[idx]
    )
    gradients[idx, ] <- ghs$grad
    hessians[idx, ] <- ghs$hess
  }

  hess <- matrix(colSums(hessians), nrow = n_phi, ncol = n_phi)
  vc <- safe_solve(hess) %*% crossprod(gradients) %*% safe_solve(hess)
  dimnames(vc) <- list(names(phi), names(phi))
  vc
}

#' Summary method for a \code{coxmtl} object
#'
#' @param object An object of class \code{coxmtl}.
#' @param conf.int A numeric value between 0 and 1 indicating the confidence
#'   level of the confidence interval. Default is 0.95.
#' @param ... Additional arguments (unused).
#' @return An object of class \code{summary.coxmtl} with model size, event
#'   count, log-likelihood, coefficient table, confidence intervals, and center
#'   weights.
#' @export
summary.coxmtl <- function(object, conf.int = 0.95, ...) {
  beta_vec <- as.numeric(object$coefficients)
  names(beta_vec) <- rownames(object$basis)
  n_params <- length(beta_vec)
  n_events <- sum(object$status)
  loglik <- logLik(object)

  vcov_phi <- vcov(object)
  beta_vcov <- if (length(object$phi) == 0L) {
    matrix(0, nrow = n_params, ncol = n_params)
  } else {
    object$basis %*% vcov_phi %*% t(object$basis)
  }
  se <- sqrt(pmax(diag(beta_vcov), 0))
  z_scores <- ifelse(se > 0, beta_vec / se, 0)
  p_values <- ifelse(
    se > 0,
    stats::pchisq(z_scores^2, 1, lower.tail = FALSE),
    1
  )
  coef_matrix <- cbind(
    beta_vec, exp(beta_vec), se, z_scores, p_values
  )
  dimnames(coef_matrix) <- list(
    names(beta_vec), c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
  )

  z <- stats::qnorm((1 + conf.int) / 2)
  conf_int_matrix <- cbind(
    exp(beta_vec), exp(-beta_vec),
    exp(beta_vec - z * se), exp(beta_vec + z * se)
  )
  dimnames(conf_int_matrix) <- list(
    names(beta_vec), c(
      "exp(coef)", "exp(-coef)",
      stringr::str_c("lower .", round(100 * conf.int, 2), sep = ""),
      stringr::str_c("upper .", round(100 * conf.int, 2), sep = "")
    )
  )

  out <- list(
    n = nrow(object$x),
    nevent = n_events,
    logLik = loglik,
    df = length(object$phi),
    call = object$call,
    coefficients = coef_matrix,
    conf.int = conf_int_matrix,
    coefficient_matrix = object$coefficients,
    w = object$w
  )
  class(out) <- "summary.coxmtl"
  out
}

#' Print method for a \code{summary.coxmtl} object
#'
#' @param x A summary object produced from a fitted \code{coxmtl} model.
#' @param digits An integer controlling the number of significant digits to
#'   print for numeric values.
#' @param signif.stars Logical; if \code{TRUE}, significance stars are printed
#'   along with the p-values.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.summary.coxmtl <- function(
  x, digits = max(getOption("digits") - 3, 3),
  signif.stars = getOption("show.signif.stars"), ...
) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(
    "  n=", x$n, ", number of events=", x$nevent, ", df=", x$df, "\n\n",
    sep = ""
  )

  nonzero <- abs(x$coefficients[, "coef"]) > 1e-8
  if (any(nonzero)) {
    stats::printCoefmat(
      x$coefficients[nonzero, , drop = FALSE],
      digits = digits, signif.stars = signif.stars,
      cs.ind = 1:3, tst.ind = 4, P.values = TRUE, has.Pvalue = TRUE
    )
    print(
      format(x$conf.int[nonzero, , drop = FALSE], digits = digits),
      quote = FALSE
    )
  } else {
    cat("  No non-zero coefficients.\n")
  }

  cat("\nGlobal centers (w):\n")
  print(format(x$w, digits = digits), quote = FALSE)
  invisible(x)
}

#' Refit a coxmtl model with hard constraints (oracle estimator)
#'
#' Extracts the active constraint structure from a penalized \code{coxmtl} fit,
#' encodes it as hard constraints via the QR null-space basis, and fits an
#' unpenalized stratified Cox model in the reparametrized space. Returns
#' debiased MLE coefficients and valid standard errors for all groups.
#'
#' @param object A fitted \code{coxmtl} object.
#' @param ... Additional arguments (unused).
#'
#' @return An object of class \code{refit.coxmtl} with components:
#'   \item{coefficients}{Matrix (p x K) of oracle refitted coefficients.}
#'   \item{se}{Matrix (p x K) of standard errors (delta method).}
#'   \item{phi_hat}{Named numeric vector of free-parameter (phi) estimates.}
#'   \item{coxph_fit}{The underlying \code{coxph} object, or \code{NULL}.}
#'   \item{n}{Total observations used.}
#'   \item{nevent}{Number of events.}
#'   \item{active_features}{Integer vector: features non-zero in any group.}
#'   \item{basis}{The null-space basis \eqn{B} from the penalized fit.}
#'
#' @details
#' Oracle refit procedure:
#' \enumerate{
#'   \item Use the null-space basis \eqn{B} stored in the penalized fit, which
#'     encodes active sparsity, pairwise fusion, and center constraints.
#'   \item Build reparametrized design
#'     \eqn{Z = \mathrm{bdiag}(X_1,\ldots,X_K)\,B}.
#'   \item Fit \code{coxph(Surv(time,status) ~ Z + strata(group))} —
#'     unpenalized stratified Cox with all active constraints hard-coded in
#'     the design.
#'   \item Map back via \eqn{\hat\beta = B\hat\phi}; compute standard errors
#'     via the delta method.
#' }
#'
#' @rdname refit
#' @export
refit.coxmtl <- function(object, ...) {
  p <- nrow(object$coefficients)
  n_groups <- ncol(object$coefficients)
  feature_names <- rownames(object$coefficients)
  group_levels <- colnames(object$coefficients)
  n_phi <- ncol(object$basis)

  active_features <- which(apply(abs(object$coefficients) > 1e-8, 1, any))

  zero_mat <- matrix(0, p, n_groups,
    dimnames = list(feature_names, group_levels)
  )

  make_result <- function(coef_mat, se_mat = zero_mat,
                          coxph_fit = NULL, phi_hat = numeric(0),
                          n = 0L, nevent = 0L) {
    out <- list(
      coefficients = coef_mat,
      se = se_mat,
      phi_hat = phi_hat,
      coxph_fit = coxph_fit,
      n = n,
      nevent = nevent,
      active_features = active_features,
      basis = object$basis,
      group_levels = group_levels,
      feature_names = feature_names
    )
    class(out) <- "refit.coxmtl"
    out
  }

  if (n_phi == 0L) {
    return(make_result(zero_mat))
  }

  group_idxs <- lapply(group_levels, function(g) which(object$group == g))

  z_mat <- dense_block_diag(lapply(group_idxs, function(idx) {
    object$x[idx, , drop = FALSE]
  })) %*% object$basis
  colnames(z_mat) <- paste0("phi", seq_len(n_phi))

  time_s <- unlist(lapply(group_idxs, function(idx) object$time[idx]))
  status_s <- unlist(lapply(group_idxs, function(idx) object$status[idx]))
  strata_s <- factor(rep(seq_along(group_idxs), lengths(group_idxs)))

  df_refit <- as.data.frame(z_mat)
  df_refit$time <- time_s
  df_refit$status <- status_s
  df_refit$strata_var <- strata_s

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
    return(make_result(zero_mat))
  }

  phi_hat <- stats::coef(coxph_fit)
  beta_vec <- as.numeric(object$basis %*% phi_hat)
  coef_mat <- matrix(beta_vec,
    nrow = p, ncol = n_groups,
    dimnames = list(feature_names, group_levels)
  )

  vcov_phi <- stats::vcov(coxph_fit)
  se_mat <- zero_mat
  for (k in seq_len(n_groups)) {
    rows_k <- (k - 1L) * p + seq_len(p)
    l_k <- object$basis[rows_k, , drop = FALSE]
    se_mat[, k] <- sqrt(pmax(diag(l_k %*% vcov_phi %*% t(l_k)), 0))
  }

  make_result(coef_mat, se_mat, coxph_fit, phi_hat,
    n = nrow(df_refit), nevent = sum(status_s)
  )
}

#' Print method for a \code{refit.coxmtl} object
#'
#' @param x A \code{refit.coxmtl} object.
#' @param digits Number of significant digits.
#' @param signif.stars Logical; print significance stars.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.refit.coxmtl <- function(
  x, digits = max(getOption("digits") - 3, 3),
  signif.stars = getOption("show.signif.stars"), ...
) {
  cat("Refit coxmtl (unpenalized MLE, hard constraints)\n")
  cat("  n=", x$n, ", events=", x$nevent, "\n\n", sep = "")

  if (length(x$active_features) == 0L) {
    cat("  No active features.\n")
    return(invisible(x))
  }

  jj <- x$active_features
  for (k in seq_along(x$group_levels)) {
    cat("Group:", x$group_levels[k], "\n")
    beta_k <- x$coefficients[jj, k]
    se_k <- x$se[jj, k]
    z_k <- ifelse(se_k > 0, beta_k / se_k, 0)
    p_k <- ifelse(se_k > 0, stats::pchisq(z_k^2, 1, lower.tail = FALSE), 1)
    cm <- cbind(beta_k, exp(beta_k), se_k, z_k, p_k)
    colnames(cm) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
    rownames(cm) <- x$feature_names[jj]
    stats::printCoefmat(
      cm,
      digits = digits, signif.stars = signif.stars,
      cs.ind = 1:3, tst.ind = 4, P.values = TRUE, has.Pvalue = TRUE
    )
    cat("\n")
  }
  invisible(x)
}

#' Bayesian Information Criterion for \code{coxmtl} objects
#'
#' @param object A \code{coxmtl} object.
#' @param ... Additional arguments (unused).
#' @export
BIC.coxmtl <- function(object, ...) {
  nevent <- sum(object$status)
  if (nevent <= 0) {
    stop("BIC is undefined when the number of events is zero")
  }
  -2 * logLik(object) + log(nevent) * length(object$phi)
}
