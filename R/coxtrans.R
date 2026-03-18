#' Transfer Learning for Cox Model with Global and Local Shrinkage
#'
#' @param formula A formula object, with the response on the left of a \code{~}
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function.
#' @param data A data frame containing the variables in the model.
#' @param group A factor variable indicating the group of each observation.
#' @param target A character string specifying the target group.
#' @param lambda1 A non-negative value specifying the sparse penalty parameter.
#' The default is 0.
#' @param lambda2 A non-negative value specifying the global biased penalty
#' parameter. The default is 0.
#' @param lambda3 A non-negative value specifying the local biased penalty
#' parameter. The default is 0.
#' @param penalty A character string specifying the penalty function. The
#' default is "lasso". Other options are "MCP" and "SCAD".
#' @param gamma A non-negative value specifying the penalty parameter. The
#' default is 3.7 for SCAD and 3.0 for MCP.
#' @param vartheta A positive value specifying the fixed penalty parameter in
#' the augmented Lagrangian. Following Wang, Yin & Zeng (2019), this is kept
#' constant (not adaptive) to guarantee convergence with non-convex penalties.
#' The default is 1.0.
#' @param control An object of class \link{survtrans_control} containing control
#' parameters for the fitting algorithm. Default is
#' \code{survtrans_control(...)}.
#' @param ... Additional arguments to be passed to the fitting algorithm.
#'
#' @return An object of class \code{coxtrans}.
#' @export
#'
#' @examples
#' formula <- Surv(time, status) ~ . - group - id
#' fit <- coxtrans(
#'   formula, sim2, sim2$group, 1,
#'   lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
#' )
#' summary(fit)
coxtrans <- function(
  formula, data, group, target,
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
  if (lambda1 < 0 || lambda2 < 0 || lambda3 < 0) {
    stop("Lambda parameters must be non-negative")
  }

  # Preprocess data
  data <- preprocess(formula, data, group = group)
  x <- data$x
  x_scale <- attr(x, "scaled:scale")
  time <- data$time
  status <- data$status
  group <- data$group

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
  n_samples_group <- sapply(group_idxs, length)
  n_samples_total <- nrow(x)
  n_features <- ncol(x)
  n_parameters <- n_features * (n_groups + 1)

  # Construct constraint matrices
  contr_sum <- cbind(
    kronecker(
      matrix(n_samples_group / n_samples_total, nrow = 1),
      diag(n_features)
    ),
    0 * diag(n_features)
  )
  contr_pen <- matrix(0, n_parameters, n_parameters)
  contr_pen[cbind(
    c(
      seq_len(n_features), seq_len(n_features),
      n_features + seq_len(n_features * n_groups),
      n_features + seq_len(n_features * (n_groups - 1))
    ),
    c(
      seq_len(n_features), (n_features * n_groups) + seq_len(n_features),
      rep(1:n_features, times = n_groups),
      n_features + seq_len(n_features * (n_groups - 1))
    )
  )] <- c(
    rep(1, n_features * (n_groups + 2)), rep(-1, n_features * (n_groups - 1))
  )
  contr_cross <- crossprod(contr_sum) + crossprod(contr_pen)
  n_constraints <- nrow(contr_sum) + nrow(contr_pen)
  n_constraints_penalty <- nrow(contr_pen)

  sparse_idx <- seq_len(n_features)
  local_idx <- n_features + seq_len(n_features * (n_groups - 1))
  global_idx <- n_features * n_groups + seq_len(n_features)

  # Per-group data (stacked by group order)
  x_by_group <- lapply(group_idxs, function(idx) x[idx, , drop = FALSE])
  time_stacked <- unlist(lapply(group_idxs, function(idx) time[idx]))
  status_stacked <- unlist(lapply(group_idxs, function(idx) status[idx]))

  stacked_group_idxs <- vector("list", n_groups)
  n_passes <- 0L
  for (k in seq_len(n_groups)) {
    nk <- length(group_idxs[[k]])
    stacked_group_idxs[[k]] <- n_passes + seq_len(nk)
    n_passes <- n_passes + nk
  }

  # Initialize ADMM variables
  theta <- matrix(0, nrow = n_parameters, ncol = 1)
  eta <- numeric(n_constraints_penalty)
  mu <- numeric(n_features)
  nu <- numeric(n_constraints_penalty)

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
  lag_aug_prev <- Inf
  loss_total_prev <- Inf

  # Block index ranges for assembling xwx/xwz
  blk_idx <- lapply(seq_len(n_groups), function(k) {
    ((k - 1) * n_features + 1):(k * n_features)
  })
  blk_center <- (n_groups * n_features + 1):((n_groups + 1) * n_features)

  repeat {
    n_iterations <- n_iterations + 1

    # Compute per-group offset from theta
    offset <- calc_offset(
      theta, n_features, n_groups, x_by_group, stacked_group_idxs
    )

    # Calculate the weights and working response (IRLS)
    for (k in seq_len(n_groups)) {
      idx <- stacked_group_idxs[[k]]
      wls <- approx_likelihood(
        offset = offset[idx], time = time_stacked[idx], status = status_stacked[idx]
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
      A_k <- crossprod(x_by_group[[k]], w_k * x_by_group[[k]]) / n_samples_total
      b_k <- crossprod(x_by_group[[k]], w_k * z_k) / n_samples_total
      xwx[blk_idx[[k]], blk_idx[[k]]] <- A_k
      xwx[blk_idx[[k]], blk_center] <- A_k
      xwx[blk_center, blk_idx[[k]]] <- A_k
      xwx[blk_center, blk_center] <- xwx[blk_center, blk_center] + A_k
      xwz[blk_idx[[k]]] <- b_k
      xwz[blk_center] <- xwz[blk_center] + b_k
    }

    # Solve the linear system
    lhs <- xwx + vartheta * contr_cross
    rhs <- xwz - crossprod(contr_sum, mu) +
      vartheta * crossprod(contr_pen, eta - nu / vartheta)
    theta <- solve(lhs, rhs)

    # Constraint products
    c_theta <- as.numeric(contr_pen %*% theta)
    s_theta <- as.numeric(contr_sum %*% theta)

    # Update auxiliary variables
    eta_old <- eta
    eta <- c_theta + nu / vartheta
    eta[sparse_idx] <- threshold_prox(
      eta[sparse_idx], vartheta, penalty, lambda1, gamma
    )
    eta[global_idx] <- threshold_prox(
      eta[global_idx], vartheta, penalty, lambda2, gamma
    )
    eta[local_idx] <- threshold_prox(
      eta[local_idx], vartheta, penalty, lambda3, gamma
    )
    mu <- mu + vartheta * s_theta
    nu <- nu + vartheta * (c_theta - eta)

    # Primal and dual residuals
    r_norm <- sqrt(sum(s_theta^2) + sum((c_theta - eta)^2))
    s_norm <- sqrt(sum(crossprod(contr_pen, eta - eta_old)^2)) *
      vartheta

    eps_pri <- sqrt(n_constraints) * control$abstol +
      control$reltol * max(
        sqrt(sum(s_theta^2) + sum(c_theta^2)),
        sqrt(sum(eta^2))
      )
    dual_vec <- crossprod(contr_sum, mu) +
      crossprod(contr_pen, nu)
    eps_dual <- sqrt(n_parameters) * control$abstol +
      control$reltol * sqrt(sum(dual_vec^2))

    # Compute loss (recompute offset from updated theta)
    offset <- calc_offset(
      theta, n_features, n_groups, x_by_group, stacked_group_idxs
    )
    hazard <- exp(offset)
    risk_set <- calc_risk_set(hazard, time_stacked, stacked_group_idxs)
    loss <- sum(status_stacked * (offset - log(risk_set)))

    loss_penalty <- (
      penalty_value(c_theta[sparse_idx], penalty, lambda1, gamma) +
        penalty_value(c_theta[global_idx], penalty, lambda2, gamma) +
        penalty_value(c_theta[local_idx], penalty, lambda3, gamma)
    ) * n_samples_total
    loss_total <- loss - loss_penalty

    # Augmented Lagrangian (Lyapunov function for convergence monitoring)
    lag_aug <- loss_total +
      sum(mu * s_theta) + sum(nu * (c_theta - eta)) +
      0.5 * vartheta * (sum(s_theta^2) + sum((c_theta - eta)^2))

    # Check convergence
    if (r_norm < eps_pri && s_norm < eps_dual) {
      converged <- TRUE
      msg <- sprintf("Convergence reached at iteration %d.", n_iterations)
    } else if (is.infinite(loss) || is.nan(loss)) {
      converged <- TRUE
      msg <- "Log-likelihood is not finite. Stopping."
    } else if (n_iterations > 1 &&
      (abs(lag_aug - lag_aug_prev) / (abs(lag_aug_prev) + 1) < control$fdev ||
        abs(loss_total - loss_total_prev) / (abs(loss_total_prev) + 1) < control$fdev)) {
      converged <- TRUE
      msg <- sprintf("Objective stabilized at iteration %d.", n_iterations)
    } else if (n_iterations >= control$maxit) {
      converged <- TRUE
      msg <- sprintf("Maximum number of iterations reached (%d).", control$maxit)
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

  # Post-processing: recover coefficients from ADMM solution
  theta <- qr.solve(contr_pen, eta)
  tol <- max(abs(crossprod(contr_pen, eta - eta_old)))
  tol_local <- max(abs(eta[local_idx] - eta_old[local_idx]))
  flag_local <- matrix(abs(eta[local_idx]) <= tol_local, nrow = n_features)
  flag_local <- cbind(rep(TRUE, n_features), flag_local)

  theta <- matrix(theta, nrow = n_features, ncol = n_groups + 1)
  delta <- theta[, seq_len(n_groups)]
  center <- theta[, n_groups + 1]
  beta <- matrix(NA, nrow = n_features, ncol = n_groups)

  for (i in seq_len(n_features)) {
    delta_global <- mean(delta[i, ])
    idx <- which(flag_local[i, ])
    delta_local <- mean(delta[i, idx])
    delta[i, idx] <- ifelse(
      abs(delta_local - delta_global) < tol, 0, delta_local - delta_global
    )
    beta[i, ] <- ifelse(
      abs(delta[i, ] + center[i]) < tol, 0, delta[i, ] + center[i]
    )
  }
  center <- rowMeans(beta)

  coefficients <- cbind(beta - center, center)
  coefficients[abs(coefficients) < tol] <- 0
  colnames(coefficients) <- c(group_levels, "Center")
  rownames(coefficients) <- colnames(x)

  coefficients <- sweep(coefficients, 1, x_scale, "/")
  x <- sweep(x, 2, x_scale, "*")

  # Return the fitted model
  fit <- list(
    coefficients = coefficients,
    logLik = loss,
    iter = n_iterations,
    message = msg,
    history = history[seq_len(n_iterations), , drop = FALSE],
    penalty = penalty,
    lambda1 = lambda1,
    lambda2 = lambda2,
    lambda3 = lambda3,
    gamma = gamma,
    formula = formula,
    call = match.call(),
    time = time,
    status = status,
    group = group,
    x = x,
    control = control
  )
  class(fit) <- "coxtrans"
  fit
}

# Internal: compute per-group beta matrix and linear predictor from a
# coxtrans fit object.
calc_lp <- function(object) {
  coefficients <- object$coefficients
  n_groups <- ncol(coefficients) - 1L
  group_levels <- colnames(coefficients)[seq_len(n_groups)]
  group_idxs <- lapply(group_levels, function(g) which(object$group == g))
  beta <- coefficients[, seq_len(n_groups)] + coefficients[, n_groups + 1L]

  lp <- numeric(nrow(object$x))
  for (k in seq_len(n_groups)) {
    idx <- group_idxs[[k]]
    lp[idx] <- object$x[idx, ] %*% beta[, k]
  }
  list(beta = beta, lp = lp, group_idxs = group_idxs, group_levels = group_levels)
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

#' Extract the coefficients from a \code{coxtrans} object
#' @param object An object of class \code{coxtrans}.
#' @param ... Additional arguments (unused).
#' @return A named numeric vector containing the coefficients of the fitted
#' \code{coxtrans} object. The names indicate the group(s) to which the
#' coefficients belong. Zero coefficients are removed.
#' @export
coef.coxtrans <- function(object, ...) {
  coefficients <- object$coefficients
  n_groups <- ncol(coefficients) - 1L
  n_features <- nrow(coefficients)
  beta <- coefficients[, 1:n_groups, drop = FALSE] +
    coefficients[, n_groups + 1L]

  phi_list <- lapply(
    seq_len(n_features),
    function(j) unique(beta[j, ])
  )

  is_global <- abs(coefficients[, 1]) < .Machine$double.eps^0.5
  psi_list <- lapply(seq_len(n_features), function(j) {
    vals <- phi_list[[j]]
    if (is_global[j] && length(vals) > 1) vals[-1] else vals
  })
  psi <- unlist(psi_list)
  names(psi) <- unlist(lapply(seq_along(psi_list), function(j) {
    vals <- psi_list[[j]]
    if (is_global[j] && length(vals) > 1) {
      idx <- seq_along(vals) + 1
    } else {
      idx <- seq_along(vals)
    }
    stringr::str_c(rownames(coefficients)[j], idx, sep = ".")
  }))
  psi
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
  n_groups <- ncol(coefficients) - 1L
  # Use coefficients column order (target-first) to match link_matrix
  coef_groups <- colnames(coefficients)[seq_len(n_groups)]
  group_idxs <- lapply(coef_groups, function(g) which(group == g))

  psi <- coef(object)
  link_matrix <- build_link_matrix(coefficients)

  z <- block_diag(lapply(group_idxs, function(idx) x[idx, ])) %*% link_matrix
  time <- unlist(lapply(group_idxs, function(idx) time[idx]))
  status <- unlist(lapply(group_idxs, function(idx) status[idx]))

  is_nonzero <- as.vector(psi) != 0
  n_nonzero <- sum(is_nonzero)
  z1 <- as.matrix(z[, is_nonzero, drop = FALSE])
  psi1 <- psi[is_nonzero]
  lp <- z1 %*% psi1

  gradients <- matrix(0, nrow = n_samples, ncol = n_nonzero)
  hessians <- matrix(0, nrow = n_samples, ncol = n_nonzero^2)
  n_passes <- 0
  for (k in seq_len(n_groups)) {
    idx <- n_passes + seq_len(length(group_idxs[[k]]))
    n_passes <- n_passes + length(idx)
    ghs <- calc_grad_hess(lp[idx], z1[idx, ], time[idx], status[idx])
    gradients[idx, ] <- ghs$grad
    hessians[idx, ] <- ghs$hess
  }
  hess <- matrix(colSums(hessians), n_nonzero, n_nonzero)
  hess_inv <- solve(hess)
  grad_cov <- crossprod(gradients)
  vcov <- hess_inv %*% grad_cov %*% hess_inv
  dimnames(vcov) <- list(
    names(psi1), names(psi1)
  )
  vcov
}


#' Log-likelihood for a \code{coxtrans} object
#'
#' @param object An object of class \code{coxtrans}.
#' @param ... Additional arguments (unused).
#' @return A numeric value representing the log-likelihood of the fitted
#' \code{coxtrans} object.
#' @export
logLik.coxtrans <- function(object, ...) {
  res <- calc_lp(object)
  hazard <- exp(res$lp)
  risk_set <- calc_risk_set(hazard, object$time, res$group_idxs)
  sum(object$status * (res$lp - log(risk_set)))
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
    link_matrix <- as.matrix(build_link_matrix(object$coefficients)[
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

  # Create a summary list
  summary_list <- list(
    n = n_samples,
    nevent = n_events,
    logLik = loglik,
    call = object$call,
    coefficients = coef_matrix,
    conf.int = conf_int_matrix
  )

  class(summary_list) <- "summary.coxtrans"
  return(summary_list)
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

  invisible(x)
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

  coefficients <- object$coefficients
  n_groups <- ncol(coefficients) - 1L
  group_levels <- colnames(coefficients)[seq_len(n_groups)]
  beta <- coefficients[, seq_len(n_groups)] + coefficients[, n_groups + 1L]

  x <- stats::model.matrix(object$formula, newdata)[, -1]
  group <- factor(newgroup, levels = levels(object$group))

  lp <- numeric(nrow(x))
  for (k in seq_len(n_groups)) {
    idx <- which(group == group_levels[k])
    if (length(idx) > 0) {
      lp[idx] <- x[idx, ] %*% beta[, k]
    }
  }

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
  res <- calc_lp(object)
  n_groups <- length(res$group_idxs)
  hazard <- exp(res$lp)
  risk_set <- calc_risk_set(hazard, object$time, res$group_idxs)

  basehaz_list <- vector("list", n_groups)
  for (k in seq_len(n_groups)) {
    idx <- res$group_idxs[[k]]
    time_rev <- rev(object$time[idx])
    status_rev <- rev(object$status[idx])
    risk_set_rev <- rev(risk_set[idx])
    bh <- cumsum(status_rev / risk_set_rev)
    basehaz_list[[k]] <- data.frame(
      time = time_rev[status_rev == 1],
      basehaz = bh[status_rev == 1],
      strata = res$group_levels[k]
    )
  }
  do.call(rbind, basehaz_list)
}
