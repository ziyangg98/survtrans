#' Cross-validated tuning for coxtrans
#'
#' Selects penalty parameters by minimising the held-out partial likelihood
#' deviance over a full \code{lambda1 x lambda2 x lambda3} grid. Supports both
#' the \code{lambda.min} rule (minimum CV deviance) and the \code{lambda.1se}
#' rule (most sparse model within one standard error of the minimum).
#'
#' @param formula A formula with a \code{\link[survival]{Surv}} response.
#' @param data A data frame.
#' @param group A factor indicating group membership.
#' @param target The target group level.
#' @param prior_matrix Optional G x (K-1) transfer constraint matrix.
#' @param nlambda Number of lambda values per dimension. Default 10.
#' @param lambda.min.ratio Smallest/largest lambda ratio. Default 1e-3.
#' @param nfolds Number of CV folds. Default 10.
#' @param penalty \code{"lasso"}, \code{"MCP"}, or \code{"SCAD"}.
#' @param ncores Number of R worker processes for parallel grid evaluation.
#'   Default 1.
#' @param nthreads Number of OpenMP threads per \code{coxtrans()} call.
#'   Default 1. Total CPU usage is \code{ncores x nthreads}.
#' @param seed Random seed.
#' @param verbose Print progress.
#' @param ... Passed to \code{\link{coxtrans}}.
#'
#' @return An object of class \code{cv.coxtrans} with fields:
#'   \code{grid} (full lambda grid), \code{cvm} (mean deviance per grid
#'   point), \code{cvsd}, \code{cvup}, \code{cvlo}, \code{nzero} (mean
#'   non-zero target coefficients), \code{lambda.min} (optimal lambda list),
#'   \code{lambda.1se} (1SE-rule lambda list), \code{index.min},
#'   \code{index.1se}, \code{coxtrans.fit} (final model at lambda.min),
#'   \code{coxtrans.fit.1se} (final model at lambda.1se),
#'   \code{call}, \code{name}.
#' @importFrom stats ave predict
#' @export
#'
#' @examples
#' result <- cv.coxtrans(Surv(time, status) ~ . - group - id,
#'   sim2, sim2$group,
#'   target = 1, nlambda = 5, penalty = "SCAD"
#' )
#' result
cv.coxtrans <- function(
  formula, data, group, target,
  prior_matrix = NULL, nlambda = 10, lambda.min.ratio = 1e-3,
  nfolds = 10, penalty = c("lasso", "MCP", "SCAD"),
  ncores = 1, nthreads = 1L, seed = 42, verbose = FALSE, ...
) {
  this_call <- match.call()
  penalty <- match.arg(penalty)
  if (!is.factor(group)) group <- factor(group)
  target_level <- as.character(target)
  n_priors <- if (is.null(prior_matrix)) 1L else nrow(prior_matrix)

  lpath <- function(lmax) {
    if (!is.finite(lmax) || lmax <= 0) {
      return(0)
    }
    exp(seq(log(lmax * lambda.min.ratio), log(lmax), length.out = nlambda))
  }

  l3_cols <- paste0("lambda3_", seq_len(n_priors))
  grid <- do.call(expand.grid, c(
    list(
      lambda1 = lpath(calc_lambda1_max(formula, data, group, target)),
      lambda2 = lpath(calc_lambda2_max(formula, data, group, target))
    ),
    stats::setNames(
      lapply(
        calc_lambda3_max(formula, data, group, target, prior_matrix),
        lpath
      ),
      l3_cols
    )
  ))

  # Stratified fold assignment on target observations
  target_idx <- which(group == target_level)
  set.seed(seed)
  y_target <- stats::model.response(
    stats::model.frame(formula, data[target_idx, ])
  )
  fids <- ave(
    seq_along(target_idx), y_target[, 2],
    FUN = function(i) sample(rep(seq_len(nfolds), length.out = length(i)))
  )

  if (verbose) {
    cat(sprintf("CV: %d-fold over %d grid points\n", nfolds, nrow(grid)))
  }

  eval_point <- function(i) {
    l1 <- grid$lambda1[i]
    l2 <- grid$lambda2[i]
    l3 <- unlist(grid[i, l3_cols, drop = TRUE])
    fold_results <- lapply(seq_len(nfolds), function(k) {
      ti <- target_idx[fids == k]
      fit <- tryCatch(
        coxtrans(formula, data[-ti, ], group[-ti], target,
          lambda1 = l1, lambda2 = l2, lambda3 = l3,
          prior_matrix = prior_matrix, penalty = penalty,
          nthreads = nthreads, ...
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) return(list(dev = NA_real_, nzero = NA_real_))
      yt <- stats::model.response(stats::model.frame(formula, data[ti, ]))
      lp <- predict(fit, newdata = data[ti, ], newgroup = group[ti])
      n_events <- sum(yt[, 2])
      if (n_events == 0L) return(list(dev = NA_real_, nzero = NA_real_))
      o        <- order(yt[, 1], decreasing = TRUE)
      lp_o     <- lp[o]
      max_lp   <- max(lp_o)
      cs       <- cumsum(exp(lp_o - max_lp))
      log_risk <- max_lp + log(ave(cs, yt[o, 1], FUN = max))
      dev   <- -2 * sum(yt[o, 2] * (lp_o - log_risk)) / n_events
      nzero <- sum(abs(fit$coefficients[, 1]) > 0)
      list(dev = dev, nzero = as.numeric(nzero))
    })
    list(
      devs  = vapply(fold_results, `[[`, numeric(1), "dev"),
      nzero = mean(
        vapply(fold_results, `[[`, numeric(1), "nzero"),
        na.rm = TRUE
      )
    )
  }

  n_grid <- nrow(grid)
  results <- if (ncores > 1) {
    parallel::mclapply(seq_len(n_grid), eval_point, mc.cores = ncores)
  } else {
    lapply(seq_len(n_grid), eval_point)
  }

  devs_list <- lapply(results, `[[`, "devs")
  cvm <- vapply(devs_list, function(x) mean(x, na.rm = TRUE), numeric(1))
  cvsd <- vapply(devs_list, function(x) {
    x <- x[is.finite(x)]
    if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else NA_real_
  }, numeric(1))
  nzero <- as.integer(round(
    vapply(results, `[[`, numeric(1), "nzero")
  ))

  if (!any(is.finite(cvm))) {
    stop("All CV scores are non-finite. Check data or fold sizes.")
  }

  i_best <- which.min(cvm)
  lambda_min <- list(
    lambda1 = grid$lambda1[i_best],
    lambda2 = grid$lambda2[i_best],
    lambda3 = unlist(grid[i_best, l3_cols, drop = TRUE])
  )

  # 1SE rule: most sparse model within 1 SE of minimum deviance
  # Only select if strictly sparser than lambda.min; otherwise degenerate
  threshold_1se <- cvm[i_best] + cvsd[i_best]
  candidates_1se <- which(
    is.finite(cvm) & cvm <= threshold_1se & nzero < nzero[i_best]
  )
  if (length(candidates_1se) == 0L) {
    i_1se <- i_best
  } else {
    i_1se <- candidates_1se[which.min(nzero[candidates_1se])]
  }
  lambda_1se <- list(
    lambda1 = grid$lambda1[i_1se],
    lambda2 = grid$lambda2[i_1se],
    lambda3 = unlist(grid[i_1se, l3_cols, drop = TRUE])
  )

  if (verbose) {
    cat(sprintf(
      "lambda.min: l1=%.4f  l2=%.4f  l3=%s  (deviance=%.4f)\n",
      lambda_min$lambda1, lambda_min$lambda2,
      paste(round(lambda_min$lambda3, 4), collapse = ","),
      cvm[i_best]
    ))
    cat(sprintf(
      "lambda.1se: l1=%.4f  l2=%.4f  l3=%s  (deviance=%.4f, nzero=%d)\n",
      lambda_1se$lambda1, lambda_1se$lambda2,
      paste(round(lambda_1se$lambda3, 4), collapse = ","),
      cvm[i_1se], nzero[i_1se]
    ))
  }

  final_fit <- coxtrans(formula, data, group, target,
    lambda1 = lambda_min$lambda1,
    lambda2 = lambda_min$lambda2,
    lambda3 = lambda_min$lambda3,
    prior_matrix = prior_matrix, penalty = penalty,
    nthreads = nthreads, ...
  )

  final_fit_1se <- coxtrans(formula, data, group, target,
    lambda1 = lambda_1se$lambda1,
    lambda2 = lambda_1se$lambda2,
    lambda3 = lambda_1se$lambda3,
    prior_matrix = prior_matrix, penalty = penalty,
    nthreads = nthreads, ...
  )

  structure(list(
    grid           = grid,
    cvm            = cvm,
    cvsd           = cvsd,
    cvup           = cvm + cvsd,
    cvlo           = cvm - cvsd,
    nzero          = nzero,
    lambda.min     = lambda_min,
    lambda.1se     = lambda_1se,
    index.min      = i_best,
    index.1se      = i_1se,
    index          = i_best, # backward compat
    coxtrans.fit     = final_fit,
    coxtrans.fit.1se = final_fit_1se,
    call           = this_call,
    name           = c("Partial Likelihood Deviance" = "deviance")
  ), class = "cv.coxtrans")
}

#' @export
print.cv.coxtrans <- function(x, ...) {
  cat("cv.coxtrans\n\n")
  cat("Call: ")
  print(x$call)
  cat(sprintf("\nMeasure: %s\n\n", names(x$name)))

  i_min <- if (!is.null(x$index.min)) x$index.min else x$index
  i_1se <- x$index.1se

  cat(sprintf(
    "lambda.min:  l1=%.4f  l2=%.4f  l3=%s\n",
    x$lambda.min$lambda1, x$lambda.min$lambda2,
    paste(round(x$lambda.min$lambda3, 4), collapse = ",")
  ))
  cat(sprintf(
    "  Deviance:  %.4f (+/- %.4f)   Non-zero: %d\n",
    x$cvm[i_min], x$cvsd[i_min], x$nzero[i_min]
  ))

  if (!is.null(i_1se)) {
    cat(sprintf(
      "\nlambda.1se:  l1=%.4f  l2=%.4f  l3=%s\n",
      x$lambda.1se$lambda1, x$lambda.1se$lambda2,
      paste(round(x$lambda.1se$lambda3, 4), collapse = ",")
    ))
    cat(sprintf(
      "  Deviance:  %.4f (+/- %.4f)   Non-zero: %d\n",
      x$cvm[i_1se], x$cvsd[i_1se], x$nzero[i_1se]
    ))
  }

  invisible(x)
}

#' Plot cross-validation curve for a \code{cv.coxtrans} object
#'
#' Plots the CV deviance profile along \code{lambda1} with \code{lambda2} and
#' \code{lambda3} fixed at their optimal values, in the style of
#' \code{\link[glmnet]{plot.cv.glmnet}}.
#'
#' @param x A \code{cv.coxtrans} object.
#' @param ... Further graphical arguments passed to \code{plot}.
#'
#' @importFrom graphics abline axis legend mtext par plot points segments
#' @export
plot.cv.coxtrans <- function(x, ...) {
  l3_cols <- grep("^lambda3_", names(x$grid), value = TRUE)

  # Profile: fix lambda2 and lambda3 at optimal values, vary lambda1
  l2_opt <- x$lambda.min$lambda2
  l3_opt <- x$lambda.min$lambda3
  l3_match <- if (length(l3_cols) > 0) {
    apply(x$grid[, l3_cols, drop = FALSE], 1, function(r) {
      all(abs(unlist(r) - l3_opt) < .Machine$double.eps * 1e6)
    })
  } else {
    rep(TRUE, nrow(x$grid))
  }
  idx <- which(
    abs(x$grid$lambda2 - l2_opt) < .Machine$double.eps * 1e6 & l3_match
  )
  if (length(idx) == 0) idx <- seq_len(nrow(x$grid))

  ord  <- order(x$grid$lambda1[idx])
  idx  <- idx[ord]
  l1   <- x$grid$lambda1[idx]
  cvm  <- x$cvm[idx]
  cvlo <- x$cvlo[idx]
  cvup <- x$cvup[idx]
  nz   <- x$nzero[idx]

  xv   <- log(pmax(l1, .Machine$double.eps))
  ylim <- range(c(cvlo, cvup), na.rm = TRUE)

  old_par <- par(mar = c(4.5, 4.5, 5, 1))
  on.exit(par(old_par))

  plot(xv, cvm, type = "n", ylim = ylim,
       xlab = expression(log(lambda[1])),
       ylab = names(x$name), ...)
  segments(xv, cvlo, xv, cvup, col = "grey60")
  points(xv, cvm, pch = 20, col = "red")
  v_min <- log(max(x$lambda.min$lambda1, .Machine$double.eps))
  i_min <- if (!is.null(x$index.min)) x$index.min else x$index
  degenerate <- is.null(x$index.1se) ||
    identical(i_min, x$index.1se)

  abline(v = v_min, lty = 2, col = "darkgrey")
  if (!degenerate) {
    threshold <- x$cvm[i_min] + x$cvsd[i_min]
    abline(h = threshold, lty = 3, col = "darkgrey")
    legend("topright",
      legend = c("lambda.min", "1 SE threshold"),
      lty = c(2, 3), col = "darkgrey", bty = "n", cex = 0.8
    )
  }

  axis(3, at = xv, labels = nz, tick = FALSE, las = 2, cex.axis = 0.7)
  mtext("Non-zero coefficients", side = 3, line = 3, cex = 0.8)

  invisible(x)
}

#' @param object A \code{cv.coxtrans} object.
#' @param s Which model to extract: \code{"lambda.min"} (default) or
#'   \code{"lambda.1se"}.
#' @rdname cv.coxtrans
#' @export
coef.cv.coxtrans <- function(object, s = c("lambda.min", "lambda.1se"), ...) {
  s <- match.arg(s)
  fit <- if (s == "lambda.1se" && !is.null(object$coxtrans.fit.1se)) {
    object$coxtrans.fit.1se
  } else {
    object$coxtrans.fit
  }
  coef(fit, ...)
}

#' @param s Which model to predict from: \code{"lambda.min"} (default) or
#'   \code{"lambda.1se"}.
#' @rdname cv.coxtrans
#' @export
predict.cv.coxtrans <- function(object,
                                s = c("lambda.min", "lambda.1se"), ...) {
  s <- match.arg(s)
  fit <- if (s == "lambda.1se" && !is.null(object$coxtrans.fit.1se)) {
    object$coxtrans.fit.1se
  } else {
    object$coxtrans.fit
  }
  predict(fit, ...)
}
