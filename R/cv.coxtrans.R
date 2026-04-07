#' Cross-validated coxtrans
#'
#' Selects penalty parameters (lambda1, lambda2, lambda3) via K-fold
#' cross-validation with coordinate descent search, then applies the 1-SE
#' rule to choose the simplest model with competitive prediction.
#'
#' @param formula A formula with a \code{\link[survival]{Surv}} response.
#' @param data A data frame containing the variables in the model.
#' @param group A factor indicating the group of each observation.
#' @param target The target group level.
#' @param prior_matrix Optional G x (K-1) transfer constraint matrix.
#'   See \code{\link{coxtrans}}.
#' @param nlambda Number of lambda values per dimension. Default is 100.
#' @param lambda.min.ratio Ratio of the smallest to largest lambda value
#'   in the search path. Default is \code{1e-3}, meaning the path goes
#'   from \code{lambda_max * lambda.min.ratio} to \code{lambda_max}.
#' @param nfolds Number of cross-validation folds. Default is 10.
#' @param seed Random seed for fold assignment. Default is 42.
#' @param penalty Penalty type: \code{"lasso"}, \code{"MCP"}, or
#'   \code{"SCAD"}.
#' @param ncores Number of CPU cores for parallel fold evaluation.
#' @param verbose Print progress if \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{coxtrans}}.
#'
#' @return A list with components:
#'   \item{best_lambdas}{List with lambda1, lambda2, lambda3.}
#'   \item{cv_results}{Per-dimension data frames with columns
#'     \code{lambda}, \code{loglik}, \code{se}.}
#'   \item{fit}{The final coxtrans fit on full data.}
#'
#' @details
#' \strong{Search strategy:} Coordinate descent over lambda dimensions in
#' order lambda3 (prior) -> lambda2 (local) -> lambda1 (sparse), repeated
#' up to 3 cycles until convergence.
#'
#' \strong{Model selection:} Among all visited lambda combinations with CV
#' partial log-likelihood within one standard error of the maximum
#' (Breiman et al., 1984), the model with the smallest target effective
#' degrees of freedom is chosen.
#'
#' \strong{Target effective df} (Hodges & Sargent, 2001): For each non-zero
#' target feature, the df contribution depends on how the coefficient is
#' estimated:
#' \itemize{
#'   \item Sparse (zero): df = 0
#'   \item Prior constrained: df = n_target / (n_target + n_eff), where
#'     n_eff = 1 / sum(w_i^2 / n_source_i) is the effective sample size
#'     of the prior-weighted source combination
#'   \item Shared across groups: df = n_target / n_shared
#'   \item Group-specific: df = 1
#' }
#'
#' @export
#'
#' @examples
#' formula <- Surv(time, status) ~ . - group - id
#' cv_fit <- cv.coxtrans(formula, sim2, sim2$group,
#'   target = 1, nlambda = 5, penalty = "SCAD"
#' )
#' cv_fit$best_lambdas
#' summary(cv_fit$fit)
cv.coxtrans <- function(
  formula, data, group, target,
  prior_matrix = NULL,
  nlambda = 100,
  lambda.min.ratio = 1e-3,
  nfolds = 10,
  seed = 42,
  penalty = c("lasso", "MCP", "SCAD"),
  ncores = 1,
  verbose = FALSE,
  ...
) {
  penalty <- match.arg(penalty)
  if (!is.factor(group)) group <- factor(group)
  target_level <- as.character(target)
  target_idx <- which(group == target_level)
  n_priors <- if (is.null(prior_matrix)) 1L else nrow(prior_matrix)

  # Log-scale search path from lambda_max (like glmnet)
  lambda_max <- calc_lambda_max(formula, data, group)
  lambda_path <- exp(seq(
    log(lambda_max * lambda.min.ratio),
    log(lambda_max),
    length.out = nlambda
  ))

  # Stratified CV folds on target group (by event status)
  set.seed(seed)
  target_mf <- stats::model.frame(formula, data[target_idx, ])
  target_status <- stats::model.response(target_mf)[, 2]
  fold_ids <- integer(length(target_idx))
  for (s in unique(target_status)) {
    idx_s <- which(target_status == s)
    fold_ids[idx_s] <- sample(rep(seq_len(nfolds), length.out = length(idx_s)))
  }
  fold_data <- lapply(seq_len(nfolds), function(k) {
    test_idx <- target_idx[fold_ids == k]
    train_rows <- setdiff(seq_len(nrow(data)), test_idx)
    mf <- stats::model.frame(formula, data[test_idx, ])
    y <- stats::model.response(mf)
    list(
      train_rows = train_rows,
      x_test = stats::model.matrix(formula, data[test_idx, ])[, -1],
      time_test = y[, 1], status_test = y[, 2]
    )
  })

  # Evaluate CV log-likelihood across folds
  cv_loglik_folds <- function(lambda1, lambda2, lambda3) {
    eval_fold <- function(k) {
      fd <- fold_data[[k]]
      fit <- tryCatch(
        coxtrans(formula, data[fd$train_rows, ], group[fd$train_rows],
          target,
          lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3,
          prior_matrix = prior_matrix, penalty = penalty, ...
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(-Inf)
      }
      lp <- as.numeric(fd$x_test %*% fit$coefficients[, 1])
      ord <- order(fd$time_test, decreasing = TRUE)
      sum(fd$status_test[ord] * (lp[ord] - log(cumsum(exp(lp[ord])))))
    }
    if (ncores > 1) {
      unlist(parallel::mclapply(seq_len(nfolds), eval_fold, mc.cores = ncores))
    } else {
      vapply(seq_len(nfolds), eval_fold, numeric(1))
    }
  }

  # Search one lambda dimension
  search_dim <- function(path, make_lambdas) {
    loglik_mean <- loglik_se <- numeric(length(path))
    for (i in seq_along(path)) {
      lambdas <- make_lambdas(path[i])
      fold_logliks <- cv_loglik_folds(lambdas[1], lambdas[2], lambdas[-(1:2)])
      loglik_mean[i] <- mean(fold_logliks)
      loglik_se[i] <- stats::sd(fold_logliks) / sqrt(nfolds)
    }
    list(
      best = path[which.max(loglik_mean)],
      loglik = loglik_mean, se = loglik_se
    )
  }

  # Record visited lambda combinations for 1-SE selection
  record_visited <- function(visited, res, lambda1, lambda2, lambda3,
                             dim_idx, path) {
    for (i in seq_along(path)) {
      lam1 <- lambda1
      lam2 <- lambda2
      lam3 <- lambda3
      if (dim_idx == 1) {
        lam1 <- path[i]
      } else if (dim_idx == 2) {
        lam2 <- path[i]
      } else {
        lam3[dim_idx - 2] <- path[i]
      }
      row <- data.frame(
        lambda1 = lam1, lambda2 = lam2,
        loglik_mean = res$loglik[i], loglik_se = res$se[i]
      )
      for (g in seq_len(n_priors)) row[[paste0("lambda3_", g)]] <- lam3[g]
      visited <- rbind(visited, row)
    }
    visited
  }

  # Initialize at zero (no penalty)
  best_lambda1 <- 0
  best_lambda2 <- 0
  best_lambda3 <- rep(0, n_priors)

  dim_names <- c(
    "lambda1", "lambda2",
    if (n_priors == 1) "lambda3" else paste0("lambda3[", seq_len(n_priors), "]")
  )
  n_dims <- 2 + n_priors
  cv_results <- stats::setNames(vector("list", n_dims), dim_names)

  visited <- data.frame(
    lambda1 = numeric(0), lambda2 = numeric(0),
    loglik_mean = numeric(0), loglik_se = numeric(0)
  )
  for (g in seq_len(n_priors)) visited[[paste0("lambda3_", g)]] <- numeric(0)

  best_loglik <- -Inf

  if (verbose) {
    cat(sprintf(
      "Search: %d lambdas x %d dims x %d folds\n",
      nlambda, n_dims, nfolds
    ))
  }

  # Coordinate descent: lambda3 -> lambda2 -> lambda1
  for (cycle in seq_len(3)) {
    prev_loglik <- best_loglik
    if (verbose) cat(sprintf("Cycle %d:\n", cycle))

    # lambda3[g]
    for (g in seq_len(n_priors)) {
      res <- search_dim(lambda_path, function(v) {
        trial <- best_lambda3
        trial[g] <- v
        c(best_lambda1, best_lambda2, trial)
      })
      best_lambda3[g] <- res$best
      cv_results[[2 + g]] <- data.frame(
        lambda = lambda_path, loglik = res$loglik, se = res$se
      )
      visited <- record_visited(
        visited, res,
        best_lambda1, best_lambda2, best_lambda3, 2 + g, lambda_path
      )
      if (verbose) {
        cat(sprintf(
          "  %s: %.5f (ll=%.4f)\n",
          dim_names[2 + g], best_lambda3[g], max(res$loglik)
        ))
      }
    }

    # lambda2
    res <- search_dim(lambda_path, function(v) {
      c(best_lambda1, v, best_lambda3)
    })
    best_lambda2 <- res$best
    cv_results[[2]] <- data.frame(
      lambda = lambda_path, loglik = res$loglik, se = res$se
    )
    visited <- record_visited(
      visited, res,
      best_lambda1, best_lambda2, best_lambda3, 2, lambda_path
    )
    if (verbose) {
      cat(sprintf(
        "  lambda2: %.5f (ll=%.4f)\n",
        best_lambda2, max(res$loglik)
      ))
    }

    # lambda1
    res <- search_dim(lambda_path, function(v) {
      c(v, best_lambda2, best_lambda3)
    })
    best_lambda1 <- res$best
    cv_results[[1]] <- data.frame(
      lambda = lambda_path, loglik = res$loglik, se = res$se
    )
    visited <- record_visited(
      visited, res,
      best_lambda1, best_lambda2, best_lambda3, 1, lambda_path
    )
    if (verbose) {
      cat(sprintf(
        "  lambda1: %.5f (ll=%.4f)\n",
        best_lambda1, max(res$loglik)
      ))
    }

    best_loglik <- max(cv_results[[1]]$loglik)
    if (verbose) cat(sprintf("  => ll=%.4f\n", best_loglik))

    if (is.finite(prev_loglik) &&
      abs(best_loglik - prev_loglik) / (abs(prev_loglik) + 1) < 1e-3) {
      if (verbose) cat("  Converged.\n")
      break
    }
  }

  # --- 1-SE rule + min target effective df (Hodges & Sargent, 2001) ---
  idx_max <- which.max(visited$loglik_mean)
  threshold <- visited$loglik_mean[idx_max] - visited$loglik_se[idx_max]
  candidates <- which(visited$loglik_mean >= threshold)

  calc_target_df <- function(fit) {
    coefs <- fit$coefficients
    group_sizes <- as.numeric(table(fit$group))
    group_sizes <- group_sizes[match(colnames(coefs), levels(fit$group))]
    n_target <- group_sizes[1]
    source_sizes <- group_sizes[-1]
    pm <- fit$prior_matrix
    df <- 0
    for (j in seq_len(nrow(coefs))) {
      if (abs(coefs[j, 1]) < 1e-8) next
      active_priors <- which(fit$active_prior[j, ])
      if (length(active_priors) > 0) {
        # Prior: df = n_target / (n_target + n_eff)
        weights <- pm[active_priors[1], ]
        nonzero_w <- which(weights > 0)
        n_effective <- 1 / sum(weights[nonzero_w]^2 / source_sizes[nonzero_w])
        df <- df + n_target / (n_target + n_effective)
      } else {
        # Local shared or group-specific: df = n_target / n_shared
        shared_groups <- abs(coefs[j, ] - coefs[j, 1]) < 1e-6
        df <- df + n_target / sum(group_sizes[shared_groups])
      }
    }
    df
  }

  candidate_df <- rep(Inf, length(candidates))
  for (ci in seq_along(candidates)) {
    idx <- candidates[ci]
    lambda3_cand <- sapply(seq_len(n_priors), function(g) {
      visited[[paste0("lambda3_", g)]][idx]
    })
    fit_cand <- tryCatch(
      coxtrans(formula, data, group, target,
        lambda1 = visited$lambda1[idx], lambda2 = visited$lambda2[idx],
        lambda3 = lambda3_cand,
        prior_matrix = prior_matrix, penalty = penalty, ...
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_cand)) candidate_df[ci] <- calc_target_df(fit_cand)
  }

  # Min target_df; ties broken by max CV loglik
  min_df <- min(candidate_df)
  ties <- which(abs(candidate_df - min_df) < 1e-8)
  best_ci <- ties[which.max(visited$loglik_mean[candidates[ties]])]
  best_idx <- candidates[best_ci]

  best_lambda1 <- visited$lambda1[best_idx]
  best_lambda2 <- visited$lambda2[best_idx]
  best_lambda3 <- sapply(seq_len(n_priors), function(g) {
    visited[[paste0("lambda3_", g)]][best_idx]
  })

  if (verbose) {
    cat(sprintf(
      "1-SE rule: %d candidates, min target_df=%.2f\n",
      length(candidates), min_df
    ))
    cat(sprintf(
      "Selected: lambda=(%.5f, %.5f, %s)\n",
      best_lambda1, best_lambda2,
      paste(round(best_lambda3, 5), collapse = ", ")
    ))
  }

  # Final refit on full data
  fit <- coxtrans(formula, data, group, target,
    lambda1 = best_lambda1, lambda2 = best_lambda2,
    lambda3 = best_lambda3,
    prior_matrix = prior_matrix, penalty = penalty, ...
  )

  list(
    best_lambdas = list(
      lambda1 = best_lambda1,
      lambda2 = best_lambda2,
      lambda3 = best_lambda3
    ),
    cv_results = cv_results,
    fit = fit
  )
}
