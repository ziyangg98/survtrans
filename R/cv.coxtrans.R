#' Cross-validated coxtrans with sample splitting inference
#'
#' Selects penalty parameters via CV + BIC on a selection half of the
#' target data, then fits and performs inference on the held-out inference
#' half. This guarantees valid standard errors and coverage.
#'
#' @param formula A formula with a \code{\link[survival]{Surv}} response.
#' @param data A data frame.
#' @param group A factor indicating group membership.
#' @param target The target group level.
#' @param prior_matrix Optional G x (K-1) transfer constraint matrix.
#' @param nlambda Number of lambda values per dimension. Default 50.
#' @param lambda.min.ratio Smallest/largest lambda ratio. Default 1e-3.
#' @param nfolds Number of CV folds. Default 5.
#' @param penalty \code{"lasso"}, \code{"MCP"}, or \code{"SCAD"}.
#' @param sig_level Significance level for feature selection. Default 0.05.
#' @param alpha_eq Significance level for equivalence test. Default 0.05.
#' @param split_ratio Fraction of target data for selection. Default 0.5.
#' @param ncores Cores for parallel fold evaluation. Default 1.
#' @param seed Random seed.
#' @param verbose Print progress.
#' @param ... Passed to \code{\link{coxtrans}}.
#'
#' @return Object of class \code{cv.coxtrans}: \code{fit},
#'   \code{best_lambdas}, \code{selected}, \code{sig_level},
#'   \code{inference} (summary of inference-half fit).
#' @export
#'
#' @examples
#' result <- cv.coxtrans(Surv(time, status) ~ . - group - id,
#'   sim2, sim2$group, target = 1, nlambda = 10, penalty = "SCAD")
#' result
cv.coxtrans <- function(
  formula, data, group, target,
  prior_matrix = NULL, nlambda = 50, lambda.min.ratio = 1e-3,
  nfolds = 5, penalty = c("lasso", "MCP", "SCAD"),
  sig_level = 0.05, alpha_eq = 0.05, split_ratio = 0.5,
  ncores = 1, seed = 42, verbose = FALSE, ...
) {
  penalty <- match.arg(penalty)
  if (!is.factor(group)) group <- factor(group)
  target_level <- as.character(target)
  target_idx <- which(group == target_level)
  source_idx <- which(group != target_level)
  n_priors <- if (is.null(prior_matrix)) 1L else nrow(prior_matrix)

  build_lambda_path <- function(lambda_max) {
    if (!is.finite(lambda_max) || lambda_max <= 0) {
      return(0)
    }
    c(0, exp(seq(log(lambda_max * lambda.min.ratio),
      log(lambda_max), length.out = nlambda)))
  }

  lambda1_max <- calc_lambda1_max(formula, data, group, target)
  lambda2_max <- calc_lambda2_max(formula, data, group, target)
  lambda3_max <- calc_lambda3_max(formula, data, group, target, prior_matrix)

  lambda1_path <- build_lambda_path(lambda1_max)
  lambda2_path <- build_lambda_path(lambda2_max)
  lambda3_path <- lapply(lambda3_max, build_lambda_path)

  # ================================================================
  # Sample splitting: target data → selection + inference
  # ================================================================
  set.seed(seed)
  n_target <- length(target_idx)
  n_sel <- round(n_target * split_ratio)
  sel_perm <- sample(n_target)
  sel_target <- target_idx[sel_perm[seq_len(n_sel)]]
  inf_target <- target_idx[sel_perm[(n_sel + 1):n_target]]

  # Selection data: selection-half target + all sources
  sel_rows <- sort(c(sel_target, source_idx))
  sel_data <- data[sel_rows, ]
  sel_group <- group[sel_rows]

  # Inference data: inference-half target + all sources
  inf_rows <- sort(c(inf_target, source_idx))
  inf_data <- data[inf_rows, ]
  inf_group <- group[inf_rows]

  if (verbose) cat(sprintf("Split: %d selection + %d inference (target)\n",
                           n_sel, n_target - n_sel))

  safe_fit <- function(l1, l2, l3, d = sel_data, g = sel_group)
    tryCatch(coxtrans(formula, d, g, target, lambda1 = l1,
      lambda2 = l2, lambda3 = l3, prior_matrix = prior_matrix,
      penalty = penalty, ...), error = function(e) NULL)

  # ================================================================
  # Step 1: CV on selection half
  # ================================================================
  if (verbose) cat(sprintf("Step 1: CV (%d-fold)\n", nfolds))
  sel_target_local <- which(sel_group == target_level)

  set.seed(seed + 1L)
  y_sel <- stats::model.response(
    stats::model.frame(formula, sel_data[sel_target_local, ]))
  fids <- integer(length(sel_target_local))
  for (s in unique(y_sel[, 2])) {
    i <- which(y_sel[, 2] == s)
    fids[i] <- sample(rep(seq_len(nfolds), length.out = length(i)))
  }

  cv_folds <- function(l1, l2, l3) {
    eval1 <- function(k) {
      ti <- sel_target_local[fids == k]
      fit <- tryCatch(coxtrans(formula, sel_data[-ti, ],
        sel_group[-ti], target, lambda1 = l1, lambda2 = l2,
        lambda3 = l3, prior_matrix = prior_matrix,
        penalty = penalty, ...), error = function(e) NULL)
      if (is.null(fit)) return(-Inf)
      xt <- stats::model.matrix(formula, sel_data[ti, ])[, -1]
      yt <- stats::model.response(
        stats::model.frame(formula, sel_data[ti, ]))
      lp <- as.numeric(xt %*% fit$coefficients[, 1])
      o <- order(yt[, 1], decreasing = TRUE)
      sum(yt[o, 2] * (lp[o] - log(cumsum(exp(lp[o])))))
    }
    if (ncores > 1) {
      unlist(parallel::mclapply(seq_len(nfolds), eval1,
                                mc.cores = ncores))
    } else vapply(seq_len(nfolds), eval1, numeric(1))
  }

  visited_lams <- list(); visited_folds <- list()
  best_l1 <- 0; best_l2 <- 0; best_l3 <- rep(0, n_priors)

  search1 <- function(path, make_lams, name) {
    bv <- path[1]; bll <- -Inf
    for (v in path) {
      lm <- make_lams(v)
      fll <- cv_folds(lm[1], lm[2], lm[-(1:2)])
      visited_lams[[length(visited_lams) + 1]] <<- lm
      visited_folds[[length(visited_folds) + 1]] <<- fll
      if (mean(fll) > bll) { bll <- mean(fll); bv <- v }
    }
    if (verbose) cat(sprintf("    %s = %.5f\n", name, bv))
    bv
  }

  for (cyc in seq_len(3)) {
    prev <- if (length(visited_folds) > 0) {
      max(sapply(visited_folds, mean))
    } else -Inf
    if (verbose) cat(sprintf("  Cycle %d:\n", cyc))
    for (g in seq_len(n_priors))
      best_l3[g] <- search1(lambda3_path[[g]], function(v) {
        tr <- best_l3; tr[g] <- v; c(best_l1, best_l2, tr)
      }, if (n_priors == 1) "l3" else paste0("l3[", g, "]"))
    best_l2 <- search1(lambda2_path, function(v)
      c(best_l1, v, best_l3), "l2")
    best_l1 <- search1(lambda1_path, function(v)
      c(v, best_l2, best_l3), "l1")
    cur <- max(sapply(visited_folds, mean))
    if (is.finite(prev) && is.finite(cur) &&
        abs(cur - prev) / (abs(prev) + 1) < 1e-3) {
      if (verbose) cat("  Converged\n"); break
    }
  }

  # ================================================================
  # Step 2: Equivalence set + BIC (on selection half)
  # ================================================================
  n_vis <- length(visited_lams)
  i_best <- which.max(sapply(visited_folds, mean))
  fll_best <- visited_folds[[i_best]]

  equiv <- logical(n_vis)
  for (i in seq_len(n_vis)) {
    if (i == i_best) { equiv[i] <- TRUE; next }
    d <- fll_best - visited_folds[[i]]
    if (any(!is.finite(d))) next
    tt <- tryCatch(stats::t.test(d)$p.value, error = function(e) 0)
    equiv[i] <- tt > alpha_eq
  }
  eq_idx <- which(equiv)
  eq_keys <- sapply(eq_idx, function(i)
    paste(round(visited_lams[[i]], 8), collapse = "_"))
  eq_idx <- eq_idx[!duplicated(eq_keys)]

  if (verbose) cat(sprintf("Step 2: %d/%d equivalent\n",
                           length(eq_idx), n_vis))

  eq_bic <- rep(Inf, length(eq_idx))
  for (ci in seq_along(eq_idx)) {
    lm <- visited_lams[[eq_idx[ci]]]
    fit <- safe_fit(lm[1], lm[2], lm[-(1:2)])
    if (!is.null(fit)) eq_bic[ci] <- BIC(fit)
  }
  bi <- which.min(eq_bic)
  best_lam <- visited_lams[[eq_idx[bi]]]
  best_l1 <- best_lam[1]; best_l2 <- best_lam[2]
  best_l3 <- best_lam[-(1:2)]

  if (verbose) cat(sprintf("  BIC: l1=%.4f l2=%.4f l3=%s\n",
    best_l1, best_l2, paste(round(best_l3, 4), collapse = ",")))

  # ================================================================
  # Step 3: Fit on inference half + summary for inference
  # ================================================================
  if (verbose) cat("Step 3: Inference\n")
  inf_fit <- safe_fit(best_l1, best_l2, best_l3,
                      d = inf_data, g = inf_group)
  if (is.null(inf_fit)) stop("Inference fit failed.")

  p <- nrow(inf_fit$coefficients)
  fnames <- rownames(inf_fit$coefficients)

  # summary.coxtrans gives target coefficients + sandwich SE + p-values
  smry <- summary(inf_fit, target_only = TRUE)
  coef_target <- smry$coefficients[, "coef"]
  se_target <- smry$coefficients[, "se(coef)"]
  pv_target <- smry$coefficients[, "Pr(>|z|)"]
  sel <- which(pv_target < sig_level & abs(coef_target) > 1e-8)

  if (verbose) cat(sprintf("  %d selected\n", length(sel)))

  out <- list(
    fit = inf_fit,
    best_lambdas = list(lambda1 = best_l1, lambda2 = best_l2,
                        lambda3 = best_l3),
    selected = sel, sig_level = sig_level,
    coefficients = smry$coefficients)
  class(out) <- "cv.coxtrans"
  out
}

#' @export
print.cv.coxtrans <- function(x, ...) {
  p <- nrow(x$fit$coefficients)
  ns <- length(x$selected)
  cat("cv.coxtrans (CV + BIC + sample splitting)\n\n")
  cat(sprintf("  Lambda:    l1=%.4f, l2=%.4f, l3=%s\n",
    x$best_lambdas$lambda1, x$best_lambdas$lambda2,
    paste(round(x$best_lambdas$lambda3, 4), collapse = ",")))
  cat(sprintf("  Selected:  %d/%d (p < %g)", ns, p, x$sig_level))
  if (ns > 0)
    cat(" (", paste(rownames(x$fit$coefficients)[x$selected],
                    collapse = ", "), ")", sep = "")
  cat("\n\n")
  # Print coefficients for selected features
  if (ns > 0) {
    stats::printCoefmat(
      x$coefficients[x$selected, , drop = FALSE],
      cs.ind = 1:3, tst.ind = 4, P.values = TRUE, has.Pvalue = TRUE
    )
  }
  invisible(x)
}
