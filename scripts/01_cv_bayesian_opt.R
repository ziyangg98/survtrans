#!/usr/bin/env Rscript
# Bayesian Optimization CV for coxtrans
# 依赖: ParBayesianOptimization

if (!requireNamespace("ParBayesianOptimization", quietly = TRUE)) {
  stop("请先安装: renv::install('AnotherSamWilson/ParBayesianOptimization')")
}

library(ParBayesianOptimization)
devtools::load_all()

formula <- Surv(time, status) ~ . - group - id
target <- 1
nfolds <- 5
seed <- 42

# Lambda 搜索范围 (log scale, 2 个数量级)
lambda_max <- calc_lambda_max(formula, sim2, sim2$group)
log_lam_lo <- log(lambda_max * 0.01)
log_lam_hi <- log(lambda_max)

bounds <- list(
  log_lambda1 = c(log_lam_lo, log_lam_hi),
  log_lambda2 = c(log_lam_lo, log_lam_hi),
  log_lambda3 = c(log_lam_lo, log_lam_hi)
)

# 预计算 CV folds
target_idx <- which(sim2$group == target)
set.seed(seed)
foldid <- sample(rep(seq_len(nfolds), length.out = length(target_idx)))

fold_data <- lapply(seq_len(nfolds), function(k) {
  test_global_idx <- target_idx[foldid == k]
  mf <- stats::model.frame(formula, sim2[test_global_idx, ])
  y <- stats::model.response(mf)
  list(
    train_rows = setdiff(seq_len(nrow(sim2)), test_global_idx),
    x_test = stats::model.matrix(formula, sim2[test_global_idx, ])[, -1],
    time_test = y[, 1], status_test = y[, 2]
  )
})

# CV 目标函数
cv_objective <- function(log_lambda1, log_lambda2, log_lambda3) {
  l1 <- exp(log_lambda1)
  l2 <- exp(log_lambda2)
  l3 <- exp(log_lambda3)
  deviances <- numeric(nfolds)
  for (k in seq_len(nfolds)) {
    fd <- fold_data[[k]]
    fit <- tryCatch(
      coxtrans(formula, sim2[fd$train_rows, ], sim2$group[fd$train_rows],
        target = target, lambda1 = l1, lambda2 = l2, lambda3 = l3,
        penalty = "SCAD"
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(Score = -Inf))
    }
    beta <- fit$coefficients[, as.character(target)] + fit$coefficients[, "Center"]
    lp <- as.numeric(fd$x_test %*% beta)
    ord <- order(fd$time_test, decreasing = TRUE)
    deviances[k] <- -2 * sum(fd$status_test[ord] *
      (lp[ord] - log(cumsum(exp(lp[ord])))))
  }
  list(Score = -mean(deviances))
}

t_start <- proc.time()
set.seed(seed)
opt_result <- bayesOpt(
  FUN = cv_objective, bounds = bounds,
  initPoints = 25, iters.n = 15, iters.k = 1,
  acq = "ucb", kappa = 2.576, gsPoints = 30,
  acqThresh = 0.8, errorHandling = 5L, verbose = 1
)
elapsed <- (proc.time() - t_start)[3]

# 结果
best <- getBestPars(opt_result)
best_l <- exp(unlist(best))
best_dev <- -max(opt_result$scoreSummary$Score)
n_eval <- nrow(opt_result$scoreSummary)

cat(sprintf("\nTime: %.1fs | %d evals | deviance: %.4f\n", elapsed, n_eval, best_dev))
cat(sprintf(
  "lambda1=%.5f  lambda2=%.5f  lambda3=%.5f\n\n",
  best_l[1], best_l[2], best_l[3]
))

# Refit on full data
fit_best <- coxtrans(formula, sim2, sim2$group,
  target = target, lambda1 = best_l[1], lambda2 = best_l[2],
  lambda3 = best_l[3], penalty = "SCAD"
)
beta_target <- fit_best$coefficients[, as.character(target)] +
  fit_best$coefficients[, "Center"]
cat("Selected:", paste(names(beta_target[beta_target != 0]), collapse = ", "), "\n")
cat("True:     X1, X2, X3, X4\n")
