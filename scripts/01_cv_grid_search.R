#!/usr/bin/env Rscript
# 3D Grid Search CV for coxtrans
# Target group 做 K-fold CV，其他组作为辅助数据始终参与训练
# 评估指标：partial log-likelihood deviance（越小越好）

devtools::load_all()

formula <- Surv(time, status) ~ . - group - id
target <- 1
nfolds <- 5
seed <- 42

# Lambda grid
lambda_max <- calc_lambda_max(formula, sim2, sim2$group)
lambda1_seq <- exp(seq(log(lambda_max), log(lambda_max * 1e-4), length.out = 10))
lambda2_seq <- exp(seq(log(lambda_max), log(lambda_max * 1e-4), length.out = 5))
lambda3_seq <- exp(seq(log(lambda_max), log(lambda_max * 1e-4), length.out = 5))

grid <- expand.grid(
  lambda1 = lambda1_seq,
  lambda2 = lambda2_seq,
  lambda3 = lambda3_seq
)
cat(sprintf(
  "Grid: %d combinations x %d folds = %d fits\n",
  nrow(grid), nfolds, nrow(grid) * nfolds
))

# Target group 的数据
target_idx <- which(sim2$group == target)

# 创建 folds (仅在 target group 上)
set.seed(seed)
foldid <- sample(rep(seq_len(nfolds), length.out = length(target_idx)))

# CV 矩阵: nfolds x n_grid (存 deviance，越小越好)
cv_matrix <- matrix(NA_real_, nrow = nfolds, ncol = nrow(grid))
n_errors <- 0L

t_start <- proc.time()

for (k in seq_len(nfolds)) {
  cat(sprintf("Fold %d/%d ...\n", k, nfolds))
  test_fold <- which(foldid == k)
  train_rows <- setdiff(seq_len(nrow(sim2)), target_idx[test_fold])

  # 测试集数据（从 train_rows 构造的 model.matrix 提取对应行不可行，
  # 因此直接用全数据的 model.matrix，保证列一致）
  test_global_idx <- target_idx[test_fold]
  mf_test <- stats::model.frame(formula, sim2[test_global_idx, ])
  y_test <- stats::model.response(mf_test)
  x_test <- stats::model.matrix(formula, sim2[test_global_idx, ])[, -1]
  time_test <- y_test[, 1]
  status_test <- y_test[, 2]

  for (j in seq_len(nrow(grid))) {
    fit <- tryCatch(
      coxtrans(formula, sim2[train_rows, ], sim2$group[train_rows],
        target = target,
        lambda1 = grid$lambda1[j],
        lambda2 = grid$lambda2[j],
        lambda3 = grid$lambda3[j],
        penalty = "SCAD"
      ),
      error = function(e) {
        n_errors <<- n_errors + 1L
        NULL
      }
    )

    if (is.null(fit)) {
      cv_matrix[k, j] <- Inf
      next
    }

    coefs <- fit$coefficients
    beta <- coefs[, as.character(target)] + coefs[, "Center"]
    lp <- as.numeric(x_test %*% beta)

    # Partial log-likelihood deviance on test set
    # 按时间降序排列（与 coxtrans 内部一致）
    ord <- order(time_test, decreasing = TRUE)
    lp_ord <- lp[ord]
    status_ord <- status_test[ord]

    # -2 * partial log-likelihood (deviance)
    log_risk_cumsum <- log(cumsum(exp(lp_ord)))
    cv_matrix[k, j] <- -2 * sum(status_ord * (lp_ord - log_risk_cumsum))
  }
}

elapsed <- (proc.time() - t_start)[3]
cat(sprintf("\nTotal time: %.1f seconds (%d fitting errors)\n\n", elapsed, n_errors))

# 汇总 (deviance 越小越好)
cvm <- colMeans(cv_matrix)
cvsd <- apply(cv_matrix, 2, sd) / sqrt(nfolds)

# Best (最小 deviance)
idx_best <- which.min(cvm)
cat("=== Best lambda combination ===\n")
cat(sprintf("  lambda1 = %.5f\n", grid$lambda1[idx_best]))
cat(sprintf("  lambda2 = %.5f\n", grid$lambda2[idx_best]))
cat(sprintf("  lambda3 = %.5f\n", grid$lambda3[idx_best]))
cat(sprintf("  deviance = %.4f (sd = %.4f)\n\n", cvm[idx_best], cvsd[idx_best]))

# 用 best lambda 在全数据上 refit
fit_best <- coxtrans(formula, sim2, sim2$group,
  target = target,
  lambda1 = grid$lambda1[idx_best],
  lambda2 = grid$lambda2[idx_best],
  lambda3 = grid$lambda3[idx_best],
  penalty = "SCAD"
)

cat("=== Final model (best lambda, full data) ===\n")
cat("iter:", fit_best$iter, " ", fit_best$message, "\n\n")

beta_target <- fit_best$coefficients[, as.character(target)] +
  fit_best$coefficients[, "Center"]
cat("Target group beta (non-zero):\n")
print(round(beta_target[beta_target != 0], 4))

cat("\nTrue support: X1, X2, X3, X4\n")
cat("Selected support:", paste(names(beta_target[beta_target != 0]), collapse = ", "), "\n")
