#!/usr/bin/env Rscript
# 对比 pooled coxph vs coxtrans (lambda2=lambda3=很大)
# 验证 coxtrans 在强 transfer shrinkage 下退化为 pooled 模型

devtools::load_all()

formula <- Surv(time, status) ~ . - group - id
target <- 1

# ============================================================
# 1. Pooled: 忽略 group 差异，所有数据合在一起跑 SCAD Cox
# ============================================================
library(ncvreg)
x_all <- as.matrix(sim2[, paste0("X", 1:20)])
y_all <- survival::Surv(sim2$time, sim2$status)

set.seed(42)
cv_pool <- cv.ncvsurv(x_all, y_all, penalty = "SCAD", nfolds = 5)
lambda_pool <- cv_pool$lambda.min

beta_pool <- as.numeric(coef(cv_pool, lambda = lambda_pool))
names(beta_pool) <- paste0("X", 1:20)

cat("=== Pooled SCAD (all groups, CV) ===\n")
cat(sprintf("  lambda = %.5f\n", lambda_pool))
cat("  Support:", paste(names(beta_pool[beta_pool != 0]), collapse = ", "), "\n")
print(round(beta_pool[beta_pool != 0], 4))

# ============================================================
# 2. coxtrans: lambda2=lambda3=很大，强制各组系数趋同
# ============================================================
large_lambda <- 1e2
lambda_max <- calc_lambda_max(formula, sim2, sim2$group)
lambda1_seq <- exp(seq(log(lambda_max), log(lambda_max * 1e-4), length.out = 30))

# 用 CV 选 lambda1，固定 lambda2=lambda3=large_lambda
target_idx <- which(sim2$group == target)
nfolds <- 5
set.seed(42)
foldid <- sample(rep(seq_len(nfolds), length.out = length(target_idx)))

cv_coxtrans <- matrix(NA_real_, nrow = nfolds, ncol = length(lambda1_seq))

for (k in seq_len(nfolds)) {
  cat(sprintf("Fold %d/%d ...\n", k, nfolds))
  test_fold <- which(foldid == k)
  train_rows <- setdiff(seq_len(nrow(sim2)), target_idx[test_fold])

  test_global_idx <- target_idx[test_fold]
  mf_test <- stats::model.frame(formula, sim2[test_global_idx, ])
  y_test <- stats::model.response(mf_test)
  x_test <- stats::model.matrix(formula, sim2[test_global_idx, ])[, -1]
  time_test <- y_test[, 1]
  status_test <- y_test[, 2]

  for (j in seq_along(lambda1_seq)) {
    fit <- tryCatch(
      coxtrans(formula, sim2[train_rows, ], sim2$group[train_rows],
        target = target,
        lambda1 = lambda1_seq[j], lambda2 = large_lambda, lambda3 = large_lambda,
        penalty = "SCAD"
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      cv_coxtrans[k, j] <- Inf
      next
    }

    coefs <- fit$coefficients
    beta <- coefs[, as.character(target)] + coefs[, "Center"]
    lp <- as.numeric(x_test %*% beta)
    ord <- order(time_test, decreasing = TRUE)
    lp_ord <- lp[ord]
    status_ord <- status_test[ord]
    cv_coxtrans[k, j] <- -2 * sum(status_ord * (lp_ord - log(cumsum(exp(lp_ord)))))
  }
}

idx_best <- which.min(colMeans(cv_coxtrans))
lambda1_best <- lambda1_seq[idx_best]

fit_ct <- coxtrans(formula, sim2, sim2$group,
  target = target,
  lambda1 = lambda1_best, lambda2 = large_lambda, lambda3 = large_lambda,
  penalty = "SCAD"
)

beta_coxtrans <- fit_ct$coefficients[, as.character(target)] +
  fit_ct$coefficients[, "Center"]

cat("\n=== coxtrans SCAD (lambda2=lambda3=", large_lambda, ", CV) ===\n")
cat(sprintf("  lambda1 = %.5f\n", lambda1_best))
cat("  Support:", paste(names(beta_coxtrans[beta_coxtrans != 0]), collapse = ", "), "\n")
print(round(beta_coxtrans[beta_coxtrans != 0], 4))

# 检查各组系数是否趋同
cat("\n=== Group coefficients (should be similar across groups) ===\n")
coefs <- fit_ct$coefficients
for (col in colnames(coefs)) {
  nonzero <- coefs[coefs[, col] != 0, col]
  if (length(nonzero) > 0) {
    cat(sprintf("  %s: ", col))
    print(round(nonzero, 4))
  }
}

# ============================================================
# 3. 对比表
# ============================================================
cat("\n=== Comparison ===\n")
cat(sprintf(
  "%-8s  %6s  %10s  %10s\n",
  "Feature", "truth", "pooled", "coxtrans"
))
cat(strrep("-", 40), "\n")

true_beta <- setNames(rep(0, 20), paste0("X", 1:20))
true_beta[1:4] <- 0.3

all_nonzero <- sort(unique(c(
  names(beta_pool[beta_pool != 0]),
  names(beta_coxtrans[beta_coxtrans != 0])
)))
for (feat in all_nonzero) {
  cat(sprintf(
    "%-8s  %6.1f  %10.4f  %10.4f\n",
    feat, true_beta[feat],
    beta_pool[feat], beta_coxtrans[feat]
  ))
}

cat("\nTrue beta (target): X1=0.3, X2=0.3, X3=0.3, X4=0.3\n")
cat("Note: pooled model ignores group heterogeneity,\n")
cat("  so X1/X2 estimates are biased by source groups' different effects.\n")
