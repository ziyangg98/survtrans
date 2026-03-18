#!/usr/bin/env Rscript
# 对比 ncvreg (target only) vs coxtrans (lambda2=lambda3=0)
# 各自 CV 选 lambda，比较变量选择和系数估计

devtools::load_all()
library(ncvreg)

formula <- Surv(time, status) ~ . - group - id
target <- 1
target_data <- sim2[sim2$group == target, ]
x_target <- as.matrix(target_data[, paste0("X", 1:20)])
y_target <- survival::Surv(target_data$time, target_data$status)
target_idx <- which(sim2$group == target)

nfolds <- 5
set.seed(42)
foldid <- sample(rep(seq_len(nfolds), length.out = length(target_idx)))

# ============================================================
# 1. ncvreg: SCAD penalized Cox on target only (自带 CV)
# ============================================================
set.seed(42)
cv_ncvreg <- cv.ncvsurv(x_target, y_target, penalty = "SCAD", nfolds = nfolds)
lambda_ncvreg <- cv_ncvreg$lambda.min

beta_ncvreg <- as.numeric(coef(cv_ncvreg, lambda = lambda_ncvreg))
names(beta_ncvreg) <- paste0("X", 1:20)

cat("=== ncvreg SCAD (target only, CV) ===\n")
cat(sprintf("  lambda = %.5f\n", lambda_ncvreg))
cat("  Support:", paste(names(beta_ncvreg[beta_ncvreg != 0]), collapse = ", "), "\n")
print(round(beta_ncvreg[beta_ncvreg != 0], 4))

# ============================================================
# 2. coxtrans: lambda2=lambda3=0, CV 选 lambda1
# ============================================================
lambda_max_ct <- calc_lambda_max(formula, sim2, sim2$group)
lambda1_seq <- exp(seq(log(lambda_max_ct), log(lambda_max_ct * 1e-4), length.out = 30))

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
        lambda1 = lambda1_seq[j], lambda2 = 0, lambda3 = 0,
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

idx_best_ct <- which.min(colMeans(cv_coxtrans))
lambda_ct <- lambda1_seq[idx_best_ct]
fit_ct <- coxtrans(formula, sim2, sim2$group,
  target = target,
  lambda1 = lambda_ct, lambda2 = 0, lambda3 = 0, penalty = "SCAD"
)
beta_coxtrans <- fit_ct$coefficients[, as.character(target)] +
  fit_ct$coefficients[, "Center"]

cat("\n=== coxtrans SCAD (lambda2=lambda3=0, CV) ===\n")
cat(sprintf("  lambda1 = %.5f\n", lambda_ct))
cat("  Support:", paste(names(beta_coxtrans[beta_coxtrans != 0]), collapse = ", "), "\n")
print(round(beta_coxtrans[beta_coxtrans != 0], 4))

# ============================================================
# 3. 对比表
# ============================================================
cat("\n=== Comparison ===\n")
cat(sprintf(
  "%-8s  %6s  %10s  %10s\n",
  "Feature", "truth", "ncvreg", "coxtrans"
))
cat(strrep("-", 40), "\n")

true_beta <- setNames(rep(0, 20), paste0("X", 1:20))
true_beta[1:4] <- 0.3

all_nonzero <- sort(unique(c(
  names(beta_ncvreg[beta_ncvreg != 0]),
  names(beta_coxtrans[beta_coxtrans != 0])
)))
for (feat in all_nonzero) {
  cat(sprintf(
    "%-8s  %6.1f  %10.4f  %10.4f\n",
    feat, true_beta[feat],
    beta_ncvreg[feat], beta_coxtrans[feat]
  ))
}

cat("\nTrue beta (target): X1=0.3, X2=0.3, X3=0.3, X4=0.3\n")
