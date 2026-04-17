#!/usr/bin/env Rscript
# Compare BIC.coxtrans vs cv.coxtrans vs ncvreg on sim2

devtools::load_all()
library(ncvreg)

formula <- Surv(time, status) ~ . - group - id
target <- 1
penalty <- "SCAD"
nlambda <- 20
p <- 20

true_beta <- setNames(rep(0, p), paste0("X", 1:p))
true_beta[1:4] <- 0.3

# === 1. BIC selection (no CV) ===
lambda_max <- calc_lambda_max(formula, sim2, sim2$group)
lambda_path <- exp(seq(log(lambda_max * 1e-3), log(lambda_max),
  length.out = nlambda))

best_l1 <- 0; best_l2 <- 0; best_l3 <- 0
cat("=== BIC.coxtrans ===\n")
for (cycle in 1:3) {
  bic_vals <- numeric(nlambda)
  # lambda3
  for (i in seq_along(lambda_path)) {
    fit <- tryCatch(coxtrans(formula, sim2, sim2$group, target,
      lambda1 = best_l1, lambda2 = best_l2, lambda3 = lambda_path[i],
      penalty = penalty), error = function(e) NULL)
    bic_vals[i] <- if (is.null(fit)) Inf else BIC(fit)
  }
  best_l3 <- lambda_path[which.min(bic_vals)]
  # lambda2
  for (i in seq_along(lambda_path)) {
    fit <- tryCatch(coxtrans(formula, sim2, sim2$group, target,
      lambda1 = best_l1, lambda2 = lambda_path[i], lambda3 = best_l3,
      penalty = penalty), error = function(e) NULL)
    bic_vals[i] <- if (is.null(fit)) Inf else BIC(fit)
  }
  best_l2 <- lambda_path[which.min(bic_vals)]
  # lambda1
  for (i in seq_along(lambda_path)) {
    fit <- tryCatch(coxtrans(formula, sim2, sim2$group, target,
      lambda1 = lambda_path[i], lambda2 = best_l2, lambda3 = best_l3,
      penalty = penalty), error = function(e) NULL)
    bic_vals[i] <- if (is.null(fit)) Inf else BIC(fit)
  }
  best_l1 <- lambda_path[which.min(bic_vals)]
  cat(sprintf("  Cycle %d: lambda=(%.4f, %.4f, %.4f) BIC=%.2f\n",
    cycle, best_l1, best_l2, best_l3, min(bic_vals)))
}
fit_bic <- coxtrans(formula, sim2, sim2$group, target,
  lambda1 = best_l1, lambda2 = best_l2, lambda3 = best_l3, penalty = penalty)
summary(fit_bic)

# === 2. CV 1-SE selection ===
cat("\n=== cv.coxtrans ===\n")
cv_fit <- cv.coxtrans(formula, sim2, sim2$group, target = target,
  nlambda = nlambda, penalty = penalty, verbose = TRUE)
summary(cv_fit$fit)

# === 3. ncvreg (target only) ===
cat("\n=== ncvreg (target only) ===\n")
target_data <- sim2[sim2$group == target, ]
x_target <- as.matrix(target_data[, paste0("X", 1:p)])
y_target <- Surv(target_data$time, target_data$status)
set.seed(42)
cv_ncvreg <- cv.ncvsurv(x_target, y_target, penalty = "SCAD", nfolds = 10)
beta_ncvreg <- as.numeric(coef(cv_ncvreg, lambda = cv_ncvreg$lambda.min))
names(beta_ncvreg) <- paste0("X", 1:p)
cat("  Support:", paste(names(beta_ncvreg[beta_ncvreg != 0]), collapse = ", "), "\n")

# === Summary ===
beta_bic <- fit_bic$coefficients[, 1]
beta_cv <- cv_fit$fit$coefficients[, 1]
true_nz <- names(true_beta[true_beta != 0])

eval_method <- function(name, beta) {
  sel <- names(beta[abs(beta) > 1e-8])
  tp <- length(intersect(sel, true_nz))
  fp <- length(setdiff(sel, true_nz))
  fn <- length(setdiff(true_nz, sel))
  mse <- mean((beta - true_beta)^2)
  cat(sprintf("  %-12s MSE=%.6f  TP=%d  FP=%d  FN=%d\n", name, mse, tp, fp, fn))
}

cat("\n=== Comparison ===\n")
eval_method("BIC", beta_bic)
eval_method("CV-1SE", beta_cv)
eval_method("ncvreg", beta_ncvreg)

cat(sprintf("\n%-6s  %5s  %8s  %8s  %8s\n", "Feat", "truth", "BIC", "CV-1SE", "ncvreg"))
cat(strrep("-", 45), "\n")
all_nz <- sort(unique(c(true_nz,
  names(beta_bic[abs(beta_bic) > 1e-8]),
  names(beta_cv[abs(beta_cv) > 1e-8]),
  names(beta_ncvreg[abs(beta_ncvreg) > 1e-8]))))
for (f in all_nz) {
  cat(sprintf("%-6s  %5.1f  %8.4f  %8.4f  %8.4f\n",
    f, true_beta[f], beta_bic[f], beta_cv[f], beta_ncvreg[f]))
}
