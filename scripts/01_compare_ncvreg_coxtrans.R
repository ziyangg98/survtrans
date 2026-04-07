#!/usr/bin/env Rscript
# Compare ncvreg (target only) vs coxtrans (with transfer)
# on sim2 data to demonstrate the benefit of transfer learning.

devtools::load_all()
library(ncvreg)

formula <- Surv(time, status) ~ . - group - id
target <- 1
target_data <- sim2[sim2$group == target, ]

# === 1. ncvreg: SCAD on target only ===
x_target <- as.matrix(target_data[, paste0("X", 1:20)])
y_target <- survival::Surv(target_data$time, target_data$status)

set.seed(42)
cv_ncvreg <- cv.ncvsurv(x_target, y_target, penalty = "SCAD", nfolds = 10)
beta_ncvreg <- as.numeric(coef(cv_ncvreg, lambda = cv_ncvreg$lambda.min))
names(beta_ncvreg) <- paste0("X", 1:20)

cat("=== ncvreg SCAD (target only) ===\n")
cat("  Support:", paste(names(beta_ncvreg[beta_ncvreg != 0]), collapse = ", "), "\n")
print(round(beta_ncvreg[beta_ncvreg != 0], 4))

# === 2. coxtrans: with transfer (auto-tuned) ===
cv_fit <- cv.coxtrans(formula, sim2, sim2$group, target = target,
  nlambda = 20, penalty = "SCAD", verbose = TRUE)

cat("\n=== coxtrans SCAD (with transfer, cv.coxtrans) ===\n")
summary(cv_fit$fit)

# === 3. Comparison ===
beta_coxtrans <- cv_fit$fit$coefficients[, as.character(target)]

cat("\n=== Comparison ===\n")
cat(sprintf("%-8s  %6s  %10s  %10s\n", "Feature", "truth", "ncvreg", "coxtrans"))
cat(strrep("-", 40), "\n")

true_beta <- setNames(rep(0, 20), paste0("X", 1:20))
true_beta[1:4] <- 0.3

all_nonzero <- sort(unique(c(
  names(beta_ncvreg[beta_ncvreg != 0]),
  names(beta_coxtrans[abs(beta_coxtrans) > 1e-8])
)))
for (feat in all_nonzero) {
  cat(sprintf("%-8s  %6.1f  %10.4f  %10.4f\n",
    feat, true_beta[feat], beta_ncvreg[feat], beta_coxtrans[feat]))
}
