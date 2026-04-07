#!/usr/bin/env Rscript
# Compare pooled Cox (ignores group structure) vs coxtrans (with transfer)
# on sim2 data. Shows that coxtrans correctly handles heterogeneous effects.

devtools::load_all()
library(ncvreg)

formula <- Surv(time, status) ~ . - group - id
target <- 1

# === 1. Pooled: ignore groups, fit SCAD Cox on all data ===
x_all <- as.matrix(sim2[, paste0("X", 1:20)])
y_all <- survival::Surv(sim2$time, sim2$status)

set.seed(42)
cv_pool <- cv.ncvsurv(x_all, y_all, penalty = "SCAD", nfolds = 10)
beta_pool <- as.numeric(coef(cv_pool, lambda = cv_pool$lambda.min))
names(beta_pool) <- paste0("X", 1:20)

cat("=== Pooled SCAD (all groups, ignoring heterogeneity) ===\n")
cat("  Support:", paste(names(beta_pool[beta_pool != 0]), collapse = ", "), "\n")
print(round(beta_pool[beta_pool != 0], 4))

# === 2. coxtrans: with transfer (auto-tuned) ===
cv_fit <- cv.coxtrans(formula, sim2, sim2$group, target = target,
  nlambda = 20, penalty = "SCAD", verbose = TRUE)

cat("\n=== coxtrans SCAD (with transfer, cv.coxtrans) ===\n")
summary(cv_fit$fit)

# === 3. Comparison ===
beta_coxtrans <- cv_fit$fit$coefficients[, as.character(target)]

cat("\n=== Comparison ===\n")
cat(sprintf("%-8s  %6s  %10s  %10s\n", "Feature", "truth", "pooled", "coxtrans"))
cat(strrep("-", 40), "\n")

true_beta <- setNames(rep(0, 20), paste0("X", 1:20))
true_beta[1:4] <- 0.3

all_nonzero <- sort(unique(c(
  names(beta_pool[beta_pool != 0]),
  names(beta_coxtrans[abs(beta_coxtrans) > 1e-8])
)))
for (feat in all_nonzero) {
  cat(sprintf("%-8s  %6.1f  %10.4f  %10.4f\n",
    feat, true_beta[feat], beta_pool[feat], beta_coxtrans[feat]))
}

cat("\nNote: pooled model ignores group heterogeneity,\n")
cat("  so X1/X2 estimates may be biased by source groups' different effects.\n")
