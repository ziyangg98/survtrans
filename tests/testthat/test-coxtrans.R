test_that("coxtrans basic fit returns correct structure", {
  formula <- Surv(time, status) ~ . - group - id
  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
  )

  expect_s3_class(fit, "coxtrans")
  expect_true(is.matrix(fit$coefficients))
  expect_equal(nrow(fit$coefficients), 20L)
  expect_equal(colnames(fit$coefficients)[1], "1")
  expect_true(!is.null(fit$prior_matrix))
})

test_that("coxtrans rejects invalid target", {
  formula <- Surv(time, status) ~ . - group - id
  expect_error(
    coxtrans(formula, sim2, sim2$group, "nonexistent",
      lambda1 = 0.1, lambda2 = 0.1, lambda3 = 0.1
    ),
    "target.*not found in group levels"
  )
})

test_that("coxtrans S3 methods work", {
  formula <- Surv(time, status) ~ . - group - id
  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
  )

  # coef
  cf <- coef(fit)
  expect_true(is.numeric(cf))
  expect_true(length(cf) > 0)

  # predict (needs newdata with dot formula)
  pred <- predict(fit, newdata = sim2)
  expect_true(is.numeric(pred))

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.coxtrans")

  # basehaz
  bh <- basehaz(fit)
  expect_true(is.list(bh))

  # logLik
  ll <- logLik(fit)
  expect_true(is.numeric(ll))
})

test_that("coxtrans with custom prior_matrix works", {
  formula <- Surv(time, status) ~ . - group - id
  pm <- rbind(
    tissue_A = c(0.5, 0.5, 0, 0),
    tissue_B = c(0, 0, 0.5, 0.5)
  )
  colnames(pm) <- c("2", "3", "4", "5")
  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 0.075, lambda2 = 0.04, lambda3 = c(0.04, 0.04),
    prior_matrix = pm, penalty = "SCAD"
  )

  expect_s3_class(fit, "coxtrans")
  expect_equal(nrow(fit$prior_matrix), 2L)
  expect_equal(ncol(fit$active_prior), 2L)
  expect_equal(nrow(fit$active_local), 20L)
  expect_true(is.matrix(fit$prior_effects))
})

test_that("coxtrans feature structure is correct on sim2", {
  formula <- Surv(time, status) ~ . - group - id
  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
  )

  # X1-X4 should be non-zero
  target_coefs <- fit$coefficients[, 1]
  expect_true(all(abs(target_coefs[1:4]) > 0.1))

  # X5-X20 should be zero (sparse)
  expect_true(all(target_coefs[5:20] == 0))

  # X3, X4 should be shared (local active for all sources)
  expect_true(all(fit$active_local[3, ]))
  expect_true(all(fit$active_local[4, ]))

  # X1, X2 should not be shared (local inactive)
  expect_true(!all(fit$active_local[1, ]))
  expect_true(!all(fit$active_local[2, ]))
})

test_that("coxtrans return object has expected fields", {
  formula <- Surv(time, status) ~ . - group - id
  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
  )

  expected_fields <- c(
    "coefficients", "prior_matrix", "prior_effects",
    "active_local", "active_prior", "iter", "message",
    "history", "penalty", "lambda1", "lambda2", "lambda3",
    "formula", "call", "time", "status", "group", "x"
  )
  expect_true(all(expected_fields %in% names(fit)))
})

test_that("print.coxtrans produces output", {
  formula <- Surv(time, status) ~ . - group - id
  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
  )

  output <- capture.output(print(fit))
  expect_true(any(grepl("Feature structure", output)))
  expect_true(any(grepl("Prior transfer|Shared|Sparse", output)))
})

test_that("ncvcox basic fit returns correct structure", {
  formula <- Surv(time, status) ~ . - group - id
  df <- sim2[sim2$group == 2 | sim2$group == 4, ]
  fit <- ncvcox(formula, df, df$group, lambda = 0.1, penalty = "SCAD")

  expect_s3_class(fit, "ncvcox")
  expect_true(is.numeric(fit$coefficients))
  expect_equal(length(fit$coefficients), 20L)
})

test_that("cv.coxtrans returns a cv.coxtrans object with glmnet-like fields", {
  formula <- Surv(time, status) ~ . - group - id
  result <- cv.coxtrans(
    formula, sim2, sim2$group,
    target = 1,
    nlambda = 5, penalty = "SCAD"
  )

  expect_s3_class(result, "cv.coxtrans")
  expect_s3_class(result$coxtrans.fit, "coxtrans")

  # CV diagnostics
  n_grid <- nrow(result$grid)
  expect_true(n_grid > 0)
  expect_equal(length(result$cvm), n_grid)
  expect_equal(length(result$cvsd), n_grid)
  expect_equal(length(result$nzero), n_grid)
  expect_true(any(is.finite(result$cvm)))

  # Optimal lambdas
  expect_true(
    all(c("lambda1", "lambda2", "lambda3") %in% names(result$lambda.min))
  )
  expect_true(is.numeric(result$lambda.min$lambda1))

  # index points into grid
  expect_equal(result$cvm[result$index], min(result$cvm, na.rm = TRUE))

  # coxtrans.fit is valid
  expect_true(is.matrix(result$coxtrans.fit$coefficients))
  expect_equal(nrow(result$coxtrans.fit$coefficients), 20L)

  # delegation methods
  expect_equal(coef(result), coef(result$coxtrans.fit))

  # print works
  expect_output(print(result), "cv.coxtrans")
})

test_that("lambda max values are computed on separate penalty scales", {
  formula <- Surv(time, status) ~ . - group - id
  bounds <- survtrans:::coxtrans_penalty_bounds(
    formula, sim2, sim2$group, target = 1
  )
  group_levels <- levels(factor(sim2$group))
  group_releveled <- factor(sim2$group, levels = rev(group_levels))
  global_max <- max(vapply(group_levels, function(target) {
    b <- survtrans:::coxtrans_penalty_bounds(
      formula, sim2, sim2$group, target = target
    )
    max(b$lambda1, b$lambda2, b$lambda3)
  }, numeric(1)))

  expect_true(
    is.numeric(bounds$lambda1) &&
      length(bounds$lambda1) == 1L &&
      is.finite(bounds$lambda1)
  )
  expect_true(
    is.numeric(bounds$lambda2) &&
      length(bounds$lambda2) == 1L &&
      is.finite(bounds$lambda2)
  )
  expect_true(
    is.numeric(bounds$lambda3) &&
      length(bounds$lambda3) == 1L &&
      is.finite(bounds$lambda3)
  )
  expect_true(bounds$lambda1 >= 0)
  expect_true(bounds$lambda2 >= 0)
  expect_true(all(bounds$lambda3 >= 0))
  expect_equal(
    survtrans:::calc_lambda_max(formula, sim2, sim2$group),
    global_max,
    tolerance = 1e-8
  )
  expect_equal(
    survtrans:::calc_lambda_max(formula, sim2, group_releveled),
    survtrans:::calc_lambda_max(formula, sim2, sim2$group),
    tolerance = 1e-8
  )
})

test_that("lambda1_max is large enough to zero the target sparse block", {
  formula <- Surv(time, status) ~ . - group - id
  bounds <- survtrans:::coxtrans_penalty_bounds(
    formula, sim2, sim2$group, target = 1
  )

  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = bounds$lambda1, lambda2 = 0, lambda3 = 0, penalty = "lasso"
  )

  expect_true(all(fit$coefficients[, "1"] == 0))
})

test_that(
  "lambda2_max is large enough to merge all source-target differences",
  {
    formula <- Surv(time, status) ~ . - group - id
    bounds <- survtrans:::coxtrans_penalty_bounds(
      formula, sim2, sim2$group, target = 1
    )

    fit <- coxtrans(
      formula, sim2, sim2$group, 1,
      lambda1 = 0, lambda2 = bounds$lambda2, lambda3 = 0, penalty = "lasso"
    )

    expect_true(all(fit$active_local))
  }
)

test_that("lambda3_max is large enough to activate the prior block", {
  formula <- Surv(time, status) ~ . - group - id
  bounds <- survtrans:::coxtrans_penalty_bounds(
    formula, sim2, sim2$group, target = 1
  )

  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 0, lambda2 = 0, lambda3 = bounds$lambda3, penalty = "lasso"
  )

  expect_true(all(fit$active_prior[, 1]))
})

test_that(
  "coxtrans reduces to target-only fit when lambda2 and lambda3 are zero",
  {
    formula <- Surv(time, status) ~ . - group - id
    target_data <- sim2[sim2$group == 1, ]

    fit_transfer <- coxtrans(
      formula, sim2, sim2$group, 1,
      lambda1 = 0, lambda2 = 0, lambda3 = 0, penalty = "SCAD"
    )
    fit_target <- ncvcox(
      formula, target_data, target_data$group,
      lambda = 0, penalty = "SCAD"
    )

    xt <- as.matrix(
      target_data[, paste0("X", seq_len(nrow(fit_transfer$coefficients)))]
    )
    lp_transfer <- as.numeric(xt %*% fit_transfer$coefficients[, "1"])
    lp_target <- as.numeric(xt %*% fit_target$coefficients)
    expect_true(stats::cor(lp_transfer, lp_target) > 0.98)
  }
)

test_that("large lambda2 and lambda3 push coxtrans toward a shared fit", {
  formula <- Surv(time, status) ~ . - group - id
  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 0, lambda2 = 1e3, lambda3 = 1e3, penalty = "lasso"
  )

  coef_spread <- apply(fit$coefficients, 1, function(row) max(row) - min(row))
  expect_true(max(abs(coef_spread)) < 1e-4)
})

test_that("refit.coxtrans produces debiased estimates", {
  formula <- Surv(time, status) ~ . - group - id
  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
  )

  r <- refit(fit)

  expect_s3_class(r, "refit.coxtrans")
  expect_s3_class(r$coxph_fit, "coxph")
  expect_true(is.matrix(r$coefficients))
  expect_equal(dim(r$coefficients), c(20L, 5L))
  expect_equal(
    colnames(r$coefficients),
    c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
  )

  # Active features should have p < 1
  active <- which(abs(fit$coefficients[, 1]) > 1e-6)
  expect_true(all(r$coefficients[active, "Pr(>|z|)"] < 1))

  # Inactive features should have p = 1
  inactive <- which(abs(fit$coefficients[, 1]) <= 1e-6)
  expect_true(all(r$coefficients[inactive, "Pr(>|z|)"] == 1))

  # print works
  expect_output(print(r), "Refit coxtrans")
})

test_that("refit.coxtrans handles empty active set", {
  formula <- Surv(time, status) ~ . - group - id
  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 10, lambda2 = 0, lambda3 = 0, penalty = "SCAD"
  )

  r <- refit(fit)
  expect_s3_class(r, "refit.coxtrans")
  expect_null(r$coxph_fit)
  expect_equal(length(r$active_set), 0L)
  expect_true(all(r$coefficients[, "Pr(>|z|)"] == 1))
})

test_that("cox_data returns standardized data with correct dimensions", {
  formula <- Surv(time, status) ~ . - group - id
  result <- survtrans:::cox_data(formula, sim2, group = sim2$group)

  expect_true(is.matrix(result$x))
  expect_equal(nrow(result$x), nrow(sim2))
  expect_equal(ncol(result$x), 20L)
  expect_true(!is.null(attr(result$x, "scaled:scale")))
  expect_true(is.numeric(result$time))
  expect_true(is.numeric(result$status))
  expect_true(is.factor(result$group))
})
