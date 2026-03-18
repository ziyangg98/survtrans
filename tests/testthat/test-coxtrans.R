test_that("coxtrans basic fit returns correct structure", {
  formula <- Surv(time, status) ~ . - group - id
  fit <- coxtrans(
    formula, sim2, sim2$group, 1,
    lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
  )

  expect_s3_class(fit, "coxtrans")
  expect_true(is.matrix(fit$coefficients))
  expect_equal(nrow(fit$coefficients), 20L)
  expect_true("Center" %in% colnames(fit$coefficients))
  expect_equal(colnames(fit$coefficients)[1], "1")
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

test_that("ncvcox basic fit returns correct structure", {
  formula <- Surv(time, status) ~ . - group - id
  df <- sim2[sim2$group == 2 | sim2$group == 4, ]
  fit <- ncvcox(formula, df, df$group, lambda = 0.1, penalty = "SCAD")

  expect_s3_class(fit, "ncvcox")
  expect_true(is.numeric(fit$coefficients))
  expect_equal(length(fit$coefficients), 20L)
})

test_that("preprocess returns standardized data with correct dimensions", {
  formula <- Surv(time, status) ~ . - group - id
  result <- survtrans:::preprocess(formula, sim2, group = sim2$group)

  expect_true(is.matrix(result$x))
  expect_equal(nrow(result$x), nrow(sim2))
  expect_equal(ncol(result$x), 20L)
  expect_true(!is.null(attr(result$x, "scaled:scale")))
  expect_true(is.numeric(result$time))
  expect_true(is.numeric(result$status))
  expect_true(is.factor(result$group))
})
