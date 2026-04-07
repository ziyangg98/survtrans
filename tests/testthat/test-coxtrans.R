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
