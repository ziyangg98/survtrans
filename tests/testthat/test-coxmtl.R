test_that("coxmtl basic fit returns expected structure", {
  formula <- Surv(time, status) ~ . - group - id
  w <- rbind(
    global = rep(0.2, 5),
    anchor1 = c(1, 0, 0, 0, 0)
  )
  colnames(w) <- levels(factor(sim2$group))

  fit <- coxmtl(
    formula, sim2, sim2$group, w,
    lambda1 = 0.05, lambda2 = 0.02, lambda3 = c(0.02, 0.01),
    penalty = "SCAD"
  )

  expect_s3_class(fit, "coxmtl")
  expect_true(is.matrix(fit$coefficients))
  expect_equal(dim(fit$coefficients), c(20L, 5L))
  expect_true(is.numeric(fit$phi))
  expect_true(is.matrix(fit$basis))
  expect_equal(ncol(fit$w), 5L)
  expect_equal(rownames(fit$active_sparse), rownames(fit$coefficients))
  expect_equal(colnames(fit$active_sparse), colnames(fit$coefficients))
})

test_that("coxmtl validates W", {
  formula <- Surv(time, status) ~ . - group - id
  bad_w <- rbind(c(0.5, 0.5, 0, 0, 0.2))

  expect_error(
    coxmtl(formula, sim2, sim2$group, bad_w, lambda1 = 0.1),
    "sum to 1"
  )
  expect_error(
    coxmtl(formula, sim2, sim2$group, matrix(-1, 1, 5), lambda1 = 0.1),
    "non-negative"
  )
})

test_that("coxmtl S3 methods work", {
  formula <- Surv(time, status) ~ . - group - id
  w <- matrix(rep(1 / 5, 5), nrow = 1)
  colnames(w) <- levels(factor(sim2$group))

  fit <- coxmtl(
    formula, sim2, sim2$group, w,
    lambda1 = 0.03, lambda2 = 0.01, lambda3 = 0.01, penalty = "SCAD"
  )

  expect_true(is.numeric(coef(fit)))
  expect_true(is.numeric(predict(fit, newdata = sim2, newgroup = sim2$group)))
  expect_true(is.numeric(logLik(fit)))
  expect_true(is.list(basehaz(fit)))
  expect_true(is.matrix(vcov(fit)))
  expect_s3_class(summary(fit), "summary.coxmtl")
  expect_output(print(fit), "Global centers")
})

test_that("coxmtl predict uses rhs-only newdata and replays centering", {
  formula <- Surv(time, status) ~ . - group - id
  w <- matrix(rep(1 / 5, 5), nrow = 1)
  colnames(w) <- levels(factor(sim2$group))

  fit <- coxmtl(
    formula, sim2, sim2$group, w,
    lambda1 = 0.03, lambda2 = 0.01, lambda3 = 0.01, penalty = "SCAD"
  )

  ord <- order(sim2$time, decreasing = TRUE)
  newdata <- sim2[ord, setdiff(names(sim2), c("time", "status")), drop = FALSE]
  newgroup <- sim2$group[ord]

  expect_equal(
    predict(fit, newdata = newdata, newgroup = newgroup),
    predict(fit),
    tolerance = 1e-8
  )
  expect_equal(
    predict(fit, newdata = newdata),
    predict(fit),
    tolerance = 1e-8
  )
})

test_that("coxmtl predict rejects unseen groups", {
  formula <- Surv(time, status) ~ . - group - id
  w <- matrix(rep(1 / 5, 5), nrow = 1)
  colnames(w) <- levels(factor(sim2$group))

  fit <- coxmtl(
    formula, sim2, sim2$group, w,
    lambda1 = 0.03, lambda2 = 0.01, lambda3 = 0.01, penalty = "SCAD"
  )

  newdata <- sim2[, setdiff(names(sim2), c("time", "status")), drop = FALSE]
  bad_group <- as.character(sim2$group)
  bad_group[1] <- "typo"

  expect_error(
    predict(fit, newdata = newdata, newgroup = bad_group),
    "unseen group labels"
  )
})

test_that("coxmtl drops unused factor levels before deriving tasks", {
  formula <- Surv(time, status) ~ . - group - id
  df <- sim2[sim2$group %in% c(1, 2), ]
  df$group <- factor(df$group, levels = levels(factor(sim2$group)))
  w <- matrix(c(0.5, 0.5), nrow = 1)
  colnames(w) <- c("1", "2")

  fit <- coxmtl(
    formula, df, df$group, w,
    lambda1 = 0.01, lambda2 = 0.01, lambda3 = 0.01, penalty = "SCAD"
  )

  expect_equal(colnames(fit$coefficients), c("1", "2"))
  expect_equal(ncol(fit$coefficients), 2L)
})

test_that("large lambda1 zeros every group coefficient", {
  formula <- Surv(time, status) ~ . - group - id
  w <- matrix(rep(1 / 5, 5), nrow = 1)
  colnames(w) <- levels(factor(sim2$group))

  fit <- coxmtl(
    formula, sim2, sim2$group, w,
    lambda1 = 1e3, lambda2 = 0, lambda3 = 0, penalty = "lasso"
  )

  expect_true(all(fit$coefficients == 0))
})

test_that("large lambda2 collapses all groups to a shared profile", {
  formula <- Surv(time, status) ~ . - group - id
  w <- matrix(rep(1 / 5, 5), nrow = 1)
  colnames(w) <- levels(factor(sim2$group))

  fit <- coxmtl(
    formula, sim2, sim2$group, w,
    lambda1 = 0, lambda2 = 1e3, lambda3 = 0, penalty = "lasso"
  )

  coef_spread <- apply(fit$coefficients, 1, function(row) max(row) - min(row))
  expect_true(max(abs(coef_spread)) < 1e-4)
})

test_that("large lambda3 enforces closeness to w-defined centers", {
  formula <- Surv(time, status) ~ . - group - id
  w <- matrix(c(1, 0, 0, 0, 0), nrow = 1)
  colnames(w) <- levels(factor(sim2$group))

  fit <- coxmtl(
    formula, sim2, sim2$group, w,
    lambda1 = 0, lambda2 = 0, lambda3 = 1e3, penalty = "lasso"
  )

  residual <- sweep(fit$coefficients, 1, fit$coefficients[, 1], "-")
  expect_true(max(abs(residual)) < 1e-4)
})

test_that("BIC.coxmtl uses log(nevent) times length(phi)", {
  formula <- Surv(time, status) ~ . - group - id
  w <- matrix(rep(1 / 5, 5), nrow = 1)
  colnames(w) <- levels(factor(sim2$group))

  fit <- coxmtl(
    formula, sim2, sim2$group, w,
    lambda1 = 0.02, lambda2 = 0.01, lambda3 = 0.01, penalty = "SCAD"
  )

  expect_equal(
    BIC(fit),
    -2 * logLik(fit) + log(sum(sim2$status)) * length(fit$phi),
    tolerance = 1e-8
  )
})

test_that("single-group coxmtl reduces to single-group penalized Cox", {
  formula <- Surv(time, status) ~ . - group - id
  df <- sim2[sim2$group == 1, ]
  w <- matrix(1, nrow = 1, ncol = 1)
  colnames(w) <- "1"

  fit_mtl <- coxmtl(
    formula, df, df$group, w,
    lambda1 = 0, lambda2 = 0, lambda3 = 0, penalty = "SCAD"
  )
  fit_single <- ncvcox(
    formula, df, df$group,
    lambda = 0, penalty = "SCAD"
  )

  expect_true(
    stats::cor(
      as.numeric(fit_mtl$coefficients[, 1]),
      as.numeric(fit_single$coefficients)
    ) > 0.98
  )
})
