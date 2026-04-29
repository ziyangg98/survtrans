# Cross-validated tuning for coxtrans

Selects penalty parameters by minimising the held-out partial likelihood
deviance over a full `lambda1 x lambda2 x lambda3` grid. Supports both
the `lambda.min` rule (minimum CV deviance) and the `lambda.1se` rule
(most sparse model within one standard error of the minimum).

## Usage

``` r
cv.coxtrans(
  formula,
  data,
  group,
  target,
  prior_matrix = NULL,
  nlambda = 10,
  lambda.min.ratio = 0.001,
  nfolds = 10,
  penalty = c("lasso", "MCP", "SCAD"),
  ncores = 1,
  seed = 42,
  verbose = FALSE,
  control = survtrans_control(),
  ...
)

# S3 method for class 'cv.coxtrans'
coef(object, s = c("lambda.min", "lambda.1se"), ...)

# S3 method for class 'cv.coxtrans'
predict(object, s = c("lambda.min", "lambda.1se"), ...)
```

## Arguments

- formula:

  A formula with a [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)
  response.

- data:

  A data frame.

- group:

  A factor indicating group membership.

- target:

  The target group level.

- prior_matrix:

  Optional G x (K-1) transfer constraint matrix.

- nlambda:

  Number of lambda values per dimension. Default 10.

- lambda.min.ratio:

  Smallest/largest lambda ratio. Default 1e-3.

- nfolds:

  Number of CV folds. Default 10.

- penalty:

  `"lasso"`, `"MCP"`, or `"SCAD"`.

- ncores:

  Number of R worker processes for parallel grid evaluation. Default 1.

- seed:

  Random seed.

- verbose:

  Print progress.

- control:

  A
  [survtrans_control](http://gongziyang.com/survtrans/reference/survtrans_control.md)
  object passed to each
  [`coxtrans()`](http://gongziyang.com/survtrans/reference/coxtrans.md)
  fit. Default
  [`survtrans_control()`](http://gongziyang.com/survtrans/reference/survtrans_control.md).

- ...:

  Passed to
  [`coxtrans`](http://gongziyang.com/survtrans/reference/coxtrans.md).

- object:

  A `cv.coxtrans` object.

- s:

  Which model to predict from: `"lambda.min"` (default) or
  `"lambda.1se"`.

## Value

An object of class `cv.coxtrans` with fields: `grid` (full lambda grid),
`cvm` (mean deviance per grid point), `cvsd`, `cvup`, `cvlo`, `nzero`
(mean non-zero target coefficients), `lambda.min` (optimal lambda list),
`lambda.1se` (1SE-rule lambda list), `index.min`, `index.1se`,
`coxtrans.fit` (final model at lambda.min), `coxtrans.fit.1se` (final
model at lambda.1se), `call`, `name`.

## Examples

``` r
result <- cv.coxtrans(Surv(time, status) ~ . - group - id,
  sim2, sim2$group,
  target = 1, nlambda = 5, penalty = "SCAD"
)
result
#> cv.coxtrans
#> 
#> Call: cv.coxtrans(formula = Surv(time, status) ~ . - group - id, data = sim2, 
#>     group = sim2$group, target = 1, nlambda = 5, penalty = "SCAD")
#> 
#> Measure: Partial Likelihood Deviance
#> 
#> lambda.min:  l1=0.0588  l2=0.0177  l3=0.2154
#>   Deviance:  3.2739 (+/- 0.1211)   Non-zero: 5
#> 
#> lambda.1se:  l1=0.0588  l2=0.0006  l3=2e-04
#>   Deviance:  3.3161 (+/- 0.1089)   Non-zero: 4
```
