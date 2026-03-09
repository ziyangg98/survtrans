# Calculate the maximum value of the penalty parameter lambda

Calculate the maximum value of the penalty parameter lambda

## Usage

``` r
calc_lambda_max(formula, data, group, offset)
```

## Arguments

- formula:

  A formula expression for regression models, in the form
  `response ~ predictors`. The response must be a survival object as
  returned by the [Surv](https://rdrr.io/pkg/survival/man/Surv.html)
  function.

- data:

  A data frame containing the variables in the model.

- group:

  A factor specifying the group of each sample.

- offset:

  A numeric vector specifying the offset.

## Value

The maximum value of the penalty parameter lambda, which shrinks all the
coefficients to zero.
