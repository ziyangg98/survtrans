# Summary method for a `coxmtl` object

Summary method for a `coxmtl` object

## Usage

``` r
# S3 method for class 'coxmtl'
summary(object, conf.int = 0.95, ...)
```

## Arguments

- object:

  An object of class `coxmtl`.

- conf.int:

  A numeric value between 0 and 1 indicating the confidence level of the
  confidence interval. Default is 0.95.

- ...:

  Additional arguments (unused).

## Value

An object of class `summary.coxmtl` with model size, event count,
log-likelihood, coefficient table, confidence intervals, and center
weights.
