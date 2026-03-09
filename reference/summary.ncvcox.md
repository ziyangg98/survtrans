# Summary method for a `ncvcox` object

Summary method for a `ncvcox` object

## Usage

``` r
# S3 method for class 'ncvcox'
summary(object, conf.int = 0.95, compressed = TRUE, ...)
```

## Arguments

- object:

  An object of class `ncvcox`.

- conf.int:

  A numeric value between 0 and 1 indicating the confidence level of the
  confidence interval. Default is 0.95.

- compressed:

  Logical; if `TRUE`, the summary is compressed and only includes the
  group-level coefficients. Default is `TRUE`.

- ...:

  Additional arguments (not used).

## Value

An object of class `summary.ncvcox`, with the following components:

- `n`, `nevent`:

  Number of observations and number of events, respectively, in the fit.

- `logLik`:

  The log partial likelihood at the final value.

- `BIC`:

  The Bayesian Information Criterion at the final value.

- `coefficients`:

  A matrix with one row for each coefficient, and columns containing the
  coefficient, the hazard ratio exp(coef), standard error, Wald
  statistic, and P value.

- `conf.int`:

  A matrix with one row for each coefficient, containing the confidence
  limits for exp(coef).
