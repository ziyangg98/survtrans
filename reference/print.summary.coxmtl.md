# Print method for a `summary.coxmtl` object

Print method for a `summary.coxmtl` object

## Usage

``` r
# S3 method for class 'summary.coxmtl'
print(
  x,
  digits = max(getOption("digits") - 3, 3),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- x:

  A summary object produced from a fitted `coxmtl` model.

- digits:

  An integer controlling the number of significant digits to print for
  numeric values.

- signif.stars:

  Logical; if `TRUE`, significance stars are printed along with the
  p-values.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns `x`.
