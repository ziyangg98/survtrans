# Print method for a `summary.ncvcox` object

This function prints a summary of the results of a `ncvcox` model in a
formatted and user-friendly manner, including the model call, number of
samples, number of events, coefficients, and confidence intervals. It
also includes significance stars for p-values.

## Usage

``` r
# S3 method for class 'summary.ncvcox'
print(
  x,
  digits = max(getOption("digits") - 3, 3),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- x:

  A summary object produced from a fitted `ncvcox` model. This object
  contains information such as model coefficients and confidence
  intervals.

- digits:

  An integer controlling the number of significant digits to print for
  numeric values.

- signif.stars:

  Logical; if `TRUE`, significance stars are printed along with the
  p-values.

- ...:

  Additional arguments (unused).

## Value

The function prints the summary of the `ncvcox` model and returns the
object `x` invisibly.

## Details

The function provides a formatted output that includes:

- **Call:** The original function call that produced the model.

- **n and events:** The total number of samples and the number of events
  (e.g., deaths).

- **Coefficients:** The regression coefficients, their standard errors,
  z-values, and p-values, formatted in a table. Significance stars are
  shown next to p-values if `signif.stars` is `TRUE`.

- **Confidence intervals:** The exponentiated coefficients along with
  their confidence intervals.
