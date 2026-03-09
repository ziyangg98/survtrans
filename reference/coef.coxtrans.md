# Extract the coefficients from a `coxtrans` object

Extract the coefficients from a `coxtrans` object

## Usage

``` r
# S3 method for class 'coxtrans'
coef(object, ...)
```

## Arguments

- object:

  An object of class `coxtrans`.

- ...:

  Additional arguments (unused).

## Value

A named numeric vector containing the coefficients of the fitted
`coxtrans` object. The names indicate the group(s) to which the
coefficients belong. Zero coefficients are removed.
