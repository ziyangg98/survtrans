# Predict the cumulative baseline hazard function for `coxmtl` objects

Predict the cumulative baseline hazard function for `coxmtl` objects

## Usage

``` r
# S3 method for class 'coxmtl'
basehaz(object, ...)
```

## Arguments

- object:

  An object of class `coxmtl`.

- ...:

  Additional arguments (unused).

## Value

A `data.frame` with event time, cumulative baseline hazard, and strata
columns.
