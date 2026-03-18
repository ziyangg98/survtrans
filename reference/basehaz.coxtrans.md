# Predict the cumulative baseline hazard function for `coxtrans` objects

Predict the cumulative baseline hazard function for `coxtrans` objects

## Usage

``` r
# S3 method for class 'coxtrans'
basehaz(object, ...)
```

## Arguments

- object:

  An object of class `coxtrans`.

- ...:

  Additional arguments (unused).

## Value

A `data.frame` with one row for each time point, and columns containing
the event time, the cumulative baseline hazard function, and the strata.
