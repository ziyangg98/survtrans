# Predict the cumulative baseline hazard function for `ncvcox` objects

Predict the cumulative baseline hazard function for `ncvcox` objects

## Usage

``` r
# S3 method for class 'ncvcox'
basehaz(object, ...)
```

## Arguments

- object:

  An object of class `ncvcox`.

- ...:

  Additional arguments (unused).

## Value

A `data.frame` with one row for each time point, and columns containing
the event time, the cumulative baseline hazard function, and the strata.
