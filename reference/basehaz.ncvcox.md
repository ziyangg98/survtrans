# Predict the cumulative baseline hazard function for `ncvcox` objects

Predict the cumulative baseline hazard function for `ncvcox` objects

## Usage

``` r
# S3 method for class 'ncvcox'
basehaz(object, newdata, ...)
```

## Arguments

- object:

  An object of class `ncvcox`.

- newdata:

  A numeric vector of time points at which to predict the baseline
  hazard function. If `NULL`, the function will predict the baseline
  hazard function at the unique event times in the fitted data.

- ...:

  Additional arguments (unused).

## Value

A `data.frame` with one row for each time point, and columns containing
the event time, the cumulative baseline hazard function, and the strata.
