# Prediction method for `ncvcox` objects.

Prediction method for `ncvcox` objects.

## Usage

``` r
# S3 method for class 'ncvcox'
predict(object, newdata, type = c("lp", "risk"), ...)
```

## Arguments

- object:

  An object of class `ncvcox`.

- newdata:

  New data for making predictions.

- type:

  The type of prediction to perform. Options include:

  `"lp"`

  :   The linear predictor.

  `"risk"`

  :   The risk score \\\exp(\text{lp})\\.

- ...:

  Additional arguments (unused).

## Value

A numeric vector of predictions.
