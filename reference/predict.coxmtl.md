# Prediction method for `coxmtl` objects

Prediction method for `coxmtl` objects

## Usage

``` r
# S3 method for class 'coxmtl'
predict(object, newdata = NULL, newgroup = NULL, type = c("lp", "risk"), ...)
```

## Arguments

- object:

  An object of class `coxmtl`.

- newdata:

  Optional new data for making predictions. If omitted, predictions are
  made using the data used for fitting the model.

- newgroup:

  Optional new group labels for making predictions. If omitted,
  predictions use the groups from the original data, or the `group`
  column in `newdata` when available.

- type:

  The type of prediction to perform: `"lp"` for the linear predictor or
  `"risk"` for \\\exp(\text{lp})\\.

- ...:

  Additional arguments (unused).

## Value

A numeric vector of predictions.
