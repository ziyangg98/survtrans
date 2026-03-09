# Prediction method for `coxtrans` objects.

Prediction method for `coxtrans` objects.

## Usage

``` r
# S3 method for class 'coxtrans'
predict(object, newdata = NULL, newgroup = NULL, type = c("lp", "risk"), ...)
```

## Arguments

- object:

  An object of class `coxtrans`.

- newdata:

  Optional new data for making predictions. If omitted, predictions are
  made using the data used for fitting the model.

- newgroup:

  Optional new group for making predictions. If omitted, predictions are
  made using the groups from the original data.

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
