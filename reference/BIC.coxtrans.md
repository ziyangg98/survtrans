# BIC for a `coxtrans` object

BIC for a `coxtrans` object

## Usage

``` r
# S3 method for class 'coxtrans'
BIC(object, type = c("traditional", "modified"), ...)
```

## Arguments

- object:

  An object of class `coxtrans`.

- type:

  A character string specifying the type of BIC to compute.
  "traditional" corresponds to Cn=1, and "modified" corresponds to
  Cn=log(log(d)).

- ...:

  Additional arguments (unused).

## Value

A numeric value representing the BIC of the fitted `coxtrans` object.
