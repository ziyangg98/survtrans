# Plot cross-validation curve for a `cv.coxtrans` object

Plots the CV deviance profile along `lambda1` with `lambda2` and
`lambda3` fixed at their optimal values, in the style of
`plot.cv.glmnet`.

## Usage

``` r
# S3 method for class 'cv.coxtrans'
plot(x, ...)
```

## Arguments

- x:

  A `cv.coxtrans` object.

- ...:

  Further graphical arguments passed to `plot`.
