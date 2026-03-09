# Diagnose Cox Transfer Model's Optimization Process

Diagnose Cox Transfer Model's Optimization Process

## Usage

``` r
# S3 method for class 'coxtrans'
diagnose(object, ...)
```

## Arguments

- object:

  An object of class `coxtrans`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `NULL`. Called for its side effect of producing
diagnostic plots.

## Details

This function produces two plots:

- Residuals Convergence: Plots the evolution of primal and dual
  residuals along with their tolerance levels.

- Loss Decomposition: Plots the negative log-likelihood, total loss, and
  penalty term.
