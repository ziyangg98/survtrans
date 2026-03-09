# Compute penalty for a numeric vector using specified regularization type.

This function computes the penalty for a given parameter vector `x`
based on the specified penalty approach.

## Usage

``` r
penalty(x, penalty, lambda, gamma)
```

## Arguments

- x:

  A numeric vector for which the penalty is computed.

- penalty:

  A character string specifying the penalty type. Valid options are
  "lasso", "MCP", or "SCAD".

- lambda:

  A numeric value representing the regularization parameter.

- gamma:

  A numeric value used in the penalty for MCP and SCAD.

## Value

A numeric value representing the computed penalty.

## Details

The computation differs according to the penalty type:

- lasso:

  Returns the L1 norm of `x` scaled by `lambda`.

- MCP:

  Applies a minimax concave penalty where the penalty function behaves
  linearly when the absolute values are below `lambda * gamma` and
  quadratically otherwise.

- SCAD:

  Applies the smoothly clipped absolute deviation penalty; different
  expressions are used when `abs(x)` is below `lambda`, between `lambda`
  and `gamma * lambda`, or above these thresholds.
