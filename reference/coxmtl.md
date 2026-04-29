# Symmetric Multi-Task Cox Model with Global Centers

Fits a target-free multi-task Cox model across all groups using three
symmetric penalties: group-wise sparsity, pairwise fusion, and shrinkage
toward user-defined global centers.

## Usage

``` r
coxmtl(
  formula,
  data,
  group,
  w,
  lambda1 = 0,
  lambda2 = 0,
  lambda3 = 0,
  penalty = c("lasso", "MCP", "SCAD"),
  gamma = switch(penalty, SCAD = 3.7, MCP = 3, 1),
  vartheta = 1,
  control,
  ...
)
```

## Arguments

- formula:

  A formula with a [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)
  response.

- data:

  A data frame containing the variables in the model.

- group:

  A factor indicating the group of each observation.

- w:

  A non-negative matrix with one row per global center and one column
  per group. Each row must sum to one and defines a convex combination
  of the group-specific coefficients.

- lambda1:

  Sparse penalty applied to every group-specific coefficient.

- lambda2:

  Fusion penalty applied to every pairwise group difference.

- lambda3:

  Center penalty (scalar or vector of length `nrow(w)`).

- penalty:

  Penalty type: `"lasso"`, `"MCP"`, or `"SCAD"`.

- gamma:

  Concavity parameter for MCP/SCAD. Default 3.7 (SCAD) or 3.0 (MCP).

- vartheta:

  Fixed augmented Lagrangian parameter. Default 1.0.

- control:

  A
  [survtrans_control](http://gongziyang.com/survtrans/reference/survtrans_control.md)
  object.

- ...:

  Additional arguments passed to `survtrans_control`.

## Value

An object of class `coxmtl`.
