# Ancillary arguments for controlling survtrans fitting

Ancillary arguments for controlling survtrans fitting

## Usage

``` r
survtrans_control(
  abstol = 1e-04,
  reltol = 0.001,
  fdev = 1e-05,
  maxit = 300,
  verbose = FALSE
)
```

## Arguments

- abstol:

  the absolute tolerance for ADMM primal/dual residuals. Default is
  1e-4.

- reltol:

  the relative tolerance for ADMM primal/dual residuals. Default is
  1e-3.

- fdev:

  the minimum fractional change of the augmented Lagrangian for
  convergence. The algorithm stops when \\\|L^{k} - L^{k-1}\| /
  (\|L^{k-1}\| + 1) \< fdev\\, where \\L\\ is the augmented Lagrangian.
  This provides a fallback when primal-dual residuals oscillate under
  non-convex penalties. Default is 1e-5.

- maxit:

  the maximum number of iterations for the proposed algorithm. Default
  is 300.

- verbose:

  a logical value indicating whether to print messages during the
  fitting process. Default is `FALSE`.

## Value

A list with components `abstol`, `reltol`, `fdev`, `maxit`, and
`verbose`.
