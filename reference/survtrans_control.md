# Ancillary arguments for controlling survtrans fitting

Ancillary arguments for controlling survtrans fitting

## Usage

``` r
survtrans_control(
  abstol = 1e-04,
  reltol = 0.001,
  nthreads = 1L,
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

- nthreads:

  the number of OpenMP threads used by a single model fit. Default is 1.

- maxit:

  the maximum number of iterations for the proposed algorithm. Default
  is 300.

- verbose:

  a logical value indicating whether to print messages during the
  fitting process. Default is `FALSE`.

## Value

A list with components `abstol`, `reltol`, `nthreads`, `maxit`, and
`verbose`.
