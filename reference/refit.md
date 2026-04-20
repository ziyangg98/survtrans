# Refit a coxtrans model with hard constraints

Extracts the constraint structure (active set, sharing pattern, prior
constraints) from a penalized coxtrans fit, encodes them as hard
constraints via a link matrix, and solves an unpenalized stratified Cox
model in the reparametrized space. Returns debiased MLE coefficients and
valid standard errors for target-group features.

## Usage

``` r
refit(object, ...)

# S3 method for class 'coxtrans'
refit(object, ...)
```

## Arguments

- object:

  A fitted `coxtrans` object.

- ...:

  Additional arguments (unused).

## Value

An object of class `refit.coxtrans` with components:

- coefficients:

  Matrix (p x 5) with columns `coef`, `exp(coef)`, `se(coef)`, `z`,
  `Pr(>|z|)`. Rows correspond to original features; inactive features
  have all zeros except `Pr(>|z|)` = 1.

- coxph_fit:

  The underlying `coxph` object, or `NULL`.

- n:

  Total number of observations used in refit.

- nevent:

  Number of events.

- active_set:

  Integer vector of active feature indices.

## Details

The refit procedure:

1.  Build the link matrix \\L\\ from the penalized fit, encoding shared
    coefficients, prior constraints, and group-specific parameters.

2.  Construct reparametrized design \\Z = \mathrm{bdiag}(X_1, \ldots,
    X_K) L\_{\[:, \text{nonzero}\]}\\.

3.  Fit `coxph(Surv ~ Z + strata(group))` — unpenalized stratified Cox
    with hard constraints built into the design.

4.  Map estimates back to target coefficients via \\L\\, compute
    standard errors via delta method.
