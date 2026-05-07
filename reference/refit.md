# Refit a coxmtl model with hard constraints (oracle estimator)

Extracts the active constraint structure from a penalized `coxmtl` fit,
encodes it as hard constraints via the QR null-space basis, and fits an
unpenalized stratified Cox model in the reparametrized space. Returns
debiased MLE coefficients and valid standard errors for all groups.

Extracts the constraint structure (active set, sharing pattern, prior
constraints) from a penalized coxtrans fit, encodes them as hard
constraints via a link matrix, and solves an unpenalized stratified Cox
model in the reparametrized space. Returns debiased MLE coefficients and
valid standard errors for target-group features.

## Usage

``` r
# S3 method for class 'coxmtl'
refit(object, ...)

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

An object of class `refit.coxmtl` with components:

- coefficients:

  Matrix (p x K) of oracle refitted coefficients.

- se:

  Matrix (p x K) of standard errors (delta method).

- phi_hat:

  Named numeric vector of free-parameter (phi) estimates.

- coxph_fit:

  The underlying `coxph` object, or `NULL`.

- n:

  Total observations used.

- nevent:

  Number of events.

- active_features:

  Integer vector: features non-zero in any group.

- basis:

  The null-space basis \\B\\ from the penalized fit.

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

Oracle refit procedure:

1.  Use the null-space basis \\B\\ stored in the penalized fit, which
    encodes active sparsity, pairwise fusion, and center constraints.

2.  Build reparametrized design \\Z =
    \mathrm{bdiag}(X_1,\ldots,X_K)\\B\\.

3.  Fit `coxph(Surv(time,status) ~ Z + strata(group))` — unpenalized
    stratified Cox with all active constraints hard-coded in the design.

4.  Map back via \\\hat\beta = B\hat\phi\\; compute standard errors via
    the delta method.

The refit procedure:

1.  Build the link matrix \\L\\ from the penalized fit, encoding shared
    coefficients, prior constraints, and group-specific parameters.

2.  Construct reparametrized design \\Z = \mathrm{bdiag}(X_1, \ldots,
    X_K) L\_{\[:, \text{nonzero}\]}\\.

3.  Fit `coxph(Surv ~ Z + strata(group))` — unpenalized stratified Cox
    with hard constraints built into the design.

4.  Map estimates back to target coefficients via \\L\\, compute
    standard errors via delta method.
