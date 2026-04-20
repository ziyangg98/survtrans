# Transfer Learning Cox Model with Prior Constraints

Fits a Cox proportional hazards model that transfers survival
information from source domains to a target domain using three-layer
penalization: sparse (lambda1), local sharing (lambda2), and prior
transfer (lambda3).

## Usage

``` r
coxtrans(
  formula,
  data,
  group,
  target,
  lambda1 = 0,
  lambda2 = 0,
  lambda3 = 0,
  prior_matrix = NULL,
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

- target:

  The target group level.

- lambda1:

  Sparse penalty (non-negative scalar). Default 0.

- lambda2:

  Local penalty for source-target coefficient differences (non-negative
  scalar). Default 0.

- lambda3:

  Prior transfer penalty (non-negative scalar or vector of length G).
  When `prior_matrix` has G rows, `lambda3` should have length G; a
  scalar is recycled. Default 0.

- prior_matrix:

  Optional G x (K-1) weight matrix defining transfer constraints. Each
  row specifies a weighted combination of source coefficients that the
  target should be close to. If `NULL`, a single prior using
  sample-weighted source means is used.

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

An object of class `coxtrans` with components:

- coefficients:

  Matrix (p x K) of group-specific coefficients. Each column is the full
  beta for that group (target first).

- prior_matrix:

  The prior weight matrix used.

- prior_effects:

  Matrix of prior constraint residuals.

- active_local:

  Logical matrix (p x K-1) of active local constraints.

- active_prior:

  Logical matrix (p x G) of active prior constraints.

- iter:

  Number of ADMM iterations.

- message:

  Convergence message.

- history:

  Matrix of per-iteration diagnostics.

## Examples

``` r
formula <- Surv(time, status) ~ . - group - id

# Basic usage with default prior (sample-weighted source mean)
fit <- coxtrans(
  formula, sim2, sim2$group, 1,
  lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
)
summary(fit)
#> Call:
#> coxtrans(formula = formula, data = sim2, group = sim2$group, 
#>     target = 1, lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, 
#>     penalty = "SCAD")
#> 
#>   n=500, number of events=422
#> 
#>       coef exp(coef) se(coef)     z Pr(>|z|)    
#> X1 0.34587   1.41322  0.05340 6.477 9.34e-11 ***
#> X2 0.35601   1.42762  0.05403 6.589 4.44e-11 ***
#> X3 0.34327   1.40956  0.05396 6.362 1.99e-10 ***
#> X4 0.32658   1.38622  0.05155 6.335 2.37e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.4132    0.7076     1.2728    1.5691   
#> X2 1.4276    0.7005     1.2842    1.5871   
#> X3 1.4096    0.7094     1.2681    1.5668   
#> X4 1.3862    0.7214     1.2530    1.5336   
#> 
#> Feature structure:
#>   Prior transfer: X1, X2
#>   Shared (local): X3, X4
#>   Sparse (zero) : 16 features (X5, X6, X7, ...)

# Custom prior matrix (two tissue-based priors)
pm <- rbind(
  tissue_A = c(0.5, 0.5, 0, 0),
  tissue_B = c(0, 0, 0.5, 0.5)
)
colnames(pm) <- c("2", "3", "4", "5")
fit2 <- coxtrans(
  formula, sim2, sim2$group, 1,
  lambda1 = 0.075, lambda2 = 0.04, lambda3 = c(0.04, 0.04),
  prior_matrix = pm, penalty = "SCAD"
)
summary(fit2)
#> Call:
#> coxtrans(formula = formula, data = sim2, group = sim2$group, 
#>     target = 1, lambda1 = 0.075, lambda2 = 0.04, lambda3 = c(0.04, 
#>         0.04), prior_matrix = pm, penalty = "SCAD")
#> 
#>   n=500, number of events=422
#> 
#>       coef exp(coef) se(coef)     z Pr(>|z|)    
#> X1 0.34289   1.40901  0.06213 5.519 3.40e-08 ***
#> X2 0.35743   1.42965  0.06567 5.443 5.24e-08 ***
#> X3 0.34416   1.41081  0.05408 6.364 1.96e-10 ***
#> X4 0.32464   1.38354  0.05097 6.370 1.90e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.4090    0.7097     1.2475    1.5915   
#> X2 1.4296    0.6995     1.2570    1.6260   
#> X3 1.4108    0.7088     1.2689    1.5685   
#> X4 1.3835    0.7228     1.2520    1.5289   
#> 
#> Feature structure:
#>   Prior [tissue_A,tissue_B]: X1, X2
#>   Shared (local)           : X3, X4
#>   Sparse (zero)            : 16 features (X5, X6, X7, ...)
```
