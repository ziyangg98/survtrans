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
  vartheta = NULL,
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

  Initial augmented Lagrangian parameter. If `NULL`, lasso uses 1.0 and
  non-convex penalties use a curvature-safe value.

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
#> X1 0.35055   1.41985  0.05352 6.550 5.74e-11 ***
#> X2 0.35887   1.43171  0.05418 6.623 3.51e-11 ***
#> X3 0.34424   1.41091  0.05394 6.381 1.76e-10 ***
#> X4 0.32870   1.38916  0.05164 6.365 1.95e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.4198    0.7043     1.2785    1.5769   
#> X2 1.4317    0.6985     1.2875    1.5921   
#> X3 1.4109    0.7088     1.2694    1.5683   
#> X4 1.3892    0.7199     1.2554    1.5371   
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
#> X1 0.34612   1.41357  0.06219 5.565 2.62e-08 ***
#> X2 0.35968   1.43287  0.06571 5.473 4.42e-08 ***
#> X3 0.34534   1.41247  0.05409 6.384 1.72e-10 ***
#> X4 0.32688   1.38663  0.05105 6.403 1.52e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.4136    0.7074     1.2514    1.5968   
#> X2 1.4329    0.6979     1.2597    1.6298   
#> X3 1.4125    0.7080     1.2704    1.5705   
#> X4 1.3866    0.7212     1.2546    1.5326   
#> 
#> Feature structure:
#>   Prior [tissue_A,tissue_B]: X1, X2
#>   Shared (local)           : X3, X4
#>   Sparse (zero)            : 16 features (X5, X6, X7, ...)
```
