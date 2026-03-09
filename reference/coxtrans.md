# Transfer Learning for Cox Model with Global and Local Shrinkage

Transfer Learning for Cox Model with Global and Local Shrinkage

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
  penalty = c("lasso", "MCP", "SCAD"),
  gamma = switch(penalty, SCAD = 3.7, MCP = 3, 1),
  rho = 2,
  tau = 10,
  init,
  control,
  ...
)
```

## Arguments

- formula:

  A formula object, with the response on the left of a `~` operator, and
  the terms on the right. The response must be a survival object as
  returned by the Surv function.

- data:

  A data frame containing the variables in the model.

- group:

  A factor variable indicating the group of each observation.

- target:

  A character string specifying the target group.

- lambda1:

  A non-negative value specifying the sparse penalty parameter. The
  default is 0.

- lambda2:

  A non-negative value specifying the global biased penalty parameter.
  The default is 0.

- lambda3:

  A non-negative value specifying the local biased penalty parameter.
  The default is 0.

- penalty:

  A character string specifying the penalty function. The default is
  "lasso". Other options are "MCP" and "SCAD".

- gamma:

  A non-negative value specifying the penalty parameter. The default is
  3.7 for SCAD and 3.0 for MCP.

- rho:

  A value larger than 1 specifying the increase/decrease factor for the
  augmented Lagrangian's penalty parameter. The default is 2.0.

- tau:

  A value larger than 1 specifying the tolerance for the trade-off
  between the primal and dual residuals. The default is 10.0.

- init:

  A numeric vector of initial values for the coefficients. The default
  is a zero vector.

- control:

  An object of class
  [survtrans_control](http://gongziyang.com/survtrans/reference/survtrans_control.md)
  containing control parameters for the fitting algorithm. Default is
  `survtrans_control(...)`.

- ...:

  Additional arguments to be passed to the fitting algorithm.

## Value

An object of class `coxtrans`.

## Examples

``` r
formula <- survival::Surv(time, status) ~ . - group - id
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
#> X1 0.35052   1.41981  0.05366 6.533 6.47e-11 ***
#> X2 0.35914   1.43209  0.05424 6.621 3.56e-11 ***
#> X3 0.34485   1.41178  0.05431 6.350 2.15e-10 ***
#> X4 0.32870   1.38916  0.05190 6.333 2.40e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.4198    0.7043     1.2781    1.5773   
#> X2 1.4321    0.6983     1.2877    1.5927   
#> X3 1.4118    0.7083     1.2692    1.5703   
#> X4 1.3892    0.7199     1.2548    1.5379   
```
