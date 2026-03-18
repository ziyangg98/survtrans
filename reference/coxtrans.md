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
  vartheta = 1,
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

- vartheta:

  A positive value specifying the fixed penalty parameter in the
  augmented Lagrangian. Following Wang, Yin & Zeng (2019), this is kept
  constant (not adaptive) to guarantee convergence with non-convex
  penalties. The default is 1.0.

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
formula <- Surv(time, status) ~ . - group - id
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
#> X1 0.33676   1.40040  0.05341 6.306 2.87e-10 ***
#> X2 0.35968   1.43287  0.05402 6.659 2.76e-11 ***
#> X3 0.34368   1.41012  0.05398 6.367 1.93e-10 ***
#> X4 0.32553   1.38476  0.05157 6.313 2.75e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.4004    0.7141     1.2612    1.5549   
#> X2 1.4329    0.6979     1.2889    1.5929   
#> X3 1.4101    0.7092     1.2686    1.5675   
#> X4 1.3848    0.7221     1.2516    1.5320   
```
