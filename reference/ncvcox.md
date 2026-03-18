# Non-convex penalized Cox proportional hazards model

Non-convex penalized Cox proportional hazards model

## Usage

``` r
ncvcox(
  formula,
  data,
  group,
  lambda = 0,
  penalty = c("lasso", "MCP", "SCAD"),
  gamma = switch(penalty, SCAD = 3.7, MCP = 3, 1),
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

- lambda:

  A non-negative value specifying the penalty parameter. The default is
  0.

- penalty:

  A character string specifying the penalty function. The default is
  "lasso". Other options are "MCP" and "SCAD".

- gamma:

  A non-negative value specifying the penalty parameter. The default is
  3.7 for SCAD and 3.0 for MCP.

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

An object of class `ncvcox`.

## Examples

``` r
formula <- Surv(time, status) ~ . - group - id
df <- sim2[sim2$group == 2 | sim2$group == 4, ]
fit <- ncvcox(formula, df, df$group, lambda = 0.1, penalty = "SCAD")
summary(fit)
#> Call:
#> ncvcox(formula = formula, data = df, group = df$group, lambda = 0.1, 
#>     penalty = "SCAD")
#> 
#>   n=200, number of events=159
#> 
#>        coef exp(coef) se(coef)     z Pr(>|z|)    
#> X1  0.96399   2.62213  0.09935 9.703  < 2e-16 ***
#> X2  0.98466   2.67691  0.10030 9.817  < 2e-16 ***
#> X3  0.30559   1.35742  0.08740 3.496 0.000472 ***
#> X4  0.39330   1.48187  0.08017 4.906 9.31e-07 ***
#> X13 0.04327   1.04422  0.08511 0.508 0.611185    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#>     exp(coef) exp(-coef) lower .95 upper .95
#> X1  2.6221    0.3814     2.1582    3.1858   
#> X2  2.6769    0.3736     2.1992    3.2584   
#> X3  1.3574    0.7367     1.1437    1.6111   
#> X4  1.4819    0.6748     1.2664    1.7340   
#> X13 1.0442    0.9577     0.8838    1.2338   
```
