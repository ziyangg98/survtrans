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
formula <- survival::Surv(time, status) ~ . - group - id
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
#> X1  0.89714   2.45257  0.09238 9.712  < 2e-16 ***
#> X2  1.00976   2.74493  0.10344 9.762  < 2e-16 ***
#> X3  0.31636   1.37212  0.09135 3.463 0.000534 ***
#> X4  0.42938   1.53630  0.08829 4.863 1.15e-06 ***
#> X13 0.04271   1.04364  0.08450 0.505 0.613225    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#>     exp(coef) exp(-coef) lower .95 upper .95
#> X1  2.4526    0.4077     2.0464    2.9394   
#> X2  2.7449    0.3643     2.2412    3.3618   
#> X3  1.3721    0.7288     1.1472    1.6412   
#> X4  1.5363    0.6509     1.2922    1.8265   
#> X13 1.0436    0.9582     0.8844    1.2316   
```
