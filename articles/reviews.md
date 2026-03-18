# A Review of Transfer Learning Approaches for Survival Analysis

``` r
library(survtrans)
```

## Introduction

Consider a scenario where we have survival data from multiple source
domains and we want to transfer the survival information to a target
domain. In this vignette, we will review the transfer learning
approaches for survival analysis.

``` r
formula <- survival::Surv(time, status) ~ . - group - id
df_target <- sim2[sim2$group == 1, ]
```

## Target Learning

``` r
fit_target <- ncvcox(formula, df_target, lambda = 0.2, penalty = "SCAD")
basehaz_target <- basehaz(fit_target)
summary(fit_target)
#> Call:
#> ncvcox(formula = formula, data = df_target, lambda = 0.2, penalty = "SCAD")
#> 
#>   n=100, number of events=83
#> 
#>        coef exp(coef) se(coef)     z Pr(>|z|)
#> X1 0.053289  1.054734 0.108572 0.491    0.624
#> X2 0.131001  1.139969 0.106658 1.228    0.219
#> X3 0.002069  1.002071 0.115958 0.018    0.986
#> X4 0.060811  1.062698 0.116182 0.523    0.601
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.0547    0.9481     0.8526    1.3048   
#> X2 1.1400    0.8772     0.9249    1.4050   
#> X3 1.0021    0.9979     0.7984    1.2578   
#> X4 1.0627    0.9410     0.8463    1.3345
```

## Global Learning

``` r
fit_global <- ncvcox(formula, sim2, lambda = 0.1, penalty = "SCAD")
basehaz_global <- basehaz(fit_global)
summary(fit_global)
#> Call:
#> ncvcox(formula = formula, data = sim2, lambda = 0.1, penalty = "SCAD")
#> 
#>   n=500, number of events=422
#> 
#>       coef exp(coef) se(coef)     z Pr(>|z|)    
#> X1 0.08777   1.09173  0.04847 1.811 0.070182 .  
#> X2 0.18851   1.20745  0.04991 3.777 0.000159 ***
#> X3 0.18674   1.20532  0.05117 3.650 0.000263 ***
#> X4 0.13231   1.14147  0.04542 2.913 0.003581 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.0917    0.9160     0.9928    1.2005   
#> X2 1.2075    0.8282     1.0949    1.3315   
#> X3 1.2053    0.8297     1.0903    1.3325   
#> X4 1.1415    0.8761     1.0442    1.2478
```

### Stratified Learning

``` r
fit_strat <- ncvcox(formula, sim2, sim2$group, lambda = 0.1, penalty = "SCAD")
basehaz_strat <- basehaz(fit_strat)
summary(fit_strat)
#> Call:
#> ncvcox(formula = formula, data = sim2, group = sim2$group, lambda = 0.1, 
#>     penalty = "SCAD")
#> 
#>   n=500, number of events=422
#> 
#>       coef exp(coef) se(coef)     z Pr(>|z|)    
#> X1 0.04284   1.04377  0.05381 0.796 0.425969    
#> X2 0.08578   1.08957  0.05461 1.571 0.116210    
#> X3 0.19885   1.22000  0.05272 3.772 0.000162 ***
#> X4 0.13182   1.14091  0.04603 2.864 0.004189 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.0438    0.9581     0.9393    1.1599   
#> X2 1.0896    0.9178     0.9790    1.2127   
#> X3 1.2200    0.8197     1.1002    1.3528   
#> X4 1.1409    0.8765     1.0425    1.2486
```
