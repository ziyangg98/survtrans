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
#> X1 0.059452  1.061254 0.117669 0.505    0.613
#> X2 0.137457  1.147352 0.109477 1.256    0.209
#> X3 0.002192  1.002195 0.119920 0.018    0.985
#> X4 0.064683  1.066820 0.123495 0.524    0.600
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.0613    0.9423     0.8427    1.3365   
#> X2 1.1474    0.8716     0.9258    1.4219   
#> X3 1.0022    0.9978     0.7923    1.2677   
#> X4 1.0668    0.9374     0.8375    1.3590
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
#> X1 0.08670   1.09057  0.04759 1.822 0.068455 .  
#> X2 0.18949   1.20863  0.05004 3.787 0.000153 ***
#> X3 0.19417   1.21430  0.05302 3.662 0.000250 ***
#> X4 0.14162   1.15213  0.04860 2.914 0.003570 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.0906    0.9169     0.9935    1.1972   
#> X2 1.2086    0.8274     1.0957    1.3332   
#> X3 1.2143    0.8235     1.0945    1.3473   
#> X4 1.1521    0.8680     1.0475    1.2673
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
#> X1 0.04232   1.04323  0.05288 0.800 0.423529    
#> X2 0.08623   1.09006  0.05461 1.579 0.114362    
#> X3 0.20676   1.22969  0.05469 3.780 0.000157 ***
#> X4 0.14109   1.15153  0.04930 2.862 0.004209 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.0432    0.9586     0.9405    1.1572   
#> X2 1.0901    0.9174     0.9794    1.2132   
#> X3 1.2297    0.8132     1.1047    1.3688   
#> X4 1.1515    0.8684     1.0455    1.2683
```
