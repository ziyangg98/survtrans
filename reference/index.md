# Package index

## Fitting

- [`coxtrans()`](http://gongziyang.com/survtrans/reference/coxtrans.md)
  : Transfer Learning Cox Model with Prior Constraints
- [`coxmtl()`](http://gongziyang.com/survtrans/reference/coxmtl.md) :
  Symmetric Multi-Task Cox Model with Global Centers
- [`cv.coxtrans()`](http://gongziyang.com/survtrans/reference/cv.coxtrans.md)
  [`coef(`*`<cv.coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/cv.coxtrans.md)
  [`predict(`*`<cv.coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/cv.coxtrans.md)
  : Cross-validated tuning for coxtrans

## Inference

- [`logLik(`*`<coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/logLik.coxtrans.md)
  :

  Log-likelihood for a `coxtrans` object

- [`logLik(`*`<coxmtl>`*`)`](http://gongziyang.com/survtrans/reference/logLik.coxmtl.md)
  :

  Log-likelihood for a `coxmtl` object

- [`coef(`*`<coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/coef.coxtrans.md)
  :

  Extract the coefficients from a `coxtrans` object

- [`coef(`*`<coxmtl>`*`)`](http://gongziyang.com/survtrans/reference/coef.coxmtl.md)
  :

  Extract the coefficients from a `coxmtl` object

- [`vcov(`*`<coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/vcov.coxtrans.md)
  :

  Variance-covariance matrix for a `coxtrans` object.

- [`vcov(`*`<coxmtl>`*`)`](http://gongziyang.com/survtrans/reference/vcov.coxmtl.md)
  :

  Variance-covariance matrix for a `coxmtl` object

- [`summary(`*`<coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/summary.coxtrans.md)
  :

  Summary method for a `coxtrans` object

- [`summary(`*`<coxmtl>`*`)`](http://gongziyang.com/survtrans/reference/summary.coxmtl.md)
  :

  Summary method for a `coxmtl` object

- [`print(`*`<coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/print.coxtrans.md)
  :

  Print method for a `coxtrans` object

- [`print(`*`<coxmtl>`*`)`](http://gongziyang.com/survtrans/reference/print.coxmtl.md)
  :

  Print method for a `coxmtl` object

- [`print(`*`<summary.coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/print.summary.coxtrans.md)
  :

  Print method for a `summary.coxtrans` object

- [`print(`*`<summary.coxmtl>`*`)`](http://gongziyang.com/survtrans/reference/print.summary.coxmtl.md)
  :

  Print method for a `summary.coxmtl` object

- [`BIC(`*`<coxmtl>`*`)`](http://gongziyang.com/survtrans/reference/BIC.coxmtl.md)
  :

  Bayesian Information Criterion for `coxmtl` objects

- [`diagnose(`*`<coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/diagnose.coxtrans.md)
  : Diagnose Cox Transfer Model's Optimization Process

- [`refit()`](http://gongziyang.com/survtrans/reference/refit.md) :
  Refit a coxmtl model with hard constraints (oracle estimator)

- [`print(`*`<refit.coxmtl>`*`)`](http://gongziyang.com/survtrans/reference/print.refit.coxmtl.md)
  :

  Print method for a `refit.coxmtl` object

## Prediction

- [`basehaz(`*`<coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/basehaz.coxtrans.md)
  :

  Predict the cumulative baseline hazard function for `coxtrans` objects

- [`basehaz(`*`<coxmtl>`*`)`](http://gongziyang.com/survtrans/reference/basehaz.coxmtl.md)
  :

  Predict the cumulative baseline hazard function for `coxmtl` objects

- [`predict(`*`<coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/predict.coxtrans.md)
  :

  Prediction method for `coxtrans` objects.

- [`predict(`*`<coxmtl>`*`)`](http://gongziyang.com/survtrans/reference/predict.coxmtl.md)
  :

  Prediction method for `coxmtl` objects

- [`plot(`*`<cv.coxtrans>`*`)`](http://gongziyang.com/survtrans/reference/plot.cv.coxtrans.md)
  :

  Plot cross-validation curve for a `cv.coxtrans` object

## Utils

- [`calc_lambda_max()`](http://gongziyang.com/survtrans/reference/calc_lambda_max.md)
  : Calculate the maximum value of the penalty parameter lambda

- [`survtrans_control()`](http://gongziyang.com/survtrans/reference/survtrans_control.md)
  : Ancillary arguments for controlling survtrans fitting

- [`simsurv_tl()`](http://gongziyang.com/survtrans/reference/simsurv_tl.md)
  : Simulate survival data for a multi-source Cox model

- [`ncvcox()`](http://gongziyang.com/survtrans/reference/ncvcox.md) :
  Non-convex penalized Cox proportional hazards model

- [`logLik(`*`<ncvcox>`*`)`](http://gongziyang.com/survtrans/reference/logLik.ncvcox.md)
  :

  Log-likelihood for a `ncvcox` object

- [`coef(`*`<ncvcox>`*`)`](http://gongziyang.com/survtrans/reference/coef.ncvcox.md)
  :

  Extract the coefficients from a `ncvcox` object

- [`vcov(`*`<ncvcox>`*`)`](http://gongziyang.com/survtrans/reference/vcov.ncvcox.md)
  :

  Variance-covariance matrix for a `ncvcox` object.

- [`summary(`*`<ncvcox>`*`)`](http://gongziyang.com/survtrans/reference/summary.ncvcox.md)
  :

  Summary method for a `ncvcox` object

- [`print(`*`<summary.ncvcox>`*`)`](http://gongziyang.com/survtrans/reference/print.summary.ncvcox.md)
  :

  Print method for a `summary.ncvcox` object

- [`predict(`*`<ncvcox>`*`)`](http://gongziyang.com/survtrans/reference/predict.ncvcox.md)
  :

  Prediction method for `ncvcox` objects.

- [`basehaz(`*`<ncvcox>`*`)`](http://gongziyang.com/survtrans/reference/basehaz.ncvcox.md)
  :

  Predict the cumulative baseline hazard function for `ncvcox` objects

- [`diagnose()`](http://gongziyang.com/survtrans/reference/diagnose.md)
  : Generic function for diagnose

## Data Sets

- [`sim2`](http://gongziyang.com/survtrans/reference/sim2.md) :
  Simulated Group Survival Data: Multiple-source Transfer Learning
