# Simulated Group Survival Data: Multiple-source Transfer Learning

A dataset containing survival data for 5 groups, where part of the
effects and baseline hazards are heterogeneous. Each group consists of
100 individuals. The survival time \\T\\ for group \\i\\ is generated
according to the following model: \$\$
\lambda^{(i)}(t)=\lambda\_{0,i}(t)\exp\left(x^\top\beta^{(i)}\right),
\$\$ where \\\lambda\_{0,i}(t)\\ represents the baseline hazard for
group \\i\\, \\\beta^{(i)}\\ represents the effects for group \\i\\, and
\\x\\ represents the covariates. The covariates
\\x=(x\_{1},\ldots,x\_{20})\\ are generated from a multivariate normal
distribution with mean 0 and covariance matrix \\\Sigma=I\_{20}\\. The
baseline hazard function is defined as: \$\$
\lambda\_{0,i}(t)=\left\\\begin{array}{ll} t^{2}, & \text{if } i=1,3,5
\\ t, & \text{if } i=2,4. \end{array}\right. \$\$ The effects are
defined as: \$\$ \beta^{(i)}=\left\\\begin{array}{ll}
(0.3,0.3,0.3,0.3,0,\ldots,0), & \text{if } i=1, \\
(0.9,0.9,0.3,0.3,0,\ldots,0), & \text{if } i=2,4, \\
(-0.3,-0.3,0.3,0.3,0,\ldots,0), & \text{if } i=3,5. \end{array}\right.
\$\$ The maximum censoring time is fixed at 2, with an approximate
censoring rate of 20%.

## Format

A data frame with 500 rows and 24 variables:

- id:

  Individual identifier, 1-500

- group:

  Group indicator, 1-5

- time:

  Survival time

- status:

  Status indicator, 0=censored, 1=event

- X1:

  Covariate 1

- X2:

  Covariate 2

- X3:

  Covariate 3

- X4:

  Covariate 4

- X5:

  Covariate 5

- X6:

  Covariate 6

- X7:

  Covariate 7

- X8:

  Covariate 8

- X9:

  Covariate 9

- X10:

  Covariate 10

- X11:

  Covariate 11

- X12:

  Covariate 12

- X13:

  Covariate 13

- X14:

  Covariate 14

- X15:

  Covariate 15

- X16:

  Covariate 16

- X17:

  Covariate 17

- X18:

  Covariate 18

- X19:

  Covariate 19

- X20:

  Covariate 20

## Source

[../articles/simulate-data.html](http://gongziyang.com/survtrans/articles/simulate-data.md)
