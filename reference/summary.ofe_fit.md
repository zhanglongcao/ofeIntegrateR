# Summarise a fitted `ofe_fit`

Summarise a fitted `ofe_fit`

## Usage

``` r
# S3 method for class 'ofe_fit'
summary(object, ...)
```

## Arguments

- object:

  An `ofe_fit` from
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md).

- ...:

  Ignored.

## Value

An object of class `summary.ofe_fit`, printed as a variance component
table and a table of fixed-effect estimates. Standard errors on the
variance parameters are asymptotic, from the observed information of the
REML log-likelihood, and treat `sigma2` as independent of the
correlation parameters; read them as a guide to identifiability rather
than as exact inference.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 16, n_col = 8, seed = 1)
summary(fit_ofe(dense_response ~ treat, data = sim$grid))
#> Spatial mixed model fitted by REML (ofeIntegrateR)
#> 
#> Fixed:     dense_response ~ treat 
#> Residual:  ~ar1(row):ar1(col) 
#> 
#> Observations: 128  Fixed parameters: 3  REML logLik: -105.648 
#> 
#> Variance parameters:
#>           estimate std.error variance
#> R!row!cor  0.19610   0.09611       NA
#> R!col!cor  0.08647   0.09771       NA
#> sigma2     0.30678   0.03881   0.3068
#> 
#> Fixed effects:
#>             estimate std.error      t      p
#> (Intercept)  0.07492    0.1225 0.6117  0.542
#> treatB       0.69507    0.1560 4.4568 <1e-04
#> treatC       1.42359    0.1590 8.9545 <1e-04
```
