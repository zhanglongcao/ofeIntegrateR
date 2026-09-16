# Fit a spatial mixed model without asreml

A residual maximum likelihood (REML) fit of a linear mixed model with an
ASReml-style separable residual structure, written in base R. It exists
so that the analysis at the end of the OFE pipeline – a treatment model
with an `ar1(row):ar1(col)` residual – can be run by a collaborator who
has no ASReml-R licence, using the same formula spelling.

## Usage

``` r
fit_ofe(
  fixed,
  data,
  random = NULL,
  residual = ~ar1(row):ar1(col),
  control = ofe_control()
)
```

## Arguments

- fixed:

  Two-sided formula for the fixed effects, e.g. `yield ~ treat` or
  `yield ~ zone + zone:treat`.

- data:

  Data frame holding every variable used by `fixed`, `random` and
  `residual`. Rows with missing values in any of them are dropped.

- random:

  Optional one-sided formula of random effects, e.g. `~ rep` or
  `~ rep + rep:strip`. Each term contributes one variance component; all
  variables in a term are treated as factors.

- residual:

  One-sided formula giving the residual structure (see above). Defaults
  to `~ ar1(row):ar1(col)`.

- control:

  A list from
  [`ofe_control()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_control.md).

## Value

An object of class `ofe_fit`: a list with elements `coefficients`,
`vcov`, `sigma2`, `varcomp`, `loglik`, `fitted`, `residuals`, `n`, `p`,
`converged`, and the model frame pieces needed by the methods. See
[`summary.ofe_fit()`](https://www.zcao.space/ofeIntegrateR/reference/summary.ofe_fit.md),
[`wald_tests()`](https://www.zcao.space/ofeIntegrateR/reference/wald_tests.md)
and
[`ofe_means()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_means.md).

## Residual structures

`residual` accepts the asreml spellings, combined with `:` for a
separable (Kronecker) structure:

- `id(f)`:

  Independent; a bare variable name means the same thing.

- `ar1(f)`:

  First-order autoregressive across the ordered, equally spaced levels
  of `f`. A factor uses its level order; a numeric column is read as a
  position on a lattice whose spacing is the commonest gap between its
  values, so a row missing from the data still counts as a lag rather
  than being closed up.

- `exp(x)`:

  Exponential correlation in a numeric coordinate `x`, \\\rho =
  \exp(-d/\phi)\\. Use it when positions are irregular.

- `diag(f)`:

  Heterogeneous variance across the levels of `f`, with the first level
  as the reference.

- `dsum(~ struct | s, levels = )`:

  Independent sections: each named level of `s` gets its own copy of
  `struct` with its own parameters and its own variance. This is how a
  pseudo-environment analysis is written;
  [`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md)
  builds the formula for you.

Unlike asreml, cells with a missing response are simply dropped: the
correlation is evaluated from the row and column positions of the
observations that remain, so an incomplete lattice needs no padding.

## What it is not

This is a compact, dense implementation intended for OFE-sized problems
(a few thousand lattice cells). It does not implement the sparse average
information algorithm, `us()`/`fa()` structures, or multi-trait models.
For those, use asreml where licensed;
[`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
will dispatch to it with the same arguments.

## See also

[`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md)
to build a `dsum()` formula from the realised geometry,
[`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
to derive the sections in the first place, and
[`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
for the integration wrapper.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 16, n_col = 8, seed = 1)
fit <- fit_ofe(dense_response ~ treat, data = sim$grid,
               residual = ~ ar1(row):ar1(col))
summary(fit)
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
wald_tests(fit)
#>    term df     wald        F      p.value
#> 1 treat  2 81.44574 40.72287 2.407166e-14
ofe_means(fit, "treat")
#>   treat   estimate std.error      lower     upper
#> 1     A 0.07491886 0.1224767 -0.1674776 0.3173154
#> 2     B 0.76998756 0.1014510  0.5692035 0.9707716
#> 3     C 1.49850465 0.1015156  1.2975927 1.6994166
```
