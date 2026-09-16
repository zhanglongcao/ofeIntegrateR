# Control parameters for [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)

Control parameters for
[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)

## Usage

``` r
ofe_control(
  maxit = 500L,
  tol = 1e-08,
  trace = FALSE,
  max_n = 2500L,
  start = NULL,
  se = TRUE
)
```

## Arguments

- maxit:

  Integer; maximum iterations for each optimiser stage.

- tol:

  Numeric; convergence tolerance passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- trace:

  Logical; print the REML log-likelihood as it is optimised.

- max_n:

  Integer; refuse to fit more than this many observations. The engine is
  dense: it forms and factorises an `n` by `n` matrix at every
  iteration, so cost grows as `n^3`. Aggregate to a coarser lattice with
  [`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md)
  rather than raising this without thinking.

- start:

  Optional named numeric vector of starting values on the natural scale,
  named as in the `varcomp` table (correlations, ranges, relative
  variances). Useful when a fit has trouble converging.

- se:

  Logical; compute asymptotic standard errors for the variance
  parameters from the observed information. Costs one numerical Hessian.

## Value

A list of control settings.
