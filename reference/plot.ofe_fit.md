# Diagnostic plots for a fitted spatial model

Four panels, and the first is the one that earns its place. A spatial
model is fitted precisely to absorb the field's pattern, so the question
it has to answer is whether any pattern is left: if the residual map
still shows patches, the residual structure has not done its job and the
treatment standard errors are optimistic. The residual variogram is the
same question as a curve – a flat one is what you want, and one still
climbing at short lags says the correlation was not absorbed.

## Usage

``` r
# S3 method for class 'ofe_fit'
plot(x, which = 1:4, ...)
```

## Arguments

- x:

  An `ofe_fit` from
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md).

- which:

  Integer vector selecting panels: 1 residual map, 2 residuals against
  fitted values, 3 normal quantile plot, 4 residual variogram.

- ...:

  Ignored.

## Value

A list with the residual variogram (when drawn), invisibly.

## Details

The residual map uses a diverging scale forced symmetric about zero, so
that over- and under-prediction of the same size get the same weight of
ink.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 20, n_col = 12, seed = 1)
fit <- fit_ofe(dense_response ~ treat, data = sim$grid)
plot(fit)

```
