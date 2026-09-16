# Predicted means for a fixed term

The analogue of `predict.asreml()`: the model-based mean of each level
of a factor, averaged over the other fixed effects as they actually
occur in the data (covariates are held at their observed mean). With
`pairwise = TRUE` it returns the differences between levels instead,
each with the standard error of the contrast – which is the number a
trial report quotes.

## Usage

``` r
ofe_means(object, term, pairwise = FALSE, level = 0.95)
```

## Arguments

- object:

  An `ofe_fit` from
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md).

- term:

  Character; the name of a factor in the fixed model.

- pairwise:

  Logical; return all pairwise differences rather than the means
  themselves.

- level:

  Numeric; confidence level for the interval (default 0.95).

## Value

A data frame of predicted means (`term` level, `estimate`, `std.error`,
`lower`, `upper`) or, when `pairwise = TRUE`, of differences
(`contrast`, `estimate`, `std.error`, `statistic`, `p.value`).

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 16, n_col = 8, seed = 1)
fit <- fit_ofe(dense_response ~ treat, data = sim$grid)
ofe_means(fit, "treat")
#>   treat   estimate std.error      lower     upper
#> 1     A 0.07491886 0.1224767 -0.1674776 0.3173154
#> 2     B 0.76998756 0.1014510  0.5692035 0.9707716
#> 3     C 1.49850465 0.1015156  1.2975927 1.6994166
ofe_means(fit, "treat", pairwise = TRUE)
#>   contrast  estimate std.error statistic      p.value
#> 1    B - A 0.6950687 0.1559570  4.456796 1.826948e-05
#> 2    C - A 1.4235858 0.1589796  8.954519 3.976049e-15
#> 3    C - B 0.7285171 0.1411766  5.160324 9.408498e-07
```
