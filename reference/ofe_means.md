# Predicted means for a fixed term

The analogue of `predict.asreml()`: the model-based mean of each level
of a factor, averaged over the other fixed effects as they actually
occur in the data (covariates are held at their observed mean). With
`pairwise = TRUE` it returns the differences between levels instead,
each with the standard error of the contrast – which is the number a
trial report quotes.

## Usage

``` r
ofe_means(
  object,
  term,
  by = NULL,
  pairwise = FALSE,
  level = 0.95,
  adjust = "none"
)
```

## Arguments

- object:

  An `ofe_fit` from
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md).

- term:

  Character; the name of a factor in the fixed model.

- by:

  Character or `NULL`; a second factor to compute the means within, one
  set per level of it. Use this for a pseudo-environment analysis –
  `ofe_means(fit, "treat", by = "zone")` gives the treatment means of
  each zone. Means are averaged over the rows belonging to that level
  only, and comparisons are made within a level, never across.

- pairwise:

  Logical; return all pairwise differences rather than the means
  themselves.

- level:

  Numeric; confidence level for the interval (default 0.95).

- adjust:

  Multiplicity adjustment applied to the pairwise p-values: `"none"`
  (the default, i.e. unprotected LSD comparisons), `"tukey"`, `"sidak"`,
  or any method taken by
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Ignored
  when `pairwise = FALSE`. With `by`, the adjustment counts the
  comparisons within a level, not across all of them.

## Value

A data frame of predicted means (the `by` level where given, the `term`
level, `estimate`, `std.error`, `lower`, `upper`) or, when
`pairwise = TRUE`, of differences (`level1`, `level2`, `contrast`,
`estimate`, `std.error`, `statistic`, `p.value`).

## See also

[`ofe_lsd()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_lsd.md),
which adds the a/b/c letters to the means.

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
#>   level1 level2 contrast  estimate std.error statistic      p.value
#> 1      A      B    B - A 0.6950687 0.1559570  4.456796 1.826948e-05
#> 2      A      C    C - A 1.4235858 0.1589796  8.954519 3.976049e-15
#> 3      B      C    C - B 0.7285171 0.1411766  5.160324 9.408498e-07
ofe_means(fit, "treat", pairwise = TRUE, adjust = "tukey")
#>   level1 level2 contrast  estimate std.error statistic      p.value
#> 1      A      B    B - A 0.6950687 0.1559570  4.456796 5.405269e-05
#> 2      A      C    C - A 1.4235858 0.1589796  8.954519 2.575717e-14
#> 3      B      C    C - B 0.7285171 0.1411766  5.160324 2.805874e-06
```
