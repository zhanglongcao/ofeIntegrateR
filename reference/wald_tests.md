# Wald tests of the fixed-effect terms

The analogue of `wald.asreml()`: each term is tested by a Wald statistic
on all of its coefficients jointly, conditional on every other term in
the model. This is the test that answers "does treatment matter", as
opposed to the per-coefficient t-statistics in
[`summary.ofe_fit()`](https://www.zcao.space/ofeIntegrateR/reference/summary.ofe_fit.md).

## Usage

``` r
wald_tests(object, ...)

# S3 method for class 'ofe_fit'
wald_tests(object, ...)
```

## Arguments

- object:

  An `ofe_fit` from
  [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md).

- ...:

  Ignored.

## Value

A data frame with one row per fixed term and columns `term`, `df`,
`wald` (chi-square), `F` (`wald/df`), and `p.value`. The p-value uses an
F reference distribution with `n - p` denominator degrees of freedom,
which is approximate: it does not apply a Kenward-Roger correction, so
it is mildly anti-conservative when the variance parameters are poorly
determined.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 16, n_col = 8, seed = 1)
wald_tests(fit_ofe(dense_response ~ treat, data = sim$grid))
#>    term df     wald        F      p.value
#> 1 treat  2 81.44574 40.72287 2.407166e-14
```
