# Extract a tidy table of fixed-effect estimates

Works for `ofe_fit`, `lm`, `gls`, `lme`, `asreml`, and `mmer` model
objects, so the same downstream code can summarise treatment contrasts
regardless of which `engine` was used in
[`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
or
[`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md).

## Usage

``` r
extract_fixed_effects(model)
```

## Arguments

- model:

  A fitted model of class `ofe_fit`, `lm`, `gls`, `lme`, `asreml`, or
  `mmer`.

## Value

A data frame with columns `term`, `estimate`, `se`.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 12, seed = 1)
krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
#> Warning: No convergence after 200 iterations: try different initial values?
fit <- fit_integrated_kriged(krieged, response = "dense_response",
                              treat = "treat", covariate = "point_obs_kriged",
                              row = "row", col = "col")
extract_fixed_effects(fit)
#>               term    estimate        se
#> 1      (Intercept) -0.09683842 0.2851930
#> 2           treatB  1.01542575 0.1486771
#> 3           treatC  1.85459880 0.1453960
#> 4 point_obs_kriged  0.53862942 0.2013503
```
