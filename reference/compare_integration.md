# Compare an integrated analysis against the free baseline

Fits the treatment model twice on the same data — once using the dense
layer alone, once with the kriged point-source covariate added — and
reports the treatment contrasts side by side.

## Usage

``` r
compare_integration(
  data,
  response,
  treat,
  covariate,
  row = "row",
  col = "col",
  engine = c("ofe", "asreml", "gls", "lm"),
  ...
)
```

## Arguments

- data:

  Data frame containing the gridded dense response, the treatment factor
  and the kriged covariate(s) — typically the output of
  [`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md).

- response, treat:

  Character; names of the response and treatment columns.

- covariate:

  Character vector; the kriged covariate column(s) to add in the
  integrated fit.

- row, col:

  Character; position columns used by the spatial residual.

- engine:

  Fitting engine, passed to
  [`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md).
  Defaults to `"ofe"`, which needs no licence.

- ...:

  Further arguments passed to
  [`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md).

## Value

A data frame with one row per fixed-effect term, giving the estimate and
standard error under each model and the change between them. The two
fitted models are attached as `attr(, "models")`.

## Details

This is the comparison the method should always be judged on. The dense
layer is already collected, so a spatial model of it costs nothing;
point sampling does. Reporting only the integrated fit hides how much of
its accuracy came from the covariate and how much was available for
free. Report both.

Note that a large, highly significant covariate coefficient is **not**
evidence that integration helped. With hundreds of grid cells the
covariate clears conventional significance almost automatically,
including in situations where it does nothing for the treatment
estimate. Judge the method by how the treatment contrasts and their
standard errors move.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 20, seed = 8)
kr <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
#> Warning: No convergence after 200 iterations: try different initial values?
compare_integration(kr, response = "dense_response", treat = "treat",
                    covariate = "point_obs_kriged", engine = "lm")
#>               term baseline_estimate baseline_se integrated_estimate
#> 1      (Intercept)         0.3885048  0.07114438          -0.1409987
#> 2           treatB         0.7053552  0.10061335           0.7353865
#> 3           treatC         1.5941030  0.10061335           1.4635440
#> 4 point_obs_kriged                NA          NA           0.9137727
#>   integrated_se estimate_change    se_change
#> 1    0.09315540      -0.5295035  0.022011021
#> 2    0.09003254       0.0300313 -0.010580811
#> 3    0.09150187      -0.1305589 -0.009111483
#> 4    0.11745589              NA           NA
```
