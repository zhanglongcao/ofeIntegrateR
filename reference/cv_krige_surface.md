# Cross-validate a kriged point-source surface

Leave-one-out cross-validation of the point-source layer: each sample is
dropped in turn, predicted from the rest, and compared with its observed
value. It answers whether the samples are dense enough, and well enough
behaved, to reconstruct their own surface.

## Usage

``` r
cv_krige_surface(
  point_data,
  value,
  coords = c("col", "row"),
  vgm_model = NULL,
  nmax = 30
)
```

## Arguments

- point_data:

  Data frame of point-source observations.

- value:

  Character; name of the observed value column.

- coords:

  Character vector of length 2 naming the coordinate columns.

- vgm_model:

  Optional variogram model from
  [`gstat::vgm()`](https://r-spatial.github.io/gstat/reference/vgm.html).
  By default the same exponential model
  [`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md)
  would fit is used, so the diagnostic matches the surface actually
  entering the analysis.

- nmax:

  Maximum neighbours per prediction, passed to
  [`gstat::krige.cv()`](https://r-spatial.github.io/gstat/reference/krige.cv.html).

## Value

A list with `rmse` (root mean squared prediction error), `r2`
(proportion of variance explained, which is negative when kriging
predicts worse than the sample mean), `n` (number of points) and
`residuals`. The variogram used is attached as `attr(, "variogram")`.

## Details

Treat this as a **sampling** diagnostic, not a verdict on an analysis.
In the simulation work behind this package a better-validated surface
delivered a larger average benefit from integration, but the probability
of benefiting on any individual trial barely moved. A low or negative
`r2` is a reason to collect more samples next season; it does not by
itself condemn the current model. Note also that cross-validation
measures whether the surface is *well estimated*, not whether it is
*relevant to yield* — a soil property can be mapped perfectly and still
explain nothing about the treatment response.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 25, seed = 4)
cv <- cv_krige_surface(sim$point_samples, value = "point_obs")
cv$rmse
#> [1] 0.8526824
cv$r2
#> [1] 0.5337104
```
