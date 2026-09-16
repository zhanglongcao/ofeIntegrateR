# Report a fitted variogram

Prints the parameter table, classifies the spatial dependence from the
nugget ratio on the usual convention (Cambardella et al. 1994) and says
what that means for kriging the layer, and warns when the fitted range
falls outside the distances the samples actually cover.

## Usage

``` r
# S3 method for class 'ofe_variogram'
print(x, ...)
```

## Arguments

- x:

  An `ofe_variogram` object from
  [`ofe_variogram()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_variogram.md).

- ...:

  Ignored.

## Value

`x`, invisibly. Called for the printed report.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 30, n_col = 20, n_point_samples = 60,
                          point_range = 6, seed = 1)
pts <- sim$point_samples
pts$x <- pts$col
pts$y <- pts$row
print(ofe_variogram(pts, value = "point_obs"))
#> Variogram of `point_obs`  (60 samples, 14 lag bins)
#> 
#> Fitted model: exponential
#>                     value
#> nugget (c0)        0.3396
#> partial sill (c1)  0.9348
#> total sill         1.2745
#> range parameter    9.7201
#> practical range   29.1603
#> nugget ratio       0.2665
#> 
#>   The fitted practical range (29.2) is longer than the largest lag the
#>   samples cover (10.5). The variogram has not reached its sill inside the
#>   data, so the range is an extrapolation: curves that differ well beyond
#>   the sampled distances fit these points equally. Treat it as a lower
#>   bound, and sample over a longer transect before using it to set a
#>   sampling interval.
#> 
#> Spatial dependence: moderate (nugget ratio 0.27)
#>   Part of the variation is spatially structured. Kriging will smooth, and
#>   the surface will be less variable than the truth.
#> 
#> Models compared (weighted SSE, lower is better):
#>               sse
#> exponential 12.49
#> spherical   12.61
#> gaussian    13.25
```
