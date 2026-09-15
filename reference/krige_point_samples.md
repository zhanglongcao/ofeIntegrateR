# Krige sparse point-source samples onto a target grid

Fits an exponential variogram to sparse point-source observations (e.g.
soil cores, tissue samples) and uses ordinary kriging to predict values
across a target grid (e.g. the full trial extent or a dense covariate
grid). This is the first step of the "kriged-covariate" integration
strategy implemented in
[`fit_integrated_kriged()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/fit_integrated_kriged.md).

## Usage

``` r
krige_point_samples(
  point_data,
  newdata,
  value,
  coords = c("col", "row"),
  vgm_start = NULL,
  nmax = 30
)
```

## Arguments

- point_data:

  Data frame of sparse point-source observations, containing the
  coordinate columns and the value column.

- newdata:

  Data frame giving the target locations to predict onto (e.g. every
  cell of the trial grid), containing the same coordinate columns as
  `point_data`.

- value:

  Character; name of the column in `point_data` holding the observed
  point-source value.

- coords:

  Character vector of length 2 giving the coordinate column names
  present in both `point_data` and `newdata`.

- vgm_start:

  A variogram model object from
  [`gstat::vgm()`](https://r-spatial.github.io/gstat/reference/vgm.html)
  giving starting values for
  [`gstat::fit.variogram()`](https://r-spatial.github.io/gstat/reference/fit.variogram.html).
  Defaults to an exponential model with psill = 1, range =
  `max(spatial extent) / 4`, nugget = 0.1.

- nmax:

  Maximum number of nearest observations used for each kriging
  prediction (passed to
  [`gstat::krige()`](https://r-spatial.github.io/gstat/reference/krige.html)).

## Value

`newdata` with two columns appended: `<value>_kriged` (the kriging
prediction) and `<value>_kriged_var` (the kriging variance). The fitted
variogram model is attached as `attr(., "variogram")`.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 12, seed = 1)
krieged <- krige_point_samples(sim$point_samples, sim$grid, value = "point_obs")
#> Warning: No convergence after 200 iterations: try different initial values?
head(krieged)
#>   col row treat  point_true dense_response point_obs_kriged
#> 1   1   1     A -1.77966235   -0.689163155       -1.3403603
#> 2   2   1     A -1.02978579   -0.616111989       -0.9910913
#> 3   3   1     A -0.26136150   -0.297745820       -1.2990695
#> 4   4   1     B  0.03343362    0.608181968       -1.3127666
#> 5   5   1     B  0.22460620    1.230599408       -1.3138654
#> 6   6   1     B -0.63940352   -0.009025079       -1.3140695
#>   point_obs_kriged_var
#> 1         5.759455e-01
#> 2         2.220446e-16
#> 3         5.811386e-01
#> 4         5.865866e-01
#> 5         5.869316e-01
#> 6         5.869669e-01
attr(krieged, "variogram")
#>   model     psill     range
#> 1   Nug 0.1924131 0.0000000
#> 2   Exp 0.3488773 0.3891814
```
