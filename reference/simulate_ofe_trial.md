# Simulate a synthetic OFE strip trial with a sparse point-source covariate

Generates a rectangular grid trial with treatment strips, a spatially
correlated "dense" response (e.g. yield), and a spatially correlated
"point-source" covariate (e.g. soil, tissue, or disease measurements)
that is only observed at a sparse set of sampled locations. Intended for
testing and demonstrating
[`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md),
[`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md),
and
[`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md).

## Usage

``` r
simulate_ofe_trial(
  n_row = 40,
  n_col = 20,
  n_treat = 3,
  treat_effects = c(0, 0.8, 1.6),
  point_range = 8,
  point_psill = 1.2,
  point_nugget = 0.1,
  dense_var_weight = 0.6,
  noise_sd = 0.4,
  n_point_samples = 10,
  point_obs_noise_sd = 0.3,
  seed = NULL
)
```

## Arguments

- n_row, n_col:

  Integer grid dimensions (rows = one axis, columns = the axis along
  which treatment strips run).

- n_treat:

  Integer number of treatment strips, applied as contiguous,
  approximately equal-width column blocks.

- treat_effects:

  Numeric vector of length `n_treat` giving the true treatment effect
  added to the dense response. Recycled with treatment labels
  `LETTERS[1:n_treat]`.

- point_range, point_psill, point_nugget:

  Variogram parameters (exponential model) for the spatially correlated
  point-source surface.

- dense_var_weight:

  Numeric; coefficient relating the point-source surface to the dense
  response (i.e. how much of the dense response's spatial pattern is
  explained by the point-source variable).

- noise_sd:

  Standard deviation of i.i.d. measurement noise added to the dense
  response.

- n_point_samples:

  Integer number of sparse point-source samples to draw from the grid.

- point_obs_noise_sd:

  Standard deviation of i.i.d. measurement (e.g. lab) noise added to the
  observed point-source samples.

- seed:

  Optional integer seed for reproducibility.

## Value

A list with components:

- grid:

  Data frame, one row per grid cell, with `row`, `col`, `treat`,
  `point_true` (true point-source surface) and `dense_response` (e.g.
  simulated yield).

- point_samples:

  Data frame of `n_point_samples` sparse observations with `row`, `col`,
  `point_true`, and `point_obs` (true value plus measurement noise).

- true_effects:

  Named numeric vector of the true treatment effects used to simulate
  `dense_response`.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 20, n_col = 10, n_point_samples = 8, seed = 1)
head(sim$grid)
#>   col row treat  point_true dense_response
#> 1   1   1     A -1.77966235   -0.689163155
#> 2   2   1     A -1.02978579   -0.616111989
#> 3   3   1     A -0.26136150   -0.297745820
#> 4   4   1     B  0.03343362    0.608181968
#> 5   5   1     B  0.22460620    1.230599408
#> 6   6   1     B -0.63940352   -0.009025079
sim$point_samples
#>     row col point_true  point_obs
#> 132  14   2 -0.7245681 -0.8285594
#> 11    2   1 -1.8148498 -2.3667263
#> 105  11   5 -1.1524136 -0.8828369
#> 130  13  10 -0.1342264 -0.4980829
#> 10    1  10 -2.9087551 -2.9744444
#> 34    4   4 -0.8357637 -0.6664833
#> 2     1   2 -1.0297858 -1.1874161
#> 124  13   4 -1.3824246 -1.1591123
sim$true_effects
#>   A   B   C 
#> 0.0 0.8 1.6 
```
