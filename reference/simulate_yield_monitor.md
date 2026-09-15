# Simulate a trial as a yield monitor actually records it

[`simulate_ofe_trial()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_ofe_trial.md)
produces a tidy lattice: one observation per cell, treatments aligned to
whole columns. Real dense layers look nothing like that, and code that
only ever sees the tidy version tends to break on first contact with a
harvester file. This function generates the awkward version —
GPS-referenced points along machinery passes, position error, a paddock
that is not a rectangle, occasional missing passes, and point samples
taken wherever the sampler could reach — so a workflow can be tested end
to end through
[`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md).

## Usage

``` r
simulate_yield_monitor(
  field_x = 240,
  field_y = 108,
  swath = 9,
  along_spacing = 1.5,
  band_swaths = 2L,
  treat_effects = c(0, 0.8, 1.6),
  point_range = 25,
  point_psill = 1.2,
  point_nugget = 0.1,
  dense_var_weight = 0.6,
  noise_sd = 0.4,
  gps_jitter = 0.5,
  skip_pass_prob = 0.08,
  cut_corner = TRUE,
  n_point_samples = 30L,
  point_obs_noise_sd = 0.3,
  point_design = c("random", "strip_ends"),
  truth_res = 3,
  seed = NULL
)
```

## Arguments

- field_x, field_y:

  Paddock extent, in metres.

- swath:

  Harvester swath width (m); passes run along the `x` axis.

- along_spacing:

  Distance between recorded points within a pass (m).

- band_swaths:

  Number of adjacent swaths per treatment strip.

- treat_effects:

  Named or unnamed numeric vector of true treatment effects; treatments
  are labelled `A`, `B`, ... in order.

- point_range, point_psill, point_nugget:

  Exponential variogram parameters for the point-source surface, in
  metres.

- dense_var_weight:

  Coefficient relating the point-source surface to the dense response.

- noise_sd:

  Standard deviation of measurement noise on the dense response.

- gps_jitter:

  Standard deviation of across-track position error (m).

- skip_pass_prob:

  Probability that a pass is missing entirely.

- cut_corner:

  Clip a triangle off one corner, so the paddock is not a rectangle.

- n_point_samples:

  Number of point-source samples to draw.

- point_obs_noise_sd:

  Measurement noise on the point samples.

- point_design:

  `"random"` spreads samples over the trial; `"strip_ends"` places them
  at both ends of every strip, the pattern growers commonly use when
  each sample is expensive to collect.

- truth_res:

  Resolution (m) of the grid on which the point-source surface is
  simulated before being read off at the recorded locations.

- seed:

  Optional integer seed.

## Value

A list with `cloud` (the irregular dense observations, with `x`, `y`,
`treat`, `point_true` and `yield`), `point_samples` (sparse observations
at arbitrary coordinates) and `true_effects`.

## Details

Treatment strips run along the direction of travel, as they do in
practice, each covering `band_swaths` adjacent passes.

## Examples

``` r
sim <- simulate_yield_monitor(n_point_samples = 20, seed = 2)
nrow(sim$cloud)
#> [1] 1733

# The whole point: this needs gridding before it can be modelled
g <- grid_dense_layer(sim$cloud, response = "yield", treat = "treat",
                      cell_size = 9)
table(empty = g$n_obs == 0)
#> empty
#> FALSE  TRUE 
#>   291    33 
```
