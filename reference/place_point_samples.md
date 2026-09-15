# Choose where to put point samples

Picks `n` locations from a trial grid under a named sampling design. The
number of samples fixes the cost, so this function answers the
complementary question: given that budget, where should the cores go?

## Usage

``` r
place_point_samples(
  grid,
  n,
  design = c("random", "grid", "stratified", "nested"),
  coords = c("col", "row"),
  strata = NULL,
  spacings = NULL
)
```

## Arguments

- grid:

  Data frame of candidate locations, typically the output of
  [`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md)
  or the `grid` element of
  [`simulate_ofe_trial()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_ofe_trial.md).

- n:

  Integer; number of samples to place.

- design:

  One of:

  `"random"`

  :   simple random sample of grid cells.

  `"grid"`

  :   systematic lattice spanning the trial, spaced to respect its
      aspect ratio.

  `"stratified"`

  :   one random cell per stratum, strata formed by crossing `strata`
      (if given) with blocks along the long axis. Spatially balanced and
      executable without a lattice survey.

  `"nested"`

  :   clustered stages at geometrically increasing spacings around a few
      centres, for variogram reconnaissance rather than interpolation.

- coords:

  Character vector of length 2 naming the coordinate columns of `grid`
  (defaults to the `col`/`row` lattice indices this package uses).

- strata:

  Optional character; a column of `grid` (e.g. the treatment) whose
  levels should each be represented under `design = "stratified"`.

- spacings:

  Numeric vector of stage spacings for `design = "nested"`, in the units
  of `coords`. Defaults to `NULL`, which derives three stages from the
  trial extent so the design works at any scale; supply your own when
  the lags of interest are known from a prior variogram.

## Value

The selected rows of `grid`, with the design recorded in
`attr(, "design")`.

## Details

Simulation work for the AAGI-CU-RD-OFE project found that spreading
samples over the trial beats both simple random placement and clustered
placement when the samples are there to build a covariate surface,
though the effect is second-order next to the sample count. Clustered
`"nested"` designs concentrate pairs at short lags, which is what
identifies a variogram range; they are a reconnaissance tool and are a
poor choice for interpolation.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 20, n_col = 12, n_point_samples = 5, seed = 1)

set.seed(1)
pts_random <- place_point_samples(sim$grid, n = 12, design = "random")
pts_grid   <- place_point_samples(sim$grid, n = 12, design = "grid")
pts_strat  <- place_point_samples(sim$grid, n = 12, design = "stratified",
                                  strata = "treat")
nrow(pts_grid)
#> [1] 12
table(pts_strat$treat)
#> 
#> A B C 
#> 4 4 4 
```
