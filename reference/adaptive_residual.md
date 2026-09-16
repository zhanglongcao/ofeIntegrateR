# Build a residual structure that matches the realised geometry

Writes the `dsum()` residual formula for a pseudo-environment analysis,
one independent section per zone. The reason this needs a function
rather than a literal formula is that AR1 needs at least two levels in a
dimension, and real zones are not guaranteed to have them: a narrow zone
may be one row deep, and asking for `ar1()` there gives an
unidentifiable parameter and a failed fit. Each section therefore gets
`ar1()` only in the dimensions where it actually has extent, and `id()`
elsewhere.

## Usage

``` r
adaptive_residual(data, zone = "zone", row = "row", col = "col")
```

## Arguments

- data:

  Data frame containing the zone and position columns.

- zone:

  Character; the pseudo-environment column, e.g. from
  [`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md).

- row, col:

  Character; the two lattice dimensions.

## Value

A one-sided formula, with the per-zone geometry attached as
`attr(, "geometry")` and the number of zones that had to be demoted to
`id()` as `attr(, "n_degenerate")`.

## Details

The formula is plain text that both
[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
and `asreml::asreml()` accept, so the same analysis script runs with or
without a licence.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 40, n_col = 12, point_range = 6, seed = 3)
z <- partition_pseudo_env(sim$grid, response = "dense_response",
                          along = "row", treat = "treat")
r <- adaptive_residual(z)
r
#> ~dsum(~ar1(row):ar1(col) | zone, levels = c("1", "2"))
#> attr(,"geometry")
#>   zone   n n_row n_col            struct
#> 1    1 240    20    12 ar1(row):ar1(col)
#> 2    2 240    20    12 ar1(row):ar1(col)
#> attr(,"n_degenerate")
#> [1] 0
#> <environment: 0x556ca2e86cd0>
attr(r, "geometry")
#>   zone   n n_row n_col            struct
#> 1    1 240    20    12 ar1(row):ar1(col)
#> 2    2 240    20    12 ar1(row):ar1(col)
```
