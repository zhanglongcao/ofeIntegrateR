# Derive pseudo-environments from the spatial pattern of the dense layer

Splits a trial into contiguous pseudo-environments (PEs) along its
length, so that a treatment effect can be estimated separately within
each. The cut points come from the data rather than from a guess: the
treatment signal is removed, the residual field is collapsed to a
profile along the chosen axis, and the profile is segmented optimally by
dynamic programming.

## Usage

``` r
partition_pseudo_env(
  data,
  response,
  along = "row",
  treat = NULL,
  block = NULL,
  method = c("segment", "equal"),
  n_zones = NULL,
  max_zones = 8L,
  min_zone_width = NULL,
  name = "zone"
)
```

## Arguments

- data:

  Data frame, one row per lattice cell (e.g. from
  [`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md)).

- response:

  Character; the column the zones are derived from. Passing an
  independent layer (elevation, EM38, last season's yield) rather than
  this season's yield avoids the selection problem described above.

- along:

  Character; the coordinate the trial is cut across – normally the
  direction the strips run, so that every treatment appears in every
  zone. Typically `"row"` for strips running up the paddock.

- treat:

  Character or `NULL`; a treatment column to remove before zoning, so
  that the treatment pattern is not mistaken for a spatial one.

- block:

  Character or `NULL`; a further nuisance factor to remove (e.g.
  `"rep"`).

- method:

  `"segment"` (default) for the data-driven segmentation described
  above, or `"equal"` for equal-width slices, which is the convention
  this is meant to be compared against.

- n_zones:

  Integer or `NULL`; force a number of zones instead of choosing it by
  BIC. Required for `method = "equal"`.

- max_zones:

  Integer; upper limit on the number of zones considered.

- min_zone_width:

  Numeric or `NULL`; minimum zone extent in the units of `along`.
  Defaults to the fitted practical range.

- name:

  Character; name of the zone column added to `data`.

## Value

`data` with an added factor column (named by `name`) giving the
pseudo-environment of each cell. Details of the derivation are attached
as `attr(, "partition")`: the fitted `range`, the `min_zone_width` used,
the `breaks`, the per-zone summary, the BIC table, and the profile
itself.

## Details

The role of the spatial covariance is to stop the procedure inventing
zones. A smooth field will always look like it has "regions"; the
question is whether a region is wider than the correlation range,
because anything narrower is one realisation of the same correlated
surface rather than a distinct environment. So the practical range of an
exponential variogram fitted to the profile becomes the minimum zone
width, and it also caps the number of zones at
`floor(trial length / range)`. Within that cap the number of zones is
chosen by BIC.

Two things are worth knowing before trusting the output.

First, a smooth field will be split even when nothing discrete is there.
A paddock whose yield varies continuously genuinely does have a better
end and a worse end, and the procedure will say so. On simulated fields
with no step at all and a practical range of about a third of the trial,
the default settings return one zone about a third of the time and two
or three the rest; with a real step of two-and-a-half field standard
deviations they find it, at the right place, essentially always. Read a
zone as “this part of the paddock behaves differently”, not as evidence
of a boundary.

Second, PEs derived from the same yield data that are then used to
estimate zone-specific treatment effects will overstate those
differences, because the boundaries were placed where the residuals
already differed. Use them for the residual structure (via
[`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md))
with a clear conscience; treat a zone-by-treatment interaction fitted on
self-derived zones as exploratory, and confirm it against an independent
layer – elevation, EM38, a prior season's yield – by passing that layer
as `response` instead.

## See also

[`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md),
which turns the zones into a `dsum()` residual formula, and
[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md),
which fits it.

## Examples

``` r
sim <- simulate_ofe_trial(n_row = 40, n_col = 12, point_range = 6, seed = 3)
z <- partition_pseudo_env(sim$grid, response = "dense_response",
                          along = "row", treat = "treat")
table(z$zone)
#> 
#>   1   2 
#> 240 240 
attr(z, "partition")$range
#> [1] 20
```
