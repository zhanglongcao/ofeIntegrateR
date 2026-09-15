# Aggregate an irregular dense layer onto a regular lattice

Yield-monitor and proximal-sensor data arrive as an irregular cloud of
GPS-referenced observations along harvester or sprayer passes, not as
the rectangular lattice that a separable AR1\\\times\\AR1 residual
structure requires. This function is the first step of the integration
pipeline: it snaps the dense layer onto a regular grid of a chosen cell
size, aggregates the response within each cell, and returns a
**complete** lattice with integer `row` and `col` indices ready for
[`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
or
[`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md).

## Usage

``` r
grid_dense_layer(
  data,
  x = "x",
  y = "y",
  response,
  treat = NULL,
  cell_size = 5,
  fun = mean,
  n_min = 1L,
  keep_empty = TRUE
)
```

## Arguments

- data:

  Data frame of dense-layer observations.

- x, y:

  Character; names of the projected coordinate columns in `data`
  (metres, e.g. UTM/MGA easting and northing). Longitude/latitude in
  degrees will produce meaningless cell sizes – project first.

- response:

  Character; name of the column to aggregate (e.g. yield).

- treat:

  Character or `NULL`; name of a treatment column. When supplied, each
  cell is assigned the treatment of the majority of its observations,
  and `treat_purity` reports the majority share so that cells straddling
  a strip boundary can be identified and removed.

- cell_size:

  Numeric; the length of a (square) grid cell in the units of `x` and
  `y`. Choose it from the trial geometry: small enough to preserve the
  spatial pattern, large enough that most cells hold several
  observations.

- fun:

  Function used to aggregate `response` within a cell (default
  [`mean()`](https://rdrr.io/r/base/mean.html)); it is called with
  `na.rm = TRUE`.

- n_min:

  Integer; cells with fewer than `n_min` observations have their
  response set to `NA`. Use this to discard cells supported by one or
  two noisy yield-monitor pings.

- keep_empty:

  Logical; return the complete lattice including cells with no
  observations (default `TRUE`).

## Value

A data frame, one row per lattice cell, with columns `row`, `col`
(integer indices starting at 1), `x_centre`, `y_centre`, the aggregated
`response`, `n_obs`, and – when `treat` is given – the majority `treat`
and its `treat_purity`. The grid origin and cell size are attached as
`attr(, "grid")`.

## Details

Cells containing no observations are retained with an `NA` response when
`keep_empty = TRUE` (the default). This matters: `asreml`'s
`ar1():ar1()` residual is defined over the full row-by-column lattice,
so gaps must be present as missing values rather than dropped rows.

## Examples

``` r
# An irregular cloud of passes across a 100 x 60 m trial
set.seed(1)
pass_y <- rep(seq(2, 58, by = 4), each = 120)
cloud <- data.frame(
  x = rep(seq(1, 99, length.out = 120), times = length(unique(pass_y))),
  y = pass_y + stats::rnorm(length(pass_y), 0, 0.3)
)
cloud$treat <- cut(cloud$x, breaks = c(0, 33, 66, 100), labels = c("A", "B", "C"))
cloud$yield <- as.numeric(cloud$treat) + stats::rnorm(nrow(cloud), 0, 0.5)

g <- grid_dense_layer(cloud, response = "yield", treat = "treat", cell_size = 5)
head(g)
#>   row col x_centre y_centre     yield n_obs treat treat_purity
#> 1   1   1      3.5  3.83559 1.1536451    13     A            1
#> 2   1   2      8.5  3.83559 1.0342722    12     A            1
#> 3   1   3     13.5  3.83559 1.1570720    12     A            1
#> 4   1   4     18.5  3.83559 1.1602974    11     A            1
#> 5   1   5     23.5  3.83559 0.9833511    11     A            1
#> 6   1   6     28.5  3.83559 1.0039545    12     A            1
table(is.na(g$yield))
#> 
#> FALSE 
#>   240 
```
