# Partition a paddock into contiguous zones from environmental covariates

Divides a paddock into pseudo-environments using covariates that are
known before harvest – elevation, EM38 or gamma survey, soil test grids,
a previous season's imagery – rather than from the yield being analysed.
The zones are guaranteed to be spatially contiguous.

## Usage

``` r
partition_paddock(
  data,
  covariates,
  k = NULL,
  row = NULL,
  col = NULL,
  x = NULL,
  y = NULL,
  method = c("rectangle", "skater", "kmeans"),
  neighbours = c("auto", "rook", "queen", "knn"),
  n_neighbours = 8L,
  k_max = 8L,
  criterion = c("gain", "ch"),
  min_gain = 0.05,
  min_cells = NULL,
  min_r2 = 0.1,
  scale = TRUE,
  n_pc = NULL,
  treat = NULL,
  name = "zone"
)
```

## Arguments

- data:

  Data frame, one row per cell (e.g. from
  [`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md)).

- covariates:

  Character vector naming the columns of `data` to zone on, such as
  `c("elevation", "clay")`. Required, and there is deliberately no
  default: which layers define a zone is an agronomic decision. All must
  be numeric; rows missing any of them are dropped.

- k:

  Integer or `NULL`; the number of zones. `NULL` chooses it (see above).

- row, col:

  Character or `NULL`; integer lattice indices, used to build the grid
  neighbourhoods. Defaults to `"row"`/`"col"` when `data` has them.

- x, y:

  Character or `NULL`; projected coordinates in metres. Defaults to
  `"x_centre"`/`"y_centre"` when present, else `"x"`/`"y"`. Needed for
  the spatial range and for `neighbours = "knn"`.

- method:

  `"rectangle"` (default) for rectangular blocks, `"skater"` for
  contiguous zones of any shape, or `"kmeans"` for unconstrained
  clustering (not contiguous; for comparison). See above.

- neighbours:

  `"auto"` (grid neighbourhoods when `row`/`col` are available,
  otherwise nearest neighbours), `"rook"`, `"queen"` or `"knn"`.

- n_neighbours:

  Integer; neighbours per cell when `neighbours = "knn"`.

- k_max:

  Integer; largest number of zones considered when `k` is chosen
  automatically.

- criterion:

  `"gain"` (default) to stop when another zone would explain less than
  `min_gain` more of the covariate variance, or `"ch"` to maximise the
  Calinski-Harabasz index. See above.

- min_gain:

  Numeric; the share of covariate variance an extra zone must explain to
  be kept, under `criterion = "gain"`.

- min_cells:

  Integer or `NULL`; smallest zone allowed. Defaults to 2% of the cells,
  with a floor of 5.

- min_r2:

  Numeric; if the chosen partition explains less than this share of
  covariate variance, one zone is returned instead.

- scale:

  Logical; standardise each covariate to mean 0 and standard deviation 1
  before clustering (default `TRUE`). Leave it on unless the covariates
  are already on a common scale – otherwise elevation in metres will
  silently outrank pH.

- n_pc:

  Integer or `NULL`; cluster on the first `n_pc` principal components
  instead of the covariates themselves. Useful when several layers
  measure much the same thing.

- treat:

  Character or `NULL`; a treatment column, used only to report and warn
  about treatment coverage within zones.

- name:

  Character; name of the zone column added to `data`.

## Value

`data` with an added factor column (named by `name`) giving the zone of
each cell; cells dropped for missing covariates get `NA`. The derivation
is attached as `attr(, "partition")`: the `method`, `k`, `covariates`,
the fitted `range`, the `table` of criteria over k (`ssd`, `r2`, the
marginal `gain`, and `ch`), and a `zones` summary with each zone's cell
count, area, number of connected `patches`, corner coordinates,
covariate means and (with `treat`) treatment coverage.

## Shape of the zones

Three methods, differing only in what shape a zone is allowed to be. The
choice is a practical one, not a statistical one: a zone has to be
something a machine can drive and a sampling plan can stratify by.

- `"rectangle"` (default):

  Rectangular blocks. The paddock is cut repeatedly by lines running the
  full width or full length of the region being split, each cut placed
  where it removes most within-zone variance. Cutting a rectangle across
  always leaves two rectangles, so every zone is a rectangle however
  many cuts are made. This is the shape a variable-rate prescription, a
  soil sampling grid or a set of management blocks wants, and it is the
  two-dimensional version of what
  [`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
  does along a single axis. Where the paddock is not itself a rectangle
  – a clipped corner, a headland left out – a zone is a rectangle
  intersected with the paddock, which is still a block you can drive.

- `"skater"`:

  Contiguous zones of any shape. A minimum spanning tree is built over
  the neighbourhood graph – edges weighted by distance between cells in
  standardised covariate space – and pruned one edge at a time, always
  cutting the edge that removes most within-zone variance. Every zone is
  a subtree of a connected graph, so every zone is connected, but a zone
  may send a finger out between two others. Use it when the zones are
  there to describe the soil faithfully rather than to be worked, and
  when a boundary genuinely runs at an angle or along a contour. Being
  free to follow the boundary does not guarantee it fits better: both
  methods are greedy, and a greedy tree pruning can strand itself where
  a well-placed straight cut would not, so rectangles sometimes explain
  more of the covariate variance at the same `k`. Compare the `r2`
  column of the returned table on your own paddock rather than assuming.
  The method is SKATER (Assuncao et al. 2006).

- `"kmeans"`:

  Plain k-means on the covariates, for comparison. It is not contiguous
  at all: cells are assigned to the zone their soil resembles wherever
  they happen to sit, so zones come back as confetti scattered across
  the paddock. The `patches` column of the zone summary counts the
  connected pieces of each zone – 1 for every rectangle and every skater
  zone, more than 1 whenever k-means fragments.

What the choice costs depends entirely on the shape of the underlying
feature. On simulated paddocks with one clear feature in one covariate,
the variance explained at `k = 2` was:

|                                  |               |            |
|----------------------------------|---------------|------------|
| **feature**                      | **rectangle** | **skater** |
| block boundaries on the axes     | 0.78          | 0.78       |
| a smooth gradient up the paddock | 0.69          | 0.67       |
| a boundary running diagonally    | 0.25          | 0.98       |
| a round patch in the middle      | 0.13          | 0.97       |

So rectangles cost nothing when the structure is already blocky or a
gradient, and cost almost everything when the boundary runs at an angle
or curves – a creek line, a dune, a contour. If your covariate map has a
feature like that, either use `"skater"` and accept the awkward shapes,
or raise `k` so rectangles can approximate the boundary in steps.

## Choosing k

With `k = NULL`, zones are added while each new one earns its keep: the
default `criterion = "gain"` stops at the first `k` whose successor
would explain less than `min_gain` (5% by default) more of the covariate
variance. This is the elbow rule made explicit, and it is the one an
agronomist can argue with, because zones cost management effort and the
threshold is where that trade-off is stated.

`criterion = "ch"` maximises the Calinski-Harabasz index instead. Be
aware that on a smoothly varying covariate – a slope, a gradual texture
change – CH often has no interior maximum at all and simply picks
`k_max`, because every further split of a smooth surface really does
separate it further. That is a property of smooth fields rather than a
failure of the index, and it is why it is not the default here.

Either way, if the partition explains less than `min_r2` of the
covariate variance, one zone is returned: a paddock uniform in its
covariates should be reported as uniform. The full criterion table is
returned so `k` can be overridden knowingly.

The practical range of an exponential variogram fitted to the covariate
field is reported alongside, because it tells you the scale at which the
covariates vary and hence how far apart to sample. It deliberately does
**not** constrain `k`: a real zone boundary is a step, a step keeps the
variogram climbing without reaching a sill, and the fitted range is then
inflated by the very structure being looked for. Read it as context, not
as a limit.

## Using the zones in a trial analysis

Pass `treat` and the zone summary gains a `treatments` column counting
the treatment levels present in each zone. A zone missing a treatment
cannot support a `zone:treat` term, and you are warned rather than left
to discover it when the fit drops rank. Zones from this function are
arbitrary shapes, so a strip trial will often want
[`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
instead, which cuts across the strips and so keeps every treatment in
every zone.

## References

Assuncao, R.M., Neves, M.C., Camara, G. and da Costa Freitas, C. (2006)
Efficient regionalization techniques for socio-economic geographical
units using minimum spanning trees. *International Journal of
Geographical Information Science* **20**, 797-811.

## See also

[`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
for slicing a strip trial across its length, and
[`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md)
to turn either set of zones into a residual structure.

## Examples

``` r
# A paddock with a ridge running through it and a sandy corner
g <- expand.grid(row = 1:24, col = 1:18)
g$x <- g$col * 10
g$y <- g$row * 10
g$elevation <- 100 + 6 * exp(-((g$y - 120)^2) / 2000)
g$clay <- ifelse(g$x > 120 & g$y > 140, 18, 32)
set.seed(1)
g$elevation <- g$elevation + stats::rnorm(nrow(g), 0, 0.2)
g$clay <- g$clay + stats::rnorm(nrow(g), 0, 1)

z <- partition_paddock(g, covariates = c("elevation", "clay"))
attr(z, "partition")$zones
#>   zone   n  area patches x_min x_max y_min y_max elevation     clay
#> 1    1  48  4800       1   130   180   170   240  100.4655 18.11420
#> 2    2  96  9600       1    10   120   170   240  100.4691 31.81901
#> 3    3 126 12600       1    10   180    10    70  100.4996 31.92319
#> 4    4 162 16200       1    10   180    80   160  104.4964 30.91278
table(z$zone)
#> 
#>   1   2   3   4 
#>  48  96 126 162 
```
