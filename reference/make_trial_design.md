# Lay out an on-farm strip trial

Builds the trial map: treatment plots of a given width and length,
arranged in replicate blocks, snapped to a lattice at the resolution the
dense layer will be gridded to. The result is the plan you take to the
field, and it has the same columns as
[`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md),
so sampling locations can be chosen from it with
[`place_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/place_point_samples.md)
before any yield data exist.

## Usage

``` r
make_trial_design(
  treatments,
  n_rep = 4,
  layout = c("strip", "stack"),
  plot_width = 18,
  plot_length = 200,
  cell_size = 9,
  gap = 0,
  randomise = TRUE,
  origin = c(0, 0),
  seed = NULL
)
```

## Arguments

- treatments:

  Character vector of treatment labels, e.g. `c("N0", "N60", "N120")`.

- n_rep:

  Number of replicate blocks. Each block contains every treatment once.

- layout:

  `"strip"` or `"stack"`, as described above.

- plot_width:

  Width of one plot in metres, normally the working width of the
  applicator or seeder.

- plot_length:

  Length of one plot in metres, along the direction of travel. For
  `"stack"` this is the length of each tier, not of the trial.

- cell_size:

  Lattice resolution in metres, normally the harvester swath, so the
  design and the harvested data share a grid.

- gap:

  Buffer between the two tiers of a `"stack"` layout, in metres. Ignored
  for `"strip"`.

- randomise:

  Randomise treatment order within each replicate block (a randomised
  complete block design). `FALSE` repeats the same order in every block,
  giving a systematic arrangement.

- origin:

  Length-2 numeric giving the south-west corner of the trial in field
  coordinates, for placing the trial within a paddock.

- seed:

  Optional integer seed, so a randomised layout is reproducible.

## Value

A data frame with one row per lattice cell: `row`, `col`, `x`, `y` (cell
centres in metres), `treat`, `rep`, `plot` (plot identifier), and `tier`
for stacked layouts. Cells falling in a buffer carry `NA` treatment. The
design is recorded in `attr(, "design")`.

## Details

Two layouts are supported, following the geometry compared in the GRDC
AAGI-CU-RD-OFE project:

- `"strip"`:

  Every plot runs the full length of the trial, side by side across its
  width. The simplest layout to drive and the usual default.

- `"stack"`:

  The same plots split into two tiers separated by a buffer, halving the
  width and roughly doubling the length. Useful where the paddock is too
  narrow for a single row of strips, and it gives a squarer footprint,
  which changes how much spatial variation the trial spans.

`randomise` is not a formality. Grower strip trials have conventionally
been randomised by analogy with small-plot work, but under strong
spatial correlation a systematic arrangement can estimate treatment
contrasts more precisely, because it spreads each treatment evenly
across the field's spatial gradient. Both are produced here so the
choice can be made deliberately rather than by habit.

## Examples

``` r
# Three nitrogen rates, four replicates, 18 m applicator, 200 m runs
d <- make_trial_design(c("N0", "N60", "N120"), n_rep = 4,
                       plot_width = 18, plot_length = 200, seed = 1)
attr(d, "design")$trial_width
#> [1] 216
table(d$treat)
#> 
#>   N0  N60 N120 
#>  176  176  176 

# Treatment order differs between blocks when randomised
unique(d[, c("rep", "plot", "treat")])[1:6, ]
#>    rep plot treat
#> 1    1    1    N0
#> 3    1    2   N60
#> 5    1    3  N120
#> 7    2    4    N0
#> 9    2    5   N60
#> 11   2    6  N120

# The same plots stacked into two tiers, for a narrower paddock
s <- make_trial_design(c("N0", "N60", "N120"), n_rep = 4, layout = "stack",
                       plot_width = 18, plot_length = 200, gap = 10, seed = 1)
attr(s, "design")$trial_width
#> [1] 108

# Feeds the sampling design directly
set.seed(1)
cores <- place_point_samples(d[!is.na(d$treat), ], n = 24,
                             design = "stratified",
                             coords = c("x", "y"), strata = "treat")
table(cores$treat)
#> 
#>   N0  N60 N120 
#>    8    8    8 
```
