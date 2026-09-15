# Build a candidate grid for a planned trial

Lays out a rectangular trial as a lattice of cells with treatment
strips, so that sampling locations can be chosen with
[`place_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/place_point_samples.md)
before any yield data exist. The result has the same shape as the output
of
[`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md)
— `row`, `col`, `x`, `y` and a treatment column — so the same code works
before and after harvest.

## Usage

``` r
make_trial_grid(
  width,
  height,
  cell_size = 9,
  treatments = c("A", "B", "C"),
  strip_width = NULL,
  along = c("x", "y")
)
```

## Arguments

- width, height:

  Trial extent in metres. `width` runs along `x`.

- cell_size:

  Cell size in metres. A natural choice is the harvester swath, since
  that is the resolution the dense layer will be gridded to.

- treatments:

  Character vector of treatment labels.

- strip_width:

  Width of each treatment strip in metres. Defaults to the extent across
  the strips divided by the number of treatments, i.e. one strip per
  treatment. Give a smaller value for replicated strips, which then
  cycle through the treatments.

- along:

  Axis the strips run along: `"x"` (strips stacked up the `y` axis, the
  usual layout for machinery passes) or `"y"`.

## Value

A data frame with one row per cell: `row`, `col`, `x`, `y` (cell centres
in metres) and `treat`. The layout is recorded in `attr(, "layout")`.

## Details

This lays out **strips on a rectangle**. It does not choose the number
of treatments, the replication, or the strip width; those come from the
trial design, which this package does not perform.

## Examples

``` r
# A 240 x 108 m trial, three treatments, replicated strips two swaths wide
g <- make_trial_grid(width = 240, height = 108, cell_size = 9,
                     treatments = c("A", "B", "C"), strip_width = 18)
dim(g)
#> [1] 312   5
table(g$treat)
#> 
#>   A   B   C 
#> 104 104 104 

# Choose 30 coring locations spread across the treatments
set.seed(1)
cores <- place_point_samples(g, n = 30, design = "stratified",
                             coords = c("x", "y"), strata = "treat")
table(cores$treat)
#> 
#>  A  B  C 
#> 10 10 10 
```
