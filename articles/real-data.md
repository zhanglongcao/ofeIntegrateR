# A worked example on real trial data

The other article uses simulated trials, where the truth is known. This
one uses a real one, where it is not — and where the data arrive in the
state real data arrive in.

[`agridat::lasrosas.corn`](https://kwstat.github.io/agridat/reference/lasrosas.corn.html)
is an on-farm nitrogen experiment from Las Rosas, Córdoba, Argentina,
recorded by a yield monitor. It has what this package is built for:
GPS-referenced yield observations, a categorical treatment applied
across the paddock, and an ancillary spatial layer.

``` r

library(ofeIntegrateR)
library(agridat)

d <- subset(lasrosas.corn, year == 1999)
str(d)
#> 'data.frame':    1738 obs. of  9 variables:
#>  $ year : int  1999 1999 1999 1999 1999 1999 1999 1999 1999 1999 ...
#>  $ lat  : num  -33.1 -33.1 -33.1 -33.1 -33.1 ...
#>  $ long : num  -63.8 -63.8 -63.8 -63.8 -63.8 ...
#>  $ yield: num  72.1 73.8 77.2 76.3 75.5 ...
#>  $ nitro: num  132 132 132 132 132 ...
#>  $ topo : Factor w/ 4 levels "E","HT","LO",..: 4 4 4 4 4 4 4 4 4 4 ...
#>  $ bv   : num  163 170 168 177 171 ...
#>  $ rep  : Factor w/ 3 levels "R1","R2","R3": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ nf   : Factor w/ 6 levels "N0","N1","N2",..: 6 6 6 6 6 6 6 6 6 6 ...
```

`yield` is in quintals per hectare, `nf` is the nitrogen treatment (six
levels), `bv` is a brightness value from aerial imagery, and `topo` is a
topographic class.

## Project the coordinates first

[`grid_dense_layer()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/grid_dense_layer.md)
works in the units of the coordinates it is given, so degrees of
latitude and longitude will produce meaningless cell sizes. Project to
metres before gridding. A local equirectangular approximation is
accurate enough over a single paddock; use a proper projection
([`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html))
for anything larger.

``` r

lat0 <- mean(d$lat)
lon0 <- mean(d$long)
d$x <- (d$long - lon0) * 111320 * cos(lat0 * pi / 180)
d$y <- (d$lat - lat0) * 110540

round(c(width_m = diff(range(d$x)), height_m = diff(range(d$y))))
#>  width_m height_m 
#>      657      240
nrow(d)
#> [1] 1738
```

About 1,700 yield observations over a 660 m × 240 m block.

## Choosing a cell size is a real decision

The cell size trades two things against each other. Large cells average
more yield observations, which suppresses monitor noise. Small cells are
less likely to straddle a treatment boundary.
[`grid_dense_layer()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/grid_dense_layer.md)
reports both sides of that trade-off, so it can be made on evidence
rather than habit.

``` r

sizes <- c(5, 8, 10, 15, 20)
do.call(rbind, lapply(sizes, function(cs) {
  g <- grid_dense_layer(d, x = "x", y = "y", response = "yield", treat = "nf",
                        cell_size = cs)
  data.frame(cell_m = cs,
             cells = nrow(g),
             with_data = sum(g$n_obs > 0),
             median_obs = median(g$n_obs[g$n_obs > 0]),
             pure = round(mean(g$treat_purity == 1, na.rm = TRUE), 2))
}))
#>   cell_m cells with_data median_obs pure
#> 1      5  6336      1738          1 1.00
#> 2      8  2490      1246          1 0.81
#> 3     10  1584       815          2 0.49
#> 4     15   704       385          5 0.11
#> 5     20   396       223          9 0.07
```

The collapse in purity above 8 m is informative: it says the treatment
units in this trial are narrow. They are — all six nitrogen levels
appear across the full extent of the block, so this is a randomised
trial in small plots rather than a strip trial with wide contiguous
strips.

That points to a cell size of 5–8 m. We take 8 m and exclude the cells
that still straddle a boundary.

``` r

g <- grid_dense_layer(d, x = "x", y = "y", response = "yield", treat = "nf",
                      cell_size = 8)

# Carry the brightness layer onto the same lattice
gb <- grid_dense_layer(d, x = "x", y = "y", response = "bv", cell_size = 8)
g$bv <- gb$bv[match(paste(g$row, g$col), paste(gb$row, gb$col))]

# Exclude ambiguous cells, keeping the lattice complete
g$yield[!is.na(g$treat_purity) & g$treat_purity < 0.8] <- NA

c(cells = nrow(g), usable = sum(!is.na(g$yield)))
#>  cells usable 
#>   2490   1010
```

## Standing in for a sparse sampling campaign

This trial has no soil-core layer — `bv` was measured everywhere. To
show the integration workflow we treat it as if it had been *sampled* at
forty locations, the kind of budget a soil-coring campaign would have.
This is a demonstration device: the surface being kriged is real, but
the sparsity is imposed.

``` r

obs <- which(!is.na(g$bv) & !is.na(g$yield))
s <- sample(obs, 40)

cores <- data.frame(x = g$x_centre[s], y = g$y_centre[s], bv_obs = g$bv[s])

g$x <- g$x_centre
g$y <- g$y_centre
kr <- krige_point_samples(cores, g, value = "bv_obs", coords = c("x", "y"))
#> Warning in gstat::fit.variogram(vgm_emp, vgm_start): No convergence after 200
#> iterations: try different initial values?
```

Is forty enough to reconstruct the surface?

``` r

cv <- cv_krige_surface(cores, value = "bv_obs", coords = c("x", "y"))
round(c(rmse = cv$rmse, r2 = cv$r2), 2)
#> rmse   r2 
#> 5.64 0.57
```

A positive `r2` means kriging beats the sample mean, so the layer is
worth carrying forward. Recall that this says the surface is *well
estimated*, not that it is *relevant to yield*.

## Baseline against integrated

``` r

cmp <- compare_integration(kr, response = "yield", treat = "nf",
                           covariate = "bv_obs_kriged", engine = "gls")

cmp[grepl("^nf|bv", cmp$term), ]
#>            term baseline_estimate baseline_se integrated_estimate integrated_se
#> 2          nfN1          3.414609   0.4170008           3.4001387     0.4166863
#> 3          nfN2          5.000892   0.3601712           4.9937017     0.3600925
#> 4          nfN3          6.916330   0.3993302           6.9302149     0.3990885
#> 5          nfN4          9.811794   0.3699268           9.8093736     0.3698157
#> 6          nfN5         10.894283   0.4321121          10.9002948     0.4316057
#> 7 bv_obs_kriged                NA          NA          -0.1745328     0.1308631
#>   estimate_change     se_change
#> 2    -0.014470474 -3.144572e-04
#> 3    -0.007190486 -7.874024e-05
#> 4     0.013885072 -2.417316e-04
#> 5    -0.002420587 -1.110379e-04
#> 6     0.006011384 -5.064099e-04
#> 7              NA            NA
```

Two things to read here.

**The nitrogen response is clean.** The contrasts rise monotonically
from N1 to N5 — a textbook response curve, recovered from irregular
yield-monitor data that started as a GPS point cloud.

**Integration changed almost nothing.** The treatment estimates barely
move and `se_change` is negligible. On this trial, adding a forty-point
brightness surface did not improve the treatment estimates, even though
the surface itself was adequately estimated.

That is not a failure of the method; it is the method working as
intended. The dense yield layer, averaged over several observations per
cell and spread across hundreds of cells, already estimates these
contrasts well. A covariate has to earn its place against that, and here
it did not. Simulation work behind this package found the same pattern
whenever the baseline is strong — which is why
[`compare_integration()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/compare_integration.md)
reports both models rather than the integrated one alone.

## What to take from this

- Project your coordinates before gridding.
- Choose the cell size from `n_obs` and `treat_purity`, not from
  convention. Their trade-off is specific to your trial’s plot geometry.
- Expect integration to help when the baseline is weak, and to do little
  when it is already strong. Report both so the difference is visible.
