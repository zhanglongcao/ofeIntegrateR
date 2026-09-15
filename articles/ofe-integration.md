# Integrating point-source data with a yield map

An on-farm experiment usually produces two kinds of data. A **dense
layer** — yield monitor, EM38, NDVI — covers the whole trial at low
marginal cost but measures the outcome rather than its cause. A
**point-source layer** — soil cores, tissue tests, disease counts —
carries direct agronomic meaning but is expensive, so a trial can afford
only tens of observations.

This vignette walks through combining them to estimate a treatment
effect. Every step runs on open-source engines; ASReml-R is used
automatically when licensed but is never required.

``` r

library(ofeIntegrateR)
```

## Before the season: where should the cores go?

The sample count fixes the cost, so the design question is where to put
them.
[`place_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/place_point_samples.md)
implements the common schemes.

``` r

design_grid <- simulate_ofe_trial(n_row = 24, n_col = 12, n_point_samples = 5,
                                  seed = 1)$grid

spread <- function(p) round(mean(dist(cbind(p$col, p$row))), 1)
sapply(c("random", "grid", "stratified", "nested"),
       function(d) spread(place_point_samples(design_grid, n = 16, design = d)))
#>     random       grid stratified     nested 
#>        9.3       11.6        9.7        8.4
```

Higher numbers mean samples are further apart. Two points follow from
the simulation work behind this package:

- **Spread beats clustering when the samples will build a covariate
  surface.** A systematic `"grid"` or a `"stratified"` draw beats simple
  random placement, though the margin is second-order next to the sample
  count.
- **`"nested"` is a reconnaissance design, not an interpolation
  design.** It concentrates pairs at short lags, which is what
  identifies a variogram range. Reusing that geometry for the
  integration sample leaves much of the trial far from any observation
  and produces a worse covariate surface for the same money.

## The dense layer, as a real trial produces it

Yield-monitor data are not a tidy grid. They arrive as a cloud of
GPS-referenced pings along harvester passes, with position error, over a
paddock that is rarely a rectangle, and with the occasional pass
missing.
[`simulate_yield_monitor()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_yield_monitor.md)
generates that shape, so the workflow below is the one you would
actually run.

``` r

sim <- simulate_yield_monitor(n_point_samples = 30, seed = 11)

nrow(sim$cloud)
#> [1] 1893
head(sim$cloud, 3)
#>     x        y treat point_true     yield
#> 1 0.0 4.646274     A  -3.050483 -1.807292
#> 2 1.5 4.206177     A  -3.050483 -1.506388
#> 3 3.0 3.772419     A  -3.050483 -2.125932
sim$true_effects
#>   A   B   C 
#> 0.0 0.8 1.6
```

The contrast we are trying to recover is **C − A = 1.6**.

## Step 1 — put the dense layer on an estimable lattice

Spatial mixed models with a separable AR1 residual need a complete
rectangular lattice.
[`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md)
snaps the cloud onto one, averages within cells, and — importantly —
keeps cells that contain no observations as missing values rather than
dropping them.

``` r

g <- grid_dense_layer(sim$cloud, response = "yield", treat = "treat",
                      cell_size = 9)

dim(g)
#> [1] 324   8
table(empty_cells = g$n_obs == 0)
#> empty_cells
#> FALSE  TRUE 
#>   318     6
summary(g$n_obs[g$n_obs > 0])
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   3.000   6.000   6.000   5.953   6.000   7.000
```

Two columns are worth checking before going further.

`n_obs` is how many pings support each cell. Cells resting on one or two
noisy pings can be blanked with `n_min`. `treat_purity` is the share of
a cell’s observations belonging to its assigned treatment; cells
straddling a strip boundary are ambiguous and are usually excluded.

``` r

table(pure = g$treat_purity == 1, useNA = "ifany")
#> pure
#> FALSE  TRUE  <NA> 
#>     1   317     6

# Exclude boundary cells, but keep the lattice complete
g$yield[!is.na(g$treat_purity) & g$treat_purity < 0.8] <- NA
```

## Step 2 — krige the point samples onto the same grid

A real trial samples wherever the sampler could reach, not at cell
centres.
[`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md)
fits a variogram and predicts onto the grid, so the sparse layer becomes
a covariate defined everywhere.

``` r

g$x <- g$x_centre
g$y <- g$y_centre
kr <- krige_point_samples(sim$point_samples, g, value = "point_obs",
                          coords = c("x", "y"))

attr(kr, "variogram")
#>   model     psill    range
#> 1   Nug 0.2046432  0.00000
#> 2   Exp 1.1487721 20.78104
```

Always look at the fitted variogram. A range far larger than the trial,
or a nugget close to the sill, means the surface is not being estimated.

### Is the point layer dense enough?

[`cv_krige_surface()`](https://www.zcao.space/ofeIntegrateR/reference/cv_krige_surface.md)
leaves each sample out in turn and predicts it from the rest.

``` r

cv <- cv_krige_surface(sim$point_samples, value = "point_obs",
                       coords = c("x", "y"))
c(rmse = round(cv$rmse, 3), r2 = round(cv$r2, 3))
#>  rmse    r2 
#> 0.993 0.232
```

Read this as a **sampling** diagnostic. A negative `r2` means kriging
predicts worse than the sample mean — the layer is too sparse or too
noisy to reconstruct its own surface, and the answer is more samples
next season. It does not by itself condemn the current analysis, and a
healthy `r2` does not guarantee the surface is *relevant to yield*: a
soil property can be mapped perfectly and still explain nothing about
the treatment response.

## Step 3 — compare the integrated fit against the free baseline

The estimate to beat is the dense layer analysed on its own with a
spatial residual. It costs nothing and requires no sampling, so it is
the honest comparator for anything that does.
[`compare_integration()`](https://www.zcao.space/ofeIntegrateR/reference/compare_integration.md)
fits both and lines them up.

``` r

cmp <- compare_integration(kr, response = "yield", treat = "treat",
                           covariate = "point_obs_kriged", engine = "gls")

cmp[grepl("treat", cmp$term), ]
#>     term baseline_estimate baseline_se integrated_estimate integrated_se
#> 2 treatB         0.6639032  0.08448540            0.699747    0.07885479
#> 3 treatC         1.3666055  0.09707601            1.407805    0.08832748
#>   estimate_change    se_change
#> 2      0.03584379 -0.005630609
#> 3      0.04119903 -0.008748533
```

Compare `baseline_estimate` and `integrated_estimate` against the truth
of 1.6, and look at `se_change` — negative means the integrated fit is
more precise. Report both models in your write-up: the difference
between them is exactly what the sampling bought.

The two fitted models are kept for inspection:

``` r

names(attr(cmp, "models"))
#> [1] "baseline"   "integrated"
```

To fit either model on its own, use
[`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
directly; `covariate = NULL` gives the baseline.

## Integrating more than one point variable

Pass several kriged surfaces and they all enter the model:

``` r

fit_integrated_kriged(kr, response = "yield", treat = "treat",
                      covariate = c("soil_n_kriged", "soil_p_kriged"),
                      engine = "gls")
```

## Choosing an engine

| `engine` | Spatial residual | Licence |
|----|----|----|
| `"asreml"` | `ar1(row):ar1(col)` | commercial |
| `"gls"` | exponential, via [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html) | open source |
| `"lm"` | none | open source |

`"gls"` is the open-source workhorse. `"lm"` ignores spatial correlation
entirely and is best kept as a diagnostic reference rather than a final
model.
[`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
returns the same tidy table whichever engine fitted the model, so
downstream code does not change when a licence appears or disappears.

## When integration is not worth it

Simulation work behind this package (GRDC project AAGI-CU-RD-OFE,
Milestone 4) found the benefit is real but conditional. Three findings
are worth carrying into planning:

- **Thirty samples is a floor, not a target.** On idealised lattices
  thirty samples were enough. On realistic trial geometry the same
  thirty bought very little, because averaging many yield pings per cell
  already makes the baseline accurate — the tidy benchmark had
  understated the baseline, not overstated the method.
- **Short-range variation defeats it.** Where the point variable
  decorrelates over a distance comparable to the grid cell, no
  affordable number of samples reconstructs the surface. Check the
  fitted variogram range against your cell size before committing.
- **A significant covariate is not evidence that integration helped.**
  With hundreds of grid cells the kriged covariate clears `|t| > 2`
  almost automatically, including where it does nothing for the
  treatment estimate. Judge the method on how the treatment contrasts
  move, which is what
  [`compare_integration()`](https://www.zcao.space/ofeIntegrateR/reference/compare_integration.md)
  reports.

## The joint bivariate model

[`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md)
models the dense and point layers together with a shared spatial random
effect, estimating their cross-covariance instead of treating the kriged
surface as known. It is the right choice when that covariance is itself
of interest, but it needs an unstructured cross-covariance to inform the
treatment effect at all, and that is hard to estimate from sparse point
data. For routine treatment-effect work the kriged-covariate route is
the dependable default.
