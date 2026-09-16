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
                           covariate = "point_obs_kriged")

cmp[grepl("treat", cmp$term), ]
#>     term baseline_estimate baseline_se integrated_estimate integrated_se
#> 2 treatB         0.6142651   0.1018629            0.676549    0.08377874
#> 3 treatC         1.3400755   0.1123979            1.391838    0.09039698
#>   estimate_change   se_change
#> 2      0.06228387 -0.01808416
#> 3      0.05176207 -0.02200091
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
                      covariate = c("soil_n_kriged", "soil_p_kriged"))
```

## Choosing an engine

| `engine` | Spatial residual | Random effects | Licence |
|----|----|----|----|
| `"ofe"` (default) | `ar1(row):ar1(col)`, via [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md) | yes | open source |
| `"asreml"` | `ar1(row):ar1(col)` | yes | commercial |
| `"lme"` | exponential, via [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | yes | open source |
| `"gls"` | exponential, via [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html) | no | open source |
| `"lm"` | none | no | open source |

`"ofe"` is the default: it fits the same separable AR1 model asreml
does, in base R. `"lm"` ignores spatial correlation entirely and is best
kept as a diagnostic reference rather than a final model.
[`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
returns the same tidy table whichever engine fitted the model, so
downstream code does not change when a licence appears or disappears.

## Step 4 — fit the model yourself

[`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
is a convenience wrapper. When you want a residual structure of your
own, call
[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
directly; it takes the asreml spellings, so the formula below is the one
you would write for `asreml::asreml()`.

``` r

fit <- fit_ofe(yield ~ treat + point_obs_kriged, data = kr,
               residual = ~ ar1(row):ar1(col))
summary(fit)
#> Spatial mixed model fitted by REML (ofeIntegrateR)
#> 
#> Fixed:     yield ~ treat + point_obs_kriged 
#> Residual:  ~ar1(row):ar1(col) 
#> 
#> Observations: 318  Fixed parameters: 4  REML logLik: -151.128 
#> 
#> Variance parameters:
#>           estimate std.error variance
#> R!row!cor   0.4323   0.05410       NA
#> R!col!cor   0.3737   0.05527       NA
#> sigma2      0.2054   0.01639   0.2054
#> 
#> Fixed effects:
#>                  estimate std.error      t      p
#> (Intercept)       0.05705   0.07507  0.760  0.448
#> treatB            0.67655   0.08378  8.075 <1e-04
#> treatC            1.39184   0.09040 15.397 <1e-04
#> point_obs_kriged  0.71165   0.08835  8.055 <1e-04
#> 
#> 6 row(s) dropped for missing values.
```

[`wald_tests()`](https://www.zcao.space/ofeIntegrateR/reference/wald_tests.md)
tests each fixed term as a whole — the question a trial report asks —
rather than one coefficient at a time:

``` r

wald_tests(fit)
#>               term df      wald         F      p.value
#> 1            treat  2 237.09564 118.54782 4.418959e-39
#> 2 point_obs_kriged  1  64.88742  64.88742 1.666269e-14
```

and
[`ofe_means()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_means.md)
gives the treatment means, or the contrasts between them with the
standard error of difference:

``` r

ofe_means(fit, "treat")
#>   treat estimate  std.error        lower     upper
#> 1     A 0.135865 0.07374326 -0.009228365 0.2809584
#> 2     B 0.812414 0.07388445  0.667042833 0.9577852
#> 3     C 1.527703 0.07506330  1.380012017 1.6753932
ofe_means(fit, "treat", pairwise = TRUE)
#>   level1 level2 contrast  estimate  std.error statistic      p.value
#> 1      A      B    B - A 0.6765490 0.08377874  8.075425 1.453454e-14
#> 2      A      C    C - A 1.3918376 0.09039698 15.396948 3.044251e-40
#> 3      B      C    C - B 0.7152886 0.08428020  8.487031 8.519006e-16
```

[`ofe_lsd()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_lsd.md)
puts those together into the table a report prints: means sorted
best-first, with the a/b/c letters. Levels sharing a letter are not
separable at `alpha`.

``` r

tab <- ofe_lsd(fit, "treat")
tab
#>   treat estimate  std.error        lower     upper group
#> 1     C 1.527703 0.07506330  1.380012017 1.6753932     a
#> 2     B 0.812414 0.07388445  0.667042833 0.9577852     b
#> 3     A 0.135865 0.07374326 -0.009228365 0.2809584     c
attr(tab, "lsd")
#>   average_sed       lsd  df alpha adjust      use
#> 1  0.08615197 0.1695081 314  0.05   none pairwise
```

On a balanced design this reproduces `agricolae::LSD.test()` exactly,
and `adjust = "tukey"` reproduces `agricolae::HSD.test()`. The default
here differs from both in one respect that matters for a spatial model:
each pair is judged on its own standard error of difference, because
neighbouring strips are compared more precisely than distant ones. Pass
`use = "lsd"` for the single average LSD a published table usually
means, and quote the value from `attr(tab, "lsd")` with it.

The default `adjust = "none"` gives unprotected comparisons — the
convention these letters usually carry, and not a family-wise error
rate. With more than three or four treatments, say so in the caption or
adjust:

``` r

ofe_lsd(fit, "treat", adjust = "tukey")
#>   treat estimate  std.error        lower     upper group
#> 1     C 1.527703 0.07506330  1.380012017 1.6753932     a
#> 2     B 0.812414 0.07388445  0.667042833 0.9577852     b
#> 3     A 0.135865 0.07374326 -0.009228365 0.2809584     c
```

## Step 5 — pseudo-environments

A strip runs the length of the paddock, and the treatment effect at one
end need not be the effect at the other.
[`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
finds the zones from the data: it removes the treatment signal,
collapses what is left to a profile along the trial, and segments that
profile by dynamic programming.

The guard against inventing zones is the spatial covariance. A smooth
field always looks like it has regions; what matters is whether a region
is wider than the correlation range, because anything narrower is one
realisation of the same correlated surface. The practical range of a
fitted exponential variogram therefore sets the minimum zone width, and
caps how many zones are on offer.

``` r

z <- partition_pseudo_env(kr, response = "yield", along = "row",
                          treat = "treat")
p <- attr(z, "partition")
c(range = p$range, min_width = p$min_zone_width,
  most_zones_possible = p$max_zones_possible, chosen = p$n_zones)
#>               range           min_width most_zones_possible              chosen 
#>            3.959022            3.959022            3.000000            1.000000
p$zones
#>   zone start  end width   n
#> 1    1   0.5 12.5    12 324
```

On this trial it returns a single zone, which is the right answer:
nothing zone-like was simulated. The trial is only twelve rows deep and
the fitted range is about four, so at most three zones were ever on
offer, and BIC preferred none of them. A function that always found
zones would be worse than useless.

[`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md)
then writes the residual structure those zones imply. It gives a zone
`ar1()` only in the dimensions where it has extent: a zone one row deep
would otherwise ask for an unidentifiable parameter and fail.

``` r

r <- adaptive_residual(z)
r
#> ~dsum(~ar1(row):ar1(col) | zone, levels = c("1"))
#> attr(,"geometry")
#>   zone   n n_row n_col            struct
#> 1    1 324    12    27 ar1(row):ar1(col)
#> attr(,"n_degenerate")
#> [1] 0
attr(r, "geometry")
#>   zone   n n_row n_col            struct
#> 1    1 324    12    27 ar1(row):ar1(col)
```

If more than one zone was found, fitting a zone-specific treatment
effect is one call:

``` r

if (nlevels(z$zone) > 1) {
  fit_z <- fit_ofe(yield ~ zone + zone:treat, data = z, residual = r)
  print(wald_tests(fit_z))
  # One lettering per zone: comparisons never cross a zone boundary, and
  # each zone gets its own LSD.
  print(ofe_lsd(fit_z, "treat", by = "zone"))
}
```

### Zoning from covariates instead, when you have them

[`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
works from the response and cuts across the trial, so every treatment
stays in every zone. When the zones should instead come from what is
known before harvest – elevation, an EM38 or gamma survey, a soil test
grid – and may be any shape,
[`partition_paddock()`](https://www.zcao.space/ofeIntegrateR/reference/partition_paddock.md)
clusters those covariates directly. Here is a synthetic elevation and
soil layer over the same grid:

``` r

kr$elevation <- 100 + 4 * exp(-((kr$y_centre - 30)^2) / 400)
kr$clay <- ifelse(kr$x_centre > 120, 20, 32)
set.seed(3)
kr$elevation <- kr$elevation + rnorm(nrow(kr), 0, 0.2)
kr$clay <- kr$clay + rnorm(nrow(kr), 0, 1)

pz <- partition_paddock(kr, covariates = c("elevation", "clay"),
                        treat = "treat")
attr(pz, "partition")$zones
#>   zone  n area patches x_min x_max     y_min     y_max elevation     clay
#> 1    1 98 7938       1 121.5 238.5 52.581045 106.58104  100.2157 19.77029
#> 2    2 91 7371       1   4.5 112.5 52.581045 106.58104  100.2412 32.04948
#> 3    3 65 5265       1   4.5 112.5  7.581045  43.58104  102.7631 32.05240
#> 4    4 70 5670       1 121.5 238.5  7.581045  43.58104  102.7665 20.09846
#>   treatments
#> 1          3
#> 2          3
#> 3          3
#> 4          3
```

Every zone is a rectangle, and `x_min`…`y_max` are its corners, so a
block can be marked out in the paddock. That is what the default
`method` buys: cutting a rectangle all the way across leaves two
rectangles, so no number of cuts can produce a zone that sends a finger
out between its neighbours. `method = "skater"` relaxes that to
contiguous-but-any-shape, which suits a boundary running at an angle;
`method = "kmeans"` drops contiguity altogether and is worth running
once to see why it is not the default. `patches` counts the connected
pieces of each zone – 1 throughout here, more than 1 as soon as k-means
fragments.

The `treatments` column counts the treatment levels present in each
zone. Zones of arbitrary shape will often fail to contain every
treatment, which is exactly when
[`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
is the right tool instead: a zone that does not see all the treatments
cannot support a `zone:treat` term.

Two cautions. Zones found this way are a description of where the
paddock differs, not proof of a boundary — a smooth field with no step
in it is still split more often than not. And because the cuts were
placed where the residuals already differed, a zone-by-treatment
interaction estimated on self-derived zones is exploratory: use the
zones freely for the residual structure, and confirm an interaction
against an independent layer (elevation, EM38, last season’s yield) by
passing that layer as `response` instead.

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
