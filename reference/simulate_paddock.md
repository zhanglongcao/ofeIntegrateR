# Simulate a paddock, its covariate layers, and a trial in it

Generates a paddock whose truth is known, so that the rest of the
package can be checked rather than merely run: the zones it should find,
the treatment effect it should recover, and how much of the yield was
never attributable to treatment in the first place.

## Usage

``` r
simulate_paddock(
  n_row = 40L,
  n_col = 24L,
  cell_size = 10,
  yield_mean = 3,
  yield_range = 80,
  yield_sd = 0.6,
  covariates = NULL,
  covariate_cor = 0.6,
  n_zones = 1L,
  zone_props = NULL,
  treatments = NULL,
  n_rep = 4L,
  treat_effect = NULL,
  zone_effect = NULL,
  noise_sd = 0.3,
  randomise = TRUE,
  seed = NULL
)
```

## Arguments

- n_row, n_col:

  Lattice size of the paddock.

- cell_size:

  Cell size in metres.

- yield_mean:

  Mean yield potential, in whatever unit you are working in (t/ha, say).

- yield_range:

  Practical range of the yield-potential surface, in metres. Short
  relative to the paddock gives a patchy field; long gives a gradient.

- yield_sd:

  Standard deviation of the yield-potential surface.

- covariates:

  Named list of covariate layers to draw. Each element is a list which
  may name `mean`, `sd`, `range` and `cor` (its correlation with yield
  potential). Defaults supply `elevation` and `clay`. Use `NULL` for
  none.

- covariate_cor:

  Default correlation between a covariate and yield potential, for
  layers that do not set their own `cor`.

- n_zones:

  Integer; number of pseudo-environments across the paddock. `1` (the
  default) means no zone structure.

- zone_props:

  Optional proportions of the paddock length given to each zone;
  defaults to equal bands. Real boundaries are not tidy, so unequal
  proportions are the more realistic test.

- treatments:

  Optional treatment labels. When given, a trial is laid out and
  harvested; when `NULL`, only the paddock and its layers are returned.

- n_rep:

  Replicates, passed to
  [`make_trial_design()`](https://www.zcao.space/ofeIntegrateR/reference/make_trial_design.md).

- treat_effect:

  Numeric vector of treatment effects, one per treatment, in yield
  units. Defaults to an evenly spaced response.

- zone_effect:

  Numeric vector of multipliers, one per zone, scaling the treatment
  effect in that zone. `rep(1, n_zones)` means the treatment works
  equally everywhere.

- noise_sd:

  Standard deviation of the independent harvest noise.

- randomise:

  Passed to
  [`make_trial_design()`](https://www.zcao.space/ofeIntegrateR/reference/make_trial_design.md).

- seed:

  Optional random seed.

## Value

A data frame with one row per cell: `row`, `col`, `x`, `y`,
`yield_potential`, one column per covariate, `zone`, and – when
`treatments` was given – `treat`, `rep`, `plot` and `yield`. Classed as
`ofe_zones` so [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
maps the zones. `attr(, "truth")` holds the settings and the realised
component variances.

## What is simulated

A spatially correlated **yield potential** surface is drawn first; it is
the paddock the grower already has. Each named **covariate** is then
drawn as a mixture of that surface and its own independent structure, in
the proportion set by `covariate_cor` – which is the point of them. A
covariate unrelated to yield is not worth zoning on, and one identical
to it is an unrealistically easy test; the default puts elevation and
soil part-way, as they are.

Optional **pseudo-environments** cut the paddock into bands across the
strips. Each band gets its own treatment response through `zone_effect`,
so a zone-by-treatment interaction is genuinely there to be found – or,
with the default of one zone, genuinely is not, which is the more
important case to be able to test.

Given `treatments`, a trial is laid out with
[`make_trial_design()`](https://www.zcao.space/ofeIntegrateR/reference/make_trial_design.md)
and harvested: yield is potential, plus the treatment response for that
cell's zone, plus noise.

## Using it

The truth is attached as `attr(, "truth")` – the zone breaks, the
per-zone treatment effects, the variance of each component – so a check
can be made against what was simulated rather than against what looks
plausible. That is what separates this from a demonstration dataset.

## See also

[`simulate_ofe_trial()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_ofe_trial.md)
for a bare lattice and
[`simulate_yield_monitor()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_yield_monitor.md)
for the irregular shape a harvester records.

## Examples

``` r
p <- simulate_paddock(n_row = 30, n_col = 20, n_zones = 2,
                      treatments = c("N0", "N60", "N120"), seed = 1)
head(p)
#>   row col  x  y yield_potential elevation     clay zone treat rep plot    yield
#> 1   1   1 10 10        3.049526  100.9558 32.35223    1  N120   1    1 4.048992
#> 2   1   2 20 10        3.149135  101.5687 31.04612    1    N0   1    2 3.028595
#> 3   1   3 30 10        2.917562  102.7215 31.38791    1   N60   1    3 4.108944
#> 4   1   4 40 10        2.778075  103.6573 31.80192    1    N0   2    4 2.890744
#> 5   1   5 50 10        2.420393  101.5452 25.11243    1  N120   2    5 2.737379
#> 6   1   6 60 10        2.785883  102.2300 28.25401    1   N60   2    6 3.107974
attr(p, "truth")$zone_effect
#> [1] 1 1

# The zoning functions can now be checked against a known answer
z <- partition_paddock(p, covariates = c("elevation", "clay"))
table(z$zone, p$zone)
#>    
#>       1   2
#>   1 120 150
#>   2   0  75
#>   3  75  25
#>   4 105   0
#>   5   0  50
```
