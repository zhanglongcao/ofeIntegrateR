# ofeIntegrateR

<!-- badges: start -->
<!-- R-CMD-check is paused while the repository is private; see the note in
     .github/workflows/R-CMD-check.yaml. Restore this badge when it is
     re-enabled -- a badge showing the status of an old run would claim the
     current code had been checked when it had not.
[![R-CMD-check](https://github.com/zhanglongcao/ofeIntegrateR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zhanglongcao/ofeIntegrateR/actions/workflows/R-CMD-check.yaml)
-->
<!-- badges: end -->

A workflow for on-farm experimentation (OFE) strip trials: lay out the design,
plan the sampling, get the yield monitor's point cloud onto an estimable
lattice, find the zones, fit the spatial model, and report the means. Every step
runs on open-source dependencies alone.

Developed for the AAGI-CU-RD-OFE GRDC project.

## Install

```r
# install.packages("remotes")
remotes::install_github("zhanglongcao/ofeIntegrateR")
```

Documentation: <https://www.zcao.space/ofeIntegrateR/>

## The workflow

| Stage | Function |
|---|---|
| Lay out the trial: strip or stacked, randomised or systematic | `make_trial_design()` |
| Work out how many cores a target precision needs | `kriging_sample_interval()` |
| Decide where the cores go | `place_point_samples()` |
| Get irregular yield-monitor data onto an estimable lattice | `grid_dense_layer()` |
| Fit and draw the variogram of a sampled layer | `ofe_variogram()` |
| Interpolate that layer onto the lattice | `krige_point_samples()` |
| Check the layer is dense enough to be worth using | `cv_krige_surface()` |
| Find the zones | `partition_pseudo_env()`, `partition_paddock()` |
| Write the residual structure those zones imply | `adaptive_residual()` |
| Fit the spatial mixed model | `fit_ofe()` |
| Test the terms; get means, contrasts and LSD letters | `wald_tests()`, `ofe_means()`, `ofe_lsd()` |
| Compare the integrated fit against the free baseline | `compare_integration()` |
| Map anything; diagnose the fit | `ofe_map()`, `plot()` |

Three simulators make test data: `simulate_paddock()` for a paddock whose truth
is known (yield surface, covariate layers, zones, trial, all recorded),
`simulate_ofe_trial()` for a tidy lattice, and `simulate_yield_monitor()` for
the shape a real harvester records.

## Example

```r
library(ofeIntegrateR)

sim <- simulate_yield_monitor(n_point_samples = 30, seed = 11)

# Snap the irregular cloud onto a complete lattice
g <- grid_dense_layer(sim$cloud, response = "yield", treat = "treat",
                      cell_size = 9)
g$yield[!is.na(g$treat_purity) & g$treat_purity < 0.8] <- NA
g$x <- g$x_centre; g$y <- g$y_centre

# What scale does the sampled layer vary at?
v <- ofe_variogram(sim$point_samples, value = "point_obs")
plot(v)
kriging_sample_interval(v, target_kse = 0.8, area_ha = 2.6)

# Krige it onto the same grid, then fit
kr <- krige_point_samples(sim$point_samples, g, value = "point_obs",
                          coords = c("x", "y"))
fit <- fit_ofe(yield ~ treat + point_obs_kriged, data = kr,
               residual = ~ ar1(row):ar1(col))

wald_tests(fit)
ofe_lsd(fit, "treat")
plot(fit)
```

## Two kinds of zone, for two different questions

The distinction matters and is easy to get wrong.

**`partition_pseudo_env()`** slices the trial **across its length**, from the
pattern in the yield layer. Strips run the length of the paddock, so cutting
across them keeps **every treatment in every zone** — which is what makes
`zone:treat` estimable. This is the one for a pseudo-environment analysis.

**`partition_paddock()`** clusters **environmental covariates** — elevation,
EM38, soil tests — into rectangular blocks (or contiguous free-form zones). It
never looks at yield, and its zones are whatever shape the soil is, so a zone
may not contain every treatment; pass `treat` and it reports the coverage and
warns you. This is the one for management blocks, variable-rate prescriptions,
and stratifying a sampling plan.

Either can feed `adaptive_residual()`, which writes the matching `dsum()`
residual structure.

## Always report the baseline

The dense layer is already collected, so a spatial model of it costs nothing;
point sampling does. Reporting only the integrated fit hides how much of its
accuracy was available for free. `compare_integration()` fits both and lines up
the treatment contrasts.

A large, highly significant covariate coefficient is **not** evidence that
integration helped: with hundreds of grid cells the covariate clears
conventional significance almost automatically. Judge the method by how the
treatment contrasts and their standard errors move.

## Fitting engines

`fit_ofe()` is a REML engine written in base R. It takes separable residual
structures — `id()`, `ar1()`, `exp()`, `diag()`, and `dsum(~ struct | s)`
sections — combined with `:`, and reproduces ASReml-R's estimates to several
significant figures on the structures both support.

`fit_integrated_kriged()` and `compare_integration()` take an `engine`:

| `engine` | Fitted by | Random effects | Spatial residual |
|---|---|---|---|
| `"ofe"` (default) | `fit_ofe()` | yes | `ar1(row):ar1(col)` |
| `"asreml"` | `asreml::asreml()` | yes | `ar1(row):ar1(col)` |
| `"lme"` | `nlme::lme()` | yes | exponential + nugget |
| `"gls"` | `nlme::gls()` | no | exponential + nugget |
| `"lm"` | `stats::lm()` | no | none |

Only `"asreml"` needs a licence; it is faster on a large lattice. For any
residual structure other than `ar1(row):ar1(col)`, call `fit_ofe()` directly.
`extract_fixed_effects()` returns the same tidy table whichever engine fitted
the model.

For the joint bivariate model, `fit_integrated_joint()` uses `asreml::asreml()`
or the open-source `sommer::mmer()`.

## Articles

* `vignette("ofe-integration")` — end to end on simulated trials, where the
  truth is known.
* `vignette("real-data")` — the same workflow on `agridat::lasrosas.corn`, an
  on-farm nitrogen experiment recorded by a yield monitor, where the data
  arrive in the state real data arrive in.

## License

MIT © Zhanglong Cao
