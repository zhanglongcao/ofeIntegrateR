# ofeIntegrateR

<!-- badges: start -->
[![R-CMD-check](https://github.com/zhanglongcao/ofeIntegrateR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zhanglongcao/ofeIntegrateR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

A workflow for on-farm experimentation (OFE) strip trials, from planning a
trial to estimating its treatment effects — including integrating sparse
point-source measurements (soil cores, tissue samples, disease ratings) with
dense spatial layers (yield maps, EM38 surveys).

Developed for the AAGI-CU-RD-OFE GRDC project ("Development of processes to
integrate point-source data and high-resolution data").

**Every step has an open-source path.** ASReml-R is used automatically when
licensed, but nothing here requires a commercial licence.

## The pipeline

| Stage | Function |
|---|---|
| Lay out the trial: strip or stacked, randomised or systematic | `make_trial_design()` |
| Work out how many cores the target precision needs | `kriging_sample_interval()` |
| Decide where the cores go | `place_point_samples()` |
| Get irregular yield-monitor data onto an estimable lattice | `grid_dense_layer()` |
| Interpolate the sparse layer onto that lattice | `krige_point_samples()` |
| Check the point layer is dense enough to be worth using | `cv_krige_surface()` |
| Fit the baseline and the integrated model, and compare | `compare_integration()` |
| Fit either model on its own | `fit_integrated_kriged()` |
| Model both layers jointly instead | `fit_integrated_joint()` |
| Pull treatment contrasts out of any of them | `extract_fixed_effects()` |

Two simulators generate test data: `simulate_ofe_trial()` for a tidy lattice,
and `simulate_yield_monitor()` for the awkward shape a real harvester produces —
GPS-referenced points along passes, position error, a clipped paddock corner and
missing passes.

Two articles:

* `vignette("ofe-integration")` — an end-to-end walkthrough on simulated
  trials, where the truth is known.
* `vignette("real-data")` — the same workflow on
  `agridat::lasrosas.corn`, an on-farm nitrogen experiment from Argentina
  recorded by a yield monitor, where the data arrive in the state real data
  arrive in.

## Two integration strategies

Both have a fully open-source path as well as an
[asreml-R](https://vsni.co.uk/software/asreml-r) path for users with a licence:

- **Kriged-covariate** (`fit_integrated_kriged()`): krige the sparse
  point-source variable onto the trial grid with `krige_point_samples()`,
  then include the kriged surface as a fixed covariate in a spatial model
  of the dense response. Pass `covariate = NULL` for the dense-layer-only
  baseline, and several covariate names to integrate more than one point
  variable at once.
- **Joint bivariate model** (`fit_integrated_joint()`): fit the dense
  response and the point-source measurements *jointly* as two response
  layers/traits sharing a spatial random effect, avoiding a separate
  kriging step. `engine = "asreml"` fits this via asreml-R;
  `engine = "sommer"` (open source) fits an analogous multi-trait model
  via `sommer::mmer()`. In testing, the sommer fit reproduced the asreml
  treatment-contrast estimates almost exactly.

## Choosing an engine

Real strip trials are replicated, so the analysis usually needs a random
effect. Both `asreml` and `lme` take one, with the same spelling:

| `engine` | Fitted by | Random effects | Spatial residual | Licence |
|---|---|---|---|---|
| `"asreml"` | `asreml::asreml()` | yes | `ar1(row):ar1(col)` | commercial |
| `"lme"` | `nlme::lme()` | yes | exponential + nugget | open source |
| `"gls"` | `nlme::gls()` | no | exponential + nugget | open source |
| `"lm"` | `stats::lm()` | no | none | open source |

```r
# The same call on either engine
fit_integrated_kriged(kr, response = "yield", treat = "treat",
                      covariate = "soil_n_kriged",
                      random = ~ rep, engine = "lme")     # or "asreml"
```

`"lme"` is the open-source counterpart to the asreml fit: random effects *and*
a spatial correlation structure. `"gls"` has the spatial structure but no
random effects, and `"lm"` has neither — asking either of them for a random
effect is an error rather than something silently dropped.

For the joint bivariate model, `fit_integrated_joint()` uses
`asreml::asreml()` or the open-source `sommer::mmer()`.

`extract_fixed_effects()` returns the same tidy table for an `lm`, `gls`,
`lme`, `asreml` or `mmer` fit, so downstream code does not change when a
licence appears or disappears.

## Always report the baseline

The dense layer is already collected, so a spatial model of it costs nothing;
point sampling does. Reporting only the integrated fit hides how much of its
accuracy was available for free. `compare_integration()` fits both and reports
the treatment contrasts side by side.

A large, highly significant covariate coefficient is **not** evidence that
integration helped: with hundreds of grid cells the covariate clears
conventional significance almost automatically. Judge the method by how the
treatment contrasts and their standard errors move.

## Installation

```r
# install.packages("remotes")
remotes::install_github("zhanglongcao/ofeIntegrateR")
```

Documentation: <https://www.zcao.space/ofeIntegrateR/>

Everything works out of the box with open-source dependencies only
(`engine = "lm"` / `"gls"` for `fit_integrated_kriged()`, and
`engine = "sommer"` for `fit_integrated_joint()`). The `engine = "asreml"`
options additionally require [asreml-R](https://vsni.co.uk/software/asreml-r),
a commercial package from VSNi that is **not** installed automatically —
see VSNi for licensing.

## Example

```r
library(ofeIntegrateR)

# A trial as a yield monitor actually records it: points along passes,
# GPS jitter, a clipped corner, some passes missing
sim <- simulate_yield_monitor(n_point_samples = 30, seed = 11)

# 1. Snap the irregular cloud onto a complete lattice
g <- grid_dense_layer(sim$cloud, response = "yield", treat = "treat",
                      cell_size = 9)

# Drop cells straddling a treatment boundary, keeping the lattice complete
g$yield[!is.na(g$treat_purity) & g$treat_purity < 0.8] <- NA

# 2. Krige the sparse point layer onto the same grid
g$x <- g$x_centre; g$y <- g$y_centre
kr <- krige_point_samples(sim$point_samples, g, value = "point_obs",
                          coords = c("x", "y"))

# 3. Is the point layer dense enough to be worth using?
cv_krige_surface(sim$point_samples, value = "point_obs", coords = c("x", "y"))

# 4. Baseline vs integrated, side by side
compare_integration(kr, response = "yield", treat = "treat",
                    covariate = "point_obs_kriged", engine = "gls")
```

## Notes on asreml-R behaviour

A few asreml-R quirks this package works around, documented here so they
aren't rediscovered the hard way:

- `summary(fit)$coef.fixed` returns `NULL` in asreml-R >= 4.2; you need
  `summary(fit, coef = TRUE)$coef.fixed`. `extract_fixed_effects()`
  handles this for you.
- `ai.sing` is a *session-level* option (`asreml::asreml.options(ai.sing = TRUE)`),
  not an argument to `asreml()` — passing `ai.sing = TRUE` directly to
  `asreml()` is silently ignored. `fit_integrated_joint()` sets it
  temporarily (restored on exit) via the `ai_sing` argument (default
  `TRUE`), since the shared-layer random effect can trigger Average
  Information matrix singularities even with sparse, well-behaved data.
- An `ar1(row):ar1(col)` residual is generally unidentifiable for the
  point-source layer in the joint model, because a handful of sparse,
  irregularly placed observations cannot support a full spatial
  autocorrelation structure on their own. `fit_integrated_joint()` instead
  uses `dsum(~units | layer)` (separate i.i.d. residual variances per
  layer) and lets the shared random effect (`diag(layer)` or `us(layer)`)
  carry spatial information between layers.

## Notes on the sommer fallback

- `sommer::mmer()` defaults to `naMethodY = "exclude"`, which drops a row
  if *either* trait is missing — catastrophic here, since the
  point-source trait is `NA` almost everywhere. `fit_integrated_joint()`
  sets `naMethodY = "include"` so the dense response keeps using all its
  rows.
- sommer fits the treatment effect separately for both the dense and
  point-source traits; only the dense-trait estimate is meaningful and
  should be used (`extract_fixed_effects()` returns terms for both, named
  e.g. `dense_response:treatB` and `point_obs:treatB` — ignore the latter).
- The sommer multi-trait model needs enough point-source observations
  relative to the grid size to avoid a singular system (`vsr(unit, ...)`
  fits one random-effect level per grid cell). A sparse setup that works
  fine in asreml (e.g. 200 grid cells, 12 points) can fail outright in
  sommer; a denser one (e.g. 630 grid cells, 25 points) works reliably in
  testing.

## License

MIT © Zhanglong Cao
