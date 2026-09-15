# ofeIntegrateR

Tools for integrating sparse point-source measurements (e.g. soil cores,
tissue samples, disease ratings) with dense spatial covariates
(e.g. EM38 surveys, yield maps) in on-farm experimentation (OFE) strip
trials.

Developed for the AAGI-CU-RD-OFE GRDC project (“Development of processes
to integrate point-source data and high-resolution data”).

**Every step has an open-source path.** ASReml-R is used automatically
when licensed, but nothing here requires a commercial licence.

## The pipeline

| Stage | Function |
|----|----|
| Decide where the cores go | [`place_point_samples()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/place_point_samples.md) |
| Get irregular yield-monitor data onto an estimable lattice | [`grid_dense_layer()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/grid_dense_layer.md) |
| Interpolate the sparse layer onto that lattice | [`krige_point_samples()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/krige_point_samples.md) |
| Check the point layer is dense enough to be worth using | [`cv_krige_surface()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/cv_krige_surface.md) |
| Fit the baseline and the integrated model, and compare | [`compare_integration()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/compare_integration.md) |
| Fit either model on its own | [`fit_integrated_kriged()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/fit_integrated_kriged.md) |
| Model both layers jointly instead | [`fit_integrated_joint()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/fit_integrated_joint.md) |
| Pull treatment contrasts out of any of them | [`extract_fixed_effects()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/extract_fixed_effects.md) |

Two simulators generate test data:
[`simulate_ofe_trial()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/simulate_ofe_trial.md)
for a tidy lattice, and
[`simulate_yield_monitor()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/simulate_yield_monitor.md)
for the awkward shape a real harvester produces — GPS-referenced points
along passes, position error, a clipped paddock corner and missing
passes.

Two articles:

- [`vignette("ofe-integration")`](https://zhanglongcao.github.io/ofeIntegrateR/articles/ofe-integration.md)
  — an end-to-end walkthrough on simulated trials, where the truth is
  known.
- [`vignette("real-data")`](https://zhanglongcao.github.io/ofeIntegrateR/articles/real-data.md)
  — the same workflow on
  [`agridat::lasrosas.corn`](https://kwstat.github.io/agridat/reference/lasrosas.corn.html),
  an on-farm nitrogen experiment from Argentina recorded by a yield
  monitor, where the data arrive in the state real data arrive in.

## Two integration strategies

Both have a fully open-source path as well as an
[asreml-R](https://vsni.co.uk/software/asreml-r) path for users with a
licence:

- **Kriged-covariate**
  ([`fit_integrated_kriged()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/fit_integrated_kriged.md)):
  krige the sparse point-source variable onto the trial grid with
  [`krige_point_samples()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/krige_point_samples.md),
  then include the kriged surface as a fixed covariate in a spatial
  model of the dense response. `engine = "asreml"` fits an
  `ar1(row):ar1(col)` residual via asreml-R; `engine = "gls"` (open
  source) fits an exponential spatial correlation structure via
  [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html);
  `engine = "lm"` (open source) ignores spatial autocorrelation entirely
  — the simplest but weakest fallback. Pass `covariate = NULL` for the
  dense-layer-only baseline, and several covariate names to integrate
  more than one point variable at once.
- **Joint bivariate model**
  ([`fit_integrated_joint()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/fit_integrated_joint.md)):
  fit the dense response and the point-source measurements *jointly* as
  two response layers/traits sharing a spatial random effect, avoiding a
  separate kriging step. `engine = "asreml"` fits this via asreml-R;
  `engine = "sommer"` (open source) fits an analogous multi-trait model
  via [`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html). In
  testing, the sommer fit reproduced the asreml treatment-contrast
  estimates almost exactly.

## Always report the baseline

The dense layer is already collected, so a spatial model of it costs
nothing; point sampling does. Reporting only the integrated fit hides
how much of its accuracy was available for free.
[`compare_integration()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/compare_integration.md)
fits both and reports the treatment contrasts side by side.

A large, highly significant covariate coefficient is **not** evidence
that integration helped: with hundreds of grid cells the covariate
clears conventional significance almost automatically. Judge the method
by how the treatment contrasts and their standard errors move.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("zhanglongcao/ofeIntegrateR")
```

Everything works out of the box with open-source dependencies only
(`engine = "lm"` / `"gls"` for
[`fit_integrated_kriged()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/fit_integrated_kriged.md),
and `engine = "sommer"` for
[`fit_integrated_joint()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/fit_integrated_joint.md)).
The `engine = "asreml"` options additionally require
[asreml-R](https://vsni.co.uk/software/asreml-r), a commercial package
from VSNi that is **not** installed automatically — see VSNi for
licensing.

## Example

``` r

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
aren’t rediscovered the hard way:

- `summary(fit)$coef.fixed` returns `NULL` in asreml-R \>= 4.2; you need
  `summary(fit, coef = TRUE)$coef.fixed`.
  [`extract_fixed_effects()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/extract_fixed_effects.md)
  handles this for you.
- `ai.sing` is a *session-level* option
  (`asreml::asreml.options(ai.sing = TRUE)`), not an argument to
  `asreml()` — passing `ai.sing = TRUE` directly to `asreml()` is
  silently ignored.
  [`fit_integrated_joint()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/fit_integrated_joint.md)
  sets it temporarily (restored on exit) via the `ai_sing` argument
  (default `TRUE`), since the shared-layer random effect can trigger
  Average Information matrix singularities even with sparse,
  well-behaved data.
- An `ar1(row):ar1(col)` residual is generally unidentifiable for the
  point-source layer in the joint model, because a handful of sparse,
  irregularly placed observations cannot support a full spatial
  autocorrelation structure on their own.
  [`fit_integrated_joint()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/fit_integrated_joint.md)
  instead uses `dsum(~units | layer)` (separate i.i.d. residual
  variances per layer) and lets the shared random effect (`diag(layer)`
  or `us(layer)`) carry spatial information between layers.

## Notes on the sommer fallback

- [`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html) defaults
  to `naMethodY = "exclude"`, which drops a row if *either* trait is
  missing — catastrophic here, since the point-source trait is `NA`
  almost everywhere.
  [`fit_integrated_joint()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/fit_integrated_joint.md)
  sets `naMethodY = "include"` so the dense response keeps using all its
  rows.
- sommer fits the treatment effect separately for both the dense and
  point-source traits; only the dense-trait estimate is meaningful and
  should be used
  ([`extract_fixed_effects()`](https://zhanglongcao.github.io/ofeIntegrateR/reference/extract_fixed_effects.md)
  returns terms for both, named e.g. `dense_response:treatB` and
  `point_obs:treatB` — ignore the latter).
- The sommer multi-trait model needs enough point-source observations
  relative to the grid size to avoid a singular system (`vsr(unit, ...)`
  fits one random-effect level per grid cell). A sparse setup that works
  fine in asreml (e.g. 200 grid cells, 12 points) can fail outright in
  sommer; a denser one (e.g. 630 grid cells, 25 points) works reliably
  in testing.

## License

MIT © Zhanglong Cao
