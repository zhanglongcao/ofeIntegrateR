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

**Every step has an open-source path**, including the analysis: `fit_ofe()` is
a REML engine with ASReml-style separable residual structures written in base
R, and it reproduces `asreml::asreml()` to several significant figures. ASReml-R
is used where licensed, but nothing here requires a commercial licence.

## The pipeline

| Stage | Function |
|---|---|
| Lay out the trial: strip or stacked, randomised or systematic | `make_trial_design()` |
| Work out how many cores the target precision needs | `kriging_sample_interval()` |
| Decide where the cores go | `place_point_samples()` |
| Get irregular yield-monitor data onto an estimable lattice | `grid_dense_layer()` |
| Interpolate the sparse layer onto that lattice | `krige_point_samples()` |
| Check the point layer is dense enough to be worth using | `cv_krige_surface()` |
| Find the pseudo-environments the spatial covariance supports | `partition_pseudo_env()` |
| Write the residual structure those zones imply | `adaptive_residual()` |
| Fit the spatial mixed model by REML, no licence needed | `fit_ofe()` |
| Test the fixed terms, and get predicted treatment means | `wald_tests()`, `ofe_means()` |
| Report the means with an LSD and a/b/c letters | `ofe_lsd()` |
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

## Analysis without asreml

`fit_ofe()` fits a linear mixed model by residual maximum likelihood with a
separable, ASReml-style residual structure — written in base R, so it installs
anywhere:

```r
fit <- fit_ofe(yield ~ treat, random = ~ rep,
               residual = ~ ar1(row):ar1(col), data = g)

summary(fit)                      # variance parameters + fixed effects
wald_tests(fit)                   # like wald.asreml(): one test per term
ofe_means(fit, "treat")           # like predict.asreml(): treatment means
ofe_means(fit, "treat", pairwise = TRUE)   # contrasts with their SEDs
```

For the table that actually goes in the report — means, LSD, and the letters
next to them — `ofe_lsd()`:

```r
tab <- ofe_lsd(fit, "treat")            # means sorted best-first, with a/b/c
attr(tab, "lsd")                        # average SED, LSD, df, settings
attr(tab, "comparisons")                # every pairwise test behind the letters
```

```
  treat estimate std.error  lower  upper group
1     F   14.820     0.375 14.063 15.577     a
2     E   12.769     0.375 12.011 13.526     b
3     D   11.854     0.375 11.096 12.611    bc
4     C   11.511     0.375 10.754 12.268     c
5     B   10.294     0.375  9.536 11.051     d
6     A    9.674     0.375  8.917 10.432     d
```

On a balanced design this reproduces `agricolae::LSD.test()` exactly — the same
LSD value and the same letters — and `adjust = "tukey"` reproduces
`agricolae::HSD.test()`. Where it differs is on the data this package is for:
by default each pair is judged on **its own** standard error of difference,
because in a spatial model neighbouring strips really are compared more
precisely than distant ones. Pass `use = "lsd"` for the single average LSD a
published table usually means.

`by = "zone"` letters the means within each pseudo-environment separately —
comparisons never cross a zone, and each zone gets its own LSD. And
`compact_letters()` will letter a pairwise table from anywhere else
(`asreml::predict()`, `emmeans`, a table typed by hand), so the same annotation
can be put on a fit this package did not produce.

The `residual` argument takes the asreml spellings, combined with `:` for a
separable structure:

| Spelling | Meaning |
|---|---|
| `id(f)`, or a bare `f` | independent |
| `ar1(f)` | first-order autoregressive along the ordered levels of `f` |
| `exp(x)` | exponential correlation in a numeric coordinate |
| `diag(f)` | a separate variance for each level of `f` |
| `dsum(~ struct \| s, levels = )` | independent sections, each with its own parameters |

On the structures both support, `fit_ofe()` and `asreml::asreml()` agree to
several significant figures — variance components, fixed effects and their
standard errors alike. Where a licence is available, `engine = "asreml"` still
fits the identical model and is faster on large lattices; the engine here is
dense, and factorises an n×n matrix at every iteration, so it is meant for
OFE-sized problems (a few thousand lattice cells — aggregate with
`grid_dense_layer()` if you have more). It does not implement `us()`/`fa()`
structures or multi-trait models.

One difference in its favour: cells with a missing response are simply dropped.
The correlation is evaluated from the row and column positions of whatever
observations remain, so an incomplete lattice needs no padding — and a row
absent from the data still counts as a lag rather than being closed up.

## Pseudo-environments, derived rather than guessed

Strip trials are long, and a treatment effect at one end of the paddock need
not be the effect at the other. `partition_pseudo_env()` cuts the trial into
contiguous zones from the data: it removes the treatment signal, collapses the
residual field to a profile along the trial, and segments that profile
optimally by dynamic programming.

What stops it inventing zones is the spatial covariance. Any smooth field looks
like it has regions; the question is whether a region is wider than the
correlation range, because anything narrower is one realisation of the same
correlated surface rather than a distinct environment. So the practical range
of a fitted exponential variogram becomes the minimum zone width, and caps the
number of zones at `floor(trial length / range)`. Within that cap, BIC chooses.

```r
z <- partition_pseudo_env(g, response = "yield", along = "row", treat = "treat")
attr(z, "partition")$range      # what set the minimum zone width
attr(z, "partition")$zones      # where the cuts fell, and how big each zone is
attr(z, "partition")$bic        # the trade-off it chose from

r <- adaptive_residual(z)       # ~ dsum(~ ar1(row):ar1(col) | zone, levels = ...)
fit <- fit_ofe(yield ~ zone + zone:treat, data = z, residual = r)
```

`adaptive_residual()` exists because AR1 needs at least two levels in a
dimension and real zones do not always have them: a narrow zone may be one row
deep, and asking for `ar1()` there gives an unidentifiable parameter and a
failed fit. Each zone gets `ar1()` only where it has extent, and `id()`
elsewhere. The formula it writes is plain text that both `fit_ofe()` and
`asreml::asreml()` accept.

Two cautions, both in `?partition_pseudo_env`. A smooth field with no step at
all will still be split about two-thirds of the time under the defaults, so
read a zone as "this part of the paddock behaves differently", not as evidence
of a boundary. And zones derived from the same yield data that then estimate
zone-specific treatment effects will overstate those differences, because the
boundaries were placed where the residuals already differed — use them for the
residual structure freely, and confirm a zone-by-treatment interaction against
an independent layer (elevation, EM38, a prior season's yield) by passing that
layer as `response` instead.

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
| `"ofe"` (default) | `fit_ofe()` | yes | `ar1(row):ar1(col)` | open source |
| `"asreml"` | `asreml::asreml()` | yes | `ar1(row):ar1(col)` | commercial |
| `"lme"` | `nlme::lme()` | yes | exponential + nugget | open source |
| `"gls"` | `nlme::gls()` | no | exponential + nugget | open source |
| `"lm"` | `stats::lm()` | no | none | open source |

```r
# The same call on any engine
fit_integrated_kriged(kr, response = "yield", treat = "treat",
                      covariate = "soil_n_kriged", random = ~ rep)
```

`"ofe"` is the default because it fits the same `ar1(row):ar1(col)` model
asreml does and needs no licence. `"asreml"` fits it faster on a large lattice.
`"lme"` substitutes an exponential correlation for the separable AR1, `"gls"`
has the correlation but no random effects, and `"lm"` has neither — asking
either of the last two for a random effect is an error rather than something
silently dropped.

For a residual structure other than `ar1(row):ar1(col)` — pseudo-environment
sections, a heterogeneous variance — call `fit_ofe()` directly with your own
`residual` formula.

For the joint bivariate model, `fit_integrated_joint()` uses
`asreml::asreml()` or the open-source `sommer::mmer()`.

`extract_fixed_effects()` returns the same tidy table for an `ofe_fit`, `lm`,
`gls`, `lme`, `asreml` or `mmer` fit, so downstream code does not change when a
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

Everything works out of the box with open-source dependencies only: the
default `engine = "ofe"` for `fit_integrated_kriged()` and `fit_ofe()`, and
`engine = "sommer"` for `fit_integrated_joint()`. The `engine = "asreml"`
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
                    covariate = "point_obs_kriged")

# 5. Or fit the spatial model directly, with the structure you want
fit <- fit_ofe(yield ~ treat, data = kr, residual = ~ ar1(row):ar1(col))
wald_tests(fit)
ofe_means(fit, "treat", pairwise = TRUE)
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
