# ofeIntegrateR

A workflow for on-farm experimentation (OFE) strip trials, from planning
a trial to estimating its treatment effects — including integrating
sparse point-source measurements (soil cores, tissue samples, disease
ratings) with dense spatial layers (yield maps, EM38 surveys).

Developed for the AAGI-CU-RD-OFE GRDC project (“Development of processes
to integrate point-source data and high-resolution data”).

**Every step has an open-source path**, including the analysis:
[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
is a REML engine with ASReml-style separable residual structures written
in base R, and it reproduces `asreml::asreml()` to several significant
figures. ASReml-R is used where licensed, but nothing here requires a
commercial licence.

## The pipeline

| Stage | Function |
|----|----|
| Lay out the trial: strip or stacked, randomised or systematic | [`make_trial_design()`](https://www.zcao.space/ofeIntegrateR/reference/make_trial_design.md) |
| Work out how many cores the target precision needs | [`kriging_sample_interval()`](https://www.zcao.space/ofeIntegrateR/reference/kriging_sample_interval.md) |
| Decide where the cores go | [`place_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/place_point_samples.md) |
| Get irregular yield-monitor data onto an estimable lattice | [`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md) |
| Interpolate the sparse layer onto that lattice | [`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md) |
| Fit and draw the variogram of the sampled layer | [`ofe_variogram()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_variogram.md) |
| Check the point layer is dense enough to be worth using | [`cv_krige_surface()`](https://www.zcao.space/ofeIntegrateR/reference/cv_krige_surface.md) |
| Find the pseudo-environments the spatial covariance supports | [`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md) |
| Zone a paddock into contiguous regions from elevation and soil | [`partition_paddock()`](https://www.zcao.space/ofeIntegrateR/reference/partition_paddock.md) |
| Write the residual structure those zones imply | [`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md) |
| Fit the spatial mixed model by REML, no licence needed | [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md) |
| Test the fixed terms, and get predicted treatment means | [`wald_tests()`](https://www.zcao.space/ofeIntegrateR/reference/wald_tests.md), [`ofe_means()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_means.md) |
| Report the means with an LSD and a/b/c letters | [`ofe_lsd()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_lsd.md) |
| Fit the baseline and the integrated model, and compare | [`compare_integration()`](https://www.zcao.space/ofeIntegrateR/reference/compare_integration.md) |
| Fit either model on its own | [`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md) |
| Model both layers jointly instead | [`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md) |
| Pull treatment contrasts out of any of them | [`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md) |
| Map anything over the trial, and diagnose the fit | [`ofe_map()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_map.md), [`plot()`](https://rdrr.io/r/graphics/plot.default.html) |

Three simulators generate test data.
[`simulate_paddock()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_paddock.md)
is the one to reach for: a correlated yield-potential surface, covariate
layers related to it by a correlation you choose, pseudo-environments
with their own treatment response, and a trial laid into it — with the
truth attached as `attr(, "truth")`, so the zoning and the analysis can
be *scored* rather than merely run. A covariate unrelated to yield is
not worth zoning on and one identical to it is an unrealistically easy
test, so the default puts elevation and soil part-way, as they are.
[`simulate_ofe_trial()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_ofe_trial.md)
gives a tidy lattice for quick tests, and
[`simulate_yield_monitor()`](https://www.zcao.space/ofeIntegrateR/reference/simulate_yield_monitor.md)
the awkward shape a real harvester produces — GPS-referenced points
along passes, position error, a clipped paddock corner and missing
passes.

Two articles:

- [`vignette("ofe-integration")`](https://www.zcao.space/ofeIntegrateR/articles/ofe-integration.md)
  — an end-to-end walkthrough on simulated trials, where the truth is
  known.
- [`vignette("real-data")`](https://www.zcao.space/ofeIntegrateR/articles/real-data.md)
  — the same workflow on
  [`agridat::lasrosas.corn`](https://kwstat.github.io/agridat/reference/lasrosas.corn.html),
  an on-farm nitrogen experiment from Argentina recorded by a yield
  monitor, where the data arrive in the state real data arrive in.

## Looking at it

A trial is a thing in a paddock, and most of what goes wrong in one is
obvious on a map and invisible in a table. Every stage has a
[`plot()`](https://rdrr.io/r/graphics/plot.default.html):

``` r

plot(make_trial_design(c("N0", "N60", "N120"), n_rep = 4))  # the plan
plot(ofe_variogram(cores, value = "soil_n"))                # the sampled layer
plot(partition_paddock(g, covariates = c("elevation", "clay")))   # the zones
plot(fit)                                    # residual map, Q-Q, residual variogram
ofe_map(g, "yield")                          # or map any column yourself
```

`plot(fit)` is the one to run before quoting a standard error. A spatial
model is fitted precisely to absorb the field’s pattern, so the question
is whether any is left: if the residual map still shows patches, or the
residual variogram is still climbing at short lags, the residual
structure has not done its job and the treatment standard errors are
optimistic.

The colours are not a matter of taste and were not picked by eye. Each
scale was run through a colour-vision validator against the surface an R
device actually draws on, and only sets clearing every gate are used —
the categorical order is pairwise-separable up to six treatments, the
zone ramp up to seven zones, and
[`ofe_map()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_map.md)
will not choose a diverging scale on its own, because a variable that
merely contains negative values is not thereby a signed one.
[`ofe_palette()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_palette.md)
exposes them.

## Analysis without asreml

[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
fits a linear mixed model by residual maximum likelihood with a
separable, ASReml-style residual structure — written in base R, so it
installs anywhere:

``` r

fit <- fit_ofe(yield ~ treat, random = ~ rep,
               residual = ~ ar1(row):ar1(col), data = g)

summary(fit)                      # variance parameters + fixed effects
wald_tests(fit)                   # like wald.asreml(): one test per term
ofe_means(fit, "treat")           # like predict.asreml(): treatment means
ofe_means(fit, "treat", pairwise = TRUE)   # contrasts with their SEDs
```

For the table that actually goes in the report — means, LSD, and the
letters next to them —
[`ofe_lsd()`](https://www.zcao.space/ofeIntegrateR/reference/ofe_lsd.md):

``` r

tab <- ofe_lsd(fit, "treat")            # means sorted best-first, with a/b/c
attr(tab, "lsd")                        # average SED, LSD, df, settings
attr(tab, "comparisons")                # every pairwise test behind the letters
```

      treat estimate std.error  lower  upper group
    1     F   14.820     0.375 14.063 15.577     a
    2     E   12.769     0.375 12.011 13.526     b
    3     D   11.854     0.375 11.096 12.611    bc
    4     C   11.511     0.375 10.754 12.268     c
    5     B   10.294     0.375  9.536 11.051     d
    6     A    9.674     0.375  8.917 10.432     d

On a balanced design this reproduces `agricolae::LSD.test()` exactly —
the same LSD value and the same letters — and `adjust = "tukey"`
reproduces `agricolae::HSD.test()`. Where it differs is on the data this
package is for: by default each pair is judged on **its own** standard
error of difference, because in a spatial model neighbouring strips
really are compared more precisely than distant ones. Pass `use = "lsd"`
for the single average LSD a published table usually means.

`by = "zone"` letters the means within each pseudo-environment
separately — comparisons never cross a zone, and each zone gets its own
LSD. And
[`compact_letters()`](https://www.zcao.space/ofeIntegrateR/reference/compact_letters.md)
will letter a pairwise table from anywhere else (`asreml::predict()`,
`emmeans`, a table typed by hand), so the same annotation can be put on
a fit this package did not produce.

The `residual` argument takes the asreml spellings, combined with `:`
for a separable structure:

| Spelling | Meaning |
|----|----|
| `id(f)`, or a bare `f` | independent |
| `ar1(f)` | first-order autoregressive along the ordered levels of `f` |
| `exp(x)` | exponential correlation in a numeric coordinate |
| `diag(f)` | a separate variance for each level of `f` |
| `dsum(~ struct \| s, levels = )` | independent sections, each with its own parameters |

On the structures both support,
[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
and `asreml::asreml()` agree to several significant figures — variance
components, fixed effects and their standard errors alike. Where a
licence is available, `engine = "asreml"` still fits the identical model
and is faster on large lattices; the engine here is dense, and
factorises an n×n matrix at every iteration, so it is meant for
OFE-sized problems (a few thousand lattice cells — aggregate with
[`grid_dense_layer()`](https://www.zcao.space/ofeIntegrateR/reference/grid_dense_layer.md)
if you have more). It does not implement `us()`/`fa()` structures or
multi-trait models.

One difference in its favour: cells with a missing response are simply
dropped. The correlation is evaluated from the row and column positions
of whatever observations remain, so an incomplete lattice needs no
padding — and a row absent from the data still counts as a lag rather
than being closed up.

## Pseudo-environments, derived rather than guessed

Strip trials are long, and a treatment effect at one end of the paddock
need not be the effect at the other.
[`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
cuts the trial into contiguous zones from the data: it removes the
treatment signal, collapses the residual field to a profile along the
trial, and segments that profile optimally by dynamic programming.

What stops it inventing zones is the spatial covariance. Any smooth
field looks like it has regions; the question is whether a region is
wider than the correlation range, because anything narrower is one
realisation of the same correlated surface rather than a distinct
environment. So the practical range of a fitted exponential variogram
becomes the minimum zone width, and caps the number of zones at
`floor(trial length / range)`. Within that cap, BIC chooses.

``` r

z <- partition_pseudo_env(g, response = "yield", along = "row", treat = "treat")
attr(z, "partition")$range      # what set the minimum zone width
attr(z, "partition")$zones      # where the cuts fell, and how big each zone is
attr(z, "partition")$bic        # the trade-off it chose from

r <- adaptive_residual(z)       # ~ dsum(~ ar1(row):ar1(col) | zone, levels = ...)
fit <- fit_ofe(yield ~ zone + zone:treat, data = z, residual = r)
```

[`adaptive_residual()`](https://www.zcao.space/ofeIntegrateR/reference/adaptive_residual.md)
exists because AR1 needs at least two levels in a dimension and real
zones do not always have them: a narrow zone may be one row deep, and
asking for `ar1()` there gives an unidentifiable parameter and a failed
fit. Each zone gets `ar1()` only where it has extent, and `id()`
elsewhere. The formula it writes is plain text that both
[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
and `asreml::asreml()` accept.

### Zoning from elevation and soil instead of yield

[`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
cuts across the trial, which is what a strip trial wants, and it works
from the response. When the zones should come from what is known
*before* harvest — elevation, an EM38 or gamma survey, a soil test grid
— and may be any shape, use
[`partition_paddock()`](https://www.zcao.space/ofeIntegrateR/reference/partition_paddock.md):

``` r

z <- partition_paddock(g, covariates = c("elevation", "ec_shallow", "clay"))
attr(z, "partition")$zones   # cells, area, patches, covariate means per zone
attr(z, "partition")$table   # ssd, r2, marginal gain and CH at each k
```

**The shape of a zone is a practical decision**, so there are three
methods and they differ only in that. Same ridge-and-sandy-corner
paddock, k = 4:

     rectangle (default)     skater                  kmeans
     222222222222111111      222222222222111111      111111111111222222
     222222222222111111      222222222222111111      111111111111222222
     222222222222111111      222222222233111111      131111313133222222
     444444444444444444      222222223333111111      333333333333222222
     444444444444444444      222222333333333333      333333333333333333
     444444444444444444      444422233333333333      313111111113131333
     333333333333333333      444422222333333333      111111111111111111
     333333333333333333      222222222222222222      111111111111111111

- **`"rectangle"` (default)** — the paddock is cut by lines running the
  full width or length of the region being split, each placed where it
  removes most within-zone variance. Cutting a rectangle across leaves
  two rectangles, so every zone is a rectangle however many cuts are
  made. This is what a variable-rate prescription, a sampling grid or a
  set of management blocks wants, and the zone summary gives each
  block’s corners (`x_min`…`y_max`) so it can be marked out. It is the
  two-dimensional version of what
  [`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
  does along one axis.
- **`"skater"`** — contiguous, but any shape: a minimum spanning tree
  over the neighbourhood graph, edges weighted by distance in
  standardised covariate space, pruned one edge at a time. Every zone is
  a subtree, so every zone is connected — but a zone may send a finger
  out between two others, as zones 2 and 3 do above. Right when the
  zones are there to describe the soil rather than to be worked, or when
  a boundary really runs at an angle. Note it does not automatically fit
  better for being freer: both methods are greedy, and on this paddock
  the rectangles explain more (r² 0.78 against 0.73). Compare the `r2`
  column rather than assuming.
- **`"kmeans"`** — not contiguous at all. Cells go to the zone their
  soil resembles wherever they sit, so zones arrive as confetti. Kept
  for comparison; `patches` in the zone summary counts the connected
  pieces of each zone, which is 1 throughout for the first two methods
  and more than 1 whenever k-means fragments.

What the choice costs depends on the shape of the underlying feature.
Variance explained at k = 2, on simulated paddocks with one clear
feature:

| feature                          | `"rectangle"` | `"skater"` |
|----------------------------------|---------------|------------|
| block boundaries on the axes     | 0.78          | 0.78       |
| a smooth gradient up the paddock | 0.69          | 0.67       |
| a boundary running diagonally    | 0.25          | **0.98**   |
| a round patch in the middle      | 0.13          | **0.97**   |

Rectangles cost nothing when the structure is blocky or a gradient, and
cost almost everything when the boundary runs at an angle or curves — a
creek line, a dune, a contour. With a feature like that, either use
`"skater"` and accept the shapes, or raise `k` so rectangles can
approximate the boundary in steps.

`covariates` is required and names columns of `data` — there is no
default, because which layers define a zone is agronomy rather than
something the function can guess:

``` r

partition_paddock(dat, covariates = c("elevation", "soil"))   # those two columns
```

Leave `k` unset and zones are added while each one earns its keep — the
default stops at the first `k` whose successor would explain less than
5% more of the covariate variance (`min_gain`). A paddock uniform in its
covariates comes back as one zone rather than being carved up anyway.
`criterion = "ch"` (Calinski–Harabasz) is available, with the caveat in
[`?partition_paddock`](https://www.zcao.space/ofeIntegrateR/reference/partition_paddock.md)
that on a smoothly varying covariate it often has no interior maximum
and just picks `k_max`.

Pass `treat` and each zone’s treatment coverage is reported, with a
warning if a zone is missing a level — a zone that does not contain
every treatment cannot support a `zone:treat` term, and that is better
found here than in a rank deficiency later.

Two cautions, both in
[`?partition_pseudo_env`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md).
A smooth field with no step at all will still be split about two-thirds
of the time under the defaults, so read a zone as “this part of the
paddock behaves differently”, not as evidence of a boundary. And zones
derived from the same yield data that then estimate zone-specific
treatment effects will overstate those differences, because the
boundaries were placed where the residuals already differed — use them
for the residual structure freely, and confirm a zone-by-treatment
interaction against an independent layer (elevation, EM38, a prior
season’s yield) by passing that layer as `response` instead.

## Two integration strategies

Both have a fully open-source path as well as an
[asreml-R](https://vsni.co.uk/software/asreml-r) path for users with a
licence:

- **Kriged-covariate**
  ([`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)):
  krige the sparse point-source variable onto the trial grid with
  [`krige_point_samples()`](https://www.zcao.space/ofeIntegrateR/reference/krige_point_samples.md),
  then include the kriged surface as a fixed covariate in a spatial
  model of the dense response. Pass `covariate = NULL` for the
  dense-layer-only baseline, and several covariate names to integrate
  more than one point variable at once.
- **Joint bivariate model**
  ([`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md)):
  fit the dense response and the point-source measurements *jointly* as
  two response layers/traits sharing a spatial random effect, avoiding a
  separate kriging step. `engine = "asreml"` fits this via asreml-R;
  `engine = "sommer"` (open source) fits an analogous multi-trait model
  via [`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html). In
  testing, the sommer fit reproduced the asreml treatment-contrast
  estimates almost exactly.

## Choosing an engine

Real strip trials are replicated, so the analysis usually needs a random
effect. Both `asreml` and `lme` take one, with the same spelling:

| `engine` | Fitted by | Random effects | Spatial residual | Licence |
|----|----|----|----|----|
| `"ofe"` (default) | [`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md) | yes | `ar1(row):ar1(col)` | open source |
| `"asreml"` | `asreml::asreml()` | yes | `ar1(row):ar1(col)` | commercial |
| `"lme"` | [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) | yes | exponential + nugget | open source |
| `"gls"` | [`nlme::gls()`](https://rdrr.io/pkg/nlme/man/gls.html) | no | exponential + nugget | open source |
| `"lm"` | [`stats::lm()`](https://rdrr.io/r/stats/lm.html) | no | none | open source |

``` r

# The same call on any engine
fit_integrated_kriged(kr, response = "yield", treat = "treat",
                      covariate = "soil_n_kriged", random = ~ rep)
```

`"ofe"` is the default because it fits the same `ar1(row):ar1(col)`
model asreml does and needs no licence. `"asreml"` fits it faster on a
large lattice. `"lme"` substitutes an exponential correlation for the
separable AR1, `"gls"` has the correlation but no random effects, and
`"lm"` has neither — asking either of the last two for a random effect
is an error rather than something silently dropped.

For a residual structure other than `ar1(row):ar1(col)` —
pseudo-environment sections, a heterogeneous variance — call
[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md)
directly with your own `residual` formula.

For the joint bivariate model,
[`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md)
uses `asreml::asreml()` or the open-source
[`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html).

[`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
returns the same tidy table for an `ofe_fit`, `lm`, `gls`, `lme`,
`asreml` or `mmer` fit, so downstream code does not change when a
licence appears or disappears.

## Always report the baseline

The dense layer is already collected, so a spatial model of it costs
nothing; point sampling does. Reporting only the integrated fit hides
how much of its accuracy was available for free.
[`compare_integration()`](https://www.zcao.space/ofeIntegrateR/reference/compare_integration.md)
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

Documentation: <https://www.zcao.space/ofeIntegrateR/>

Everything works out of the box with open-source dependencies only: the
default `engine = "ofe"` for
[`fit_integrated_kriged()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_kriged.md)
and
[`fit_ofe()`](https://www.zcao.space/ofeIntegrateR/reference/fit_ofe.md),
and `engine = "sommer"` for
[`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md).
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
                    covariate = "point_obs_kriged")

# 5. Or fit the spatial model directly, with the structure you want
fit <- fit_ofe(yield ~ treat, data = kr, residual = ~ ar1(row):ar1(col))
wald_tests(fit)
ofe_means(fit, "treat", pairwise = TRUE)
```

## Notes on asreml-R behaviour

A few asreml-R quirks this package works around, documented here so they
aren’t rediscovered the hard way:

- `summary(fit)$coef.fixed` returns `NULL` in asreml-R \>= 4.2; you need
  `summary(fit, coef = TRUE)$coef.fixed`.
  [`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
  handles this for you.
- `ai.sing` is a *session-level* option
  (`asreml::asreml.options(ai.sing = TRUE)`), not an argument to
  `asreml()` — passing `ai.sing = TRUE` directly to `asreml()` is
  silently ignored.
  [`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md)
  sets it temporarily (restored on exit) via the `ai_sing` argument
  (default `TRUE`), since the shared-layer random effect can trigger
  Average Information matrix singularities even with sparse,
  well-behaved data.
- An `ar1(row):ar1(col)` residual is generally unidentifiable for the
  point-source layer in the joint model, because a handful of sparse,
  irregularly placed observations cannot support a full spatial
  autocorrelation structure on their own.
  [`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md)
  instead uses `dsum(~units | layer)` (separate i.i.d. residual
  variances per layer) and lets the shared random effect (`diag(layer)`
  or `us(layer)`) carry spatial information between layers.

## Notes on the sommer fallback

- [`sommer::mmer()`](https://rdrr.io/pkg/sommer/man/mmer.html) defaults
  to `naMethodY = "exclude"`, which drops a row if *either* trait is
  missing — catastrophic here, since the point-source trait is `NA`
  almost everywhere.
  [`fit_integrated_joint()`](https://www.zcao.space/ofeIntegrateR/reference/fit_integrated_joint.md)
  sets `naMethodY = "include"` so the dense response keeps using all its
  rows.
- sommer fits the treatment effect separately for both the dense and
  point-source traits; only the dense-trait estimate is meaningful and
  should be used
  ([`extract_fixed_effects()`](https://www.zcao.space/ofeIntegrateR/reference/extract_fixed_effects.md)
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
