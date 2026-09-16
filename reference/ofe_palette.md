# Colours used by the package's plots

Exposed so that a report can match its figures to the rest of its
styling, and so the choices can be inspected rather than taken on trust.

## Usage

``` r
ofe_palette(
  what = c("treatment", "zone", "surface", "residual", "chrome"),
  n = 5
)
```

## Arguments

- what:

  Which scale: `"treatment"`, `"zone"`, `"surface"`, `"residual"`, or
  `"chrome"` for the ink, grid and surface colours.

- n:

  Number of colours wanted. Ignored for `"chrome"`.

## Value

A character vector of hex colours; for `"chrome"`, a named vector.

## Details

Three jobs, three kinds of scale, because the job decides the scale:

- `"treatment"`:

  Categorical – treatments have identity, not magnitude. The order is
  fixed and nested, so adding a treatment never repaints the others.
  Validated pairwise for up to six levels: every pair that can appear on
  screen together is separable both to normal colour vision and under
  simulated deuteranopia, protanopia and tritanopia. Past six the
  guarantee lapses, and the plots fall back to labelling.

- `"zone"`:

  Ordinal – zones from
  [`partition_paddock()`](https://www.zcao.space/ofeIntegrateR/reference/partition_paddock.md)
  and
  [`partition_pseudo_env()`](https://www.zcao.space/ofeIntegrateR/reference/partition_pseudo_env.md)
  are ordered (by covariate mean, or along the trial), so a single hue
  running light to dark says so, where a categorical set would imply
  they are unrelated. Steps are spaced to stay distinguishable up to
  seven zones – past that the fill can no longer separate them and the
  zone label has to carry identity.

- `"surface"`:

  Sequential, for a continuous map – a kriged layer, a yield surface,
  elevation. The same hue as `"zone"` but running from a lighter start,
  because a continuous scale is read as a gradient rather than as a set
  of levels to tell apart, and the extra range buys resolution.

- `"residual"`:

  Diverging – a residual has a sign and a natural zero, so two hues meet
  at a neutral grey. Never a rainbow: a hue at the midpoint invents a
  category where the data has nothing.

## Examples

``` r
ofe_palette("treatment", 3)
#> [1] "#2a78d6" "#1baf7a" "#eda100"
ofe_palette("zone", 4)
#> [1] "#85B6EE" "#5F88C1" "#395D94" "#0C356B"
ofe_palette("surface", 6)
#> [1] "#CCE1FA" "#A7BCDC" "#8398BF" "#5F75A2" "#3B5486" "#0C356B"
ofe_palette("chrome")
#>   surface       ink      ink2     muted      grid      axis  midpoint 
#> "#fcfcfb" "#0b0b0b" "#52514e" "#898781" "#e1e0d9" "#c3c2b7" "#f0efec" 
```
