# How far apart should the samples be, and how many are needed?

Works out the grid interval at which ordinary kriging attains a target
precision, given a variogram, and converts that interval into a sample
count for a trial of known area. This is the question that has to be
answered before the sampling budget is set, and it cannot be answered by
a rule of thumb: the same sample count is generous on a paddock that
varies over hundreds of metres and useless on one that turns over every
ten.

## Usage

``` r
kriging_sample_interval(
  nugget,
  psill,
  range,
  target_kse = 0.5,
  area_ha = NULL,
  model = "Exp",
  intervals = NULL
)
```

## Arguments

- nugget, psill, range:

  Variogram parameters: nugget \\c_0\\, partial sill \\c_1\\, and
  practical range, in the units of the trial's coordinates (normally
  metres).

- target_kse:

  Target kriging standard error as a fraction of the field standard
  deviation. Must exceed `kse_floor` to be attainable.

- area_ha:

  Optional trial area in hectares. When given, the sample count implied
  by the interval is returned as well.

- model:

  Variogram model passed to
  [`gstat::vgm()`](https://r-spatial.github.io/gstat/reference/vgm.html).

- intervals:

  Optional numeric vector of intervals to evaluate. Defaults to a
  sequence spanning a twentieth of the range to twice the range.

## Value

A list with the largest `interval` meeting the target, the `rel_kse`
achieved there, `n_samples` for `area_ha` (or `NA`), the `target_kse`,
the `kse_floor` set by the nugget, and a `curve` data frame of interval
against relative kriging standard error, for plotting the trade-off.
`interval` is `NA` when the target lies below the floor.

## Details

Precision is expressed as the **kriging standard error relative to the
field standard deviation**, \\\sigma_K / \sqrt{c_0 + c_1}\\. A relative
error of 0.5 means the interpolated surface is twice as precise as
guessing the field mean everywhere.

There is a floor. At an unsampled location the nugget component cannot
be filtered out however densely you sample, so no interval attains a
relative error below \\\sqrt{c_0 / (c_0 + c_1)}\\. A property with a
nugget-to-sill ratio of 0.3 cannot be mapped to better than 0.55 of the
field standard deviation at any cost. The returned `kse_floor` reports
this, and it is often the most useful number here: it says whether the
target is worth budgeting for at all.

The calculation places four samples at the corners of a \\\Delta \times
\Delta\\ cell and predicts its centre — the worst-supported point on a
regular grid — so the answer is conservative by construction.

Variogram parameters normally come from a reconnaissance survey or from
ancillary data such as EM38 or several seasons of yield maps. Where no
prior exists, run the calculation across a plausible range and plan for
the shortest range you might encounter.

## Examples

``` r
# A medium-range soil property over a 20 ha trial
s <- kriging_sample_interval(nugget = 0.2, psill = 0.8, range = 60,
                             target_kse = 0.6, area_ha = 20)
s$interval
#> [1] 12.91525
s$n_samples
#> [1] 1200

# Short-range variation needs far denser sampling for the same precision
kriging_sample_interval(nugget = 0.2, psill = 0.8, range = 15,
                        target_kse = 0.6, area_ha = 20)$n_samples
#> [1] 19185

# The nugget sets a floor no sampling effort can beat
kriging_sample_interval(nugget = 0.5, psill = 0.5, range = 60)$kse_floor
#> Warning: A relative kriging standard error of 0.5 is unattainable at any sampling density: the nugget sets a floor of 0.707. Either accept a larger target or reduce the nugget, e.g. by compositing cores.
#> [1] 0.7071068
```
