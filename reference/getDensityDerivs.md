# Internal function: Estimate derivatives for density of `marker` for each unique subset of `subset_col`

Internal function for `get_density_gates` For each unique value in
`subset_col`, estimate the derivatives for `marker` (intensity values)

## Usage

``` r
getDensityDerivs(
  dens,
  marker,
  subset_col,
  bin_n = 512,
  peak_detect_ratio = 10,
  pos_peak_threshold = 1800
)
```

## Arguments

- dens:

  `density` object from the
  [density](https://rdrr.io/r/stats/density.html)

- marker:

  string for the marker to gate on the name needs to match exactly the
  column name in `intens_dat`

- subset_col:

  string for the column name to indicate the subsets to apply density
  gating on will perform operation on subsets corresponding to each
  unique value in column

- bin_n:

  numeric to be passed to `n` parameter of `density(n = bin_n)` for
  number of equally spaced points at which the density is to be
  estimated default is 512, which is the default of `density(n = 512)`

- peak_detect_ratio:

  numeric threshold for eliminating small peaks where a peak that is \<
  than the highest peak by `peak_detect_ratio` times will be ignored
  default = 10

- pos_peak_threshold:

  numeric for threshold to identify a positive peak default is 1800,
  which is on the biexponential scale

## Value

list of dataframe with density estimation and corresponding 1st-4th
derivatives, indicators of local peaks, plateau_pre  
each element corresponds to each unique value of `subset_col`  
for each dataframe: rows correspond to each of the bins
