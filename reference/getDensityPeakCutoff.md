# Internal function: Determine the "real peaks" and cutoff based on the density estimation and its derivs

Internal function for `getDensityGates`

## Usage

``` r
getDensityPeakCutoff(
  dens_binned_dat,
  marker,
  subset_col,
  bin_n = 512,
  peak_detect_ratio = 10,
  pos_peak_threshold = 1800,
  dens_flip = FALSE
)
```

## Arguments

- dens_binned_dat:

  list of dataframe output from the `getDensityDerivs`

- marker:

  string for the marker to gate on the name needs to match exactly the
  column name in `dens_binned_dat`

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

- dens_flip:

  logical for whether the gating should be applied "backwards" where the
  peak is a positive peak and want to gate to the left of peak instead
  of right

## Value

list of dataframe `dens_binned_dat` with additional columns added for
peak(s) identified and the cutoff each element corresponds to each
unique value of `subset_col` for each dataframe: rows correspond to each
of the bins
