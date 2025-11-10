# Internal function: Matrix of calculations for density gating of intensity values in `marker` for each unique subset of `subset_col`

Internal function for `getDensityGates` For each unique value in
`subset_col`, there is a matrix for storing calculations for density
gating contains: first to fourth derivatives of density, indicators for
local peaks, "real peaks", plateau_pre and cutoff

## Usage

``` r
getDensityMats(
  intens_dat,
  marker,
  subset_col,
  bin_n = 512,
  peak_detect_ratio = 10,
  pos_peak_threshold = 1800
)
```

## Arguments

- intens_dat:

  dataframe of pre-gated (compensated, biexp. transf, gated CD4/CD8)
  intensity values where cols = intensity value per marker, rows = each
  sample

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

  numeric for threshold to identify a positive peak ' default is 1800,
  which is on the biexponential scale

## Value

tibble of matrices for `marker` containing calculations for density
gating for each unique subset found in `subset_col`  
rows correspond to unique values in `subset_col`,  
cols correspond to the information for density gating
