# Density gating of intensity values in `marker` for each unique subset of `subset_col`

For each unique value in `subset_col`, gate using density and estimated
derivatives to identify cutoff at shoulder (i.e., point of tapering off)
relative to the peak for `marker` (intensity values). The strategy of
cutting at the shoulder mimics the strategy to gate relative to a
unimodal background negative subpopulation, which is capable of
capturing dim subpopulations.

## Usage

``` r
getDensityGates(
  intens_dat,
  marker,
  subset_col,
  bin_n = 512,
  peak_detect_ratio = 10,
  pos_peak_threshold = 1800,
  neg_intensity_threshold = -1000
)
```

## Arguments

- intens_dat:

  dataframe of pre-gated (compensated, biexp. transf, openCyto steps)
  intensity values where cols=intensity value per marker, rows=each
  sample

- marker:

  string for the marker(s) to gate on the names need to match exactly
  the column name in `intens_dat`

- subset_col:

  string for the column name to indicate the subsets to apply density
  gating on will perform operation on subsets corresponding to each
  unique value in column

- bin_n:

  numeric to be passed to `n` parameter of `density(n=bin_n)` for number
  of equally spaced points at which the density is to be estimated  
  Default is 512, which is the default of `density(n=512)`

- peak_detect_ratio:

  numeric threshold for eliminating small peaks where a peak that is \<
  than the highest peak by `peak_detect_ratio` times will be ignored  
  Default=10

- pos_peak_threshold:

  either:

  - numeric for threshold to identify a positive peak for all or

  - a dataframe if supplying multiple `marker` to gate. The dataframe
    needs to be supplied with 2 columns named `marker` and
    `pos_peak_threshold` and rows for the `marker` to gate

  Default is 1800 (note this is on the biexponential scale) for all
  `marker`

- neg_intensity_threshold:

  numeric for threshold to filter out any "very negatively" expressed
  cells in the density estimation to avoid over-compression and
  difficulty in distinguishing peaks and the gates  
  This is only applied as a filter for the density estimation, the cells
  \< `neg_intensity_threshold` are retained in the intensity matrix for
  other steps  
  Expects the `neg_intensity_threshold` is on the same scale as the
  transformed data in `intens_dat`  
  Default is `NULL`: no filters applied and density estimation based on
  all cells in corresponding subsets.  
  Suggested for biexp. transformed data is -1000 which corresponds to
  ~-3300 on the original intensity scale)

## Value

tibble of gates/cutoffs for `marker` for each unique subset found in
`subset_col` where

- rows correspond to unique values in `subset_col`

- , columns correspond to`marker`

## Examples

``` r
# Create a fake dataset
set.seed(100)
intens_dat<-tibble::tibble(
               CD3_pos=rep(c(0, 1), each=50),
               CD4=rnorm(100, 100, 10),
               CD8=rnorm(100, 100, 10)
)

# Run density gating, leaving other params at suggested defaults
# number of bins suggested is 40 but default is at `bin_n=512`,
# which is the default for the R base density() function
getDensityGates(intens_dat, marker="CD4", subset_col="CD3_pos", bin_n=40)
#> # A tibble: 2 × 2
#>   CD3_pos   CD4
#>     <dbl> <dbl>
#> 1       0  88.9
#> 2       1 113. 
```
