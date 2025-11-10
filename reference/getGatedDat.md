# Attach indicator columns to `intens_dat` based on gates provided in `cutoffs`

Adds an indicator column (0/1) to `intens_dat` for each marker in
`cutoffs` as indicated by the columns in `cutoffs`

## Usage

``` r
getGatedDat(intens_dat, cutoffs, subset_col)
```

## Arguments

- intens_dat:

  dataframe of pre-gated (compensated, biexp. transf, openCyto steps)
  intensity values where rows=each cell and cols are the intensity
  values for each marker

- cutoffs:

  tibble of gates/cutoffs for all markers to gate  
  Expects `cutoffs` to match format of output from
  [`getDensityGates()`](https://leejasme.github.io/staRgate/reference/getDensityGates.md)
  with column corresponding to a marker, and rows to the subsets defined
  in the `subset_col`

- subset_col:

  string for the column name to indicate the subsets to apply density
  gating on will perform operation on subsets corresponding to each
  unique value in column

## Value

`intens_dat` with additional columns attached for each marker in
`cutoffs`

## Details

The naming convention for the tagged on indicator columns will be
`tolower(<marker_name>_pos)` where 0 indicates negativity or intensity
\< gate provided 1 indicates positivity or intensity \> gate provided

## Examples

``` r
# Create a fake dataset
set.seed(100)
intens_dat <- tibble::tibble(
               CD3_pos=rep(c(0, 1), each=50),
               CD4=rnorm(100, 100, 10),
               CD8=rnorm(100, 100, 10)
)

# Run getDensityGates to obtain the gates
gates <- getDensityGates(intens_dat, marker="CD4", subset_col="CD3_pos", bin_n=40)

# Tag on the 0/1 on intens_dat
intens_dat_2 <- getGatedDat(intens_dat, cutoffs=gates, subset_col="CD3_pos")

# intens_dat_2 now has the cd4_pos tagged on
head(intens_dat_2)
#> # A tibble: 6 × 4
#>   CD3_pos   CD4   CD8 cd4_pos
#>     <dbl> <dbl> <dbl>   <dbl>
#> 1       0  95.0  96.7       1
#> 2       0 101.  114.        1
#> 3       0  99.2  95.3       1
#> 4       0 109.  108.        1
#> 5       0 101.   85.4       1
#> 6       0 103.   96.0       1
```
