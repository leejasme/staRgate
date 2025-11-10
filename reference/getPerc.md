# Calculate the percentage of positive cells for specific subpopulations

Expects data input same as the output from `get_gated_dat` with
indicator columns of specific naming convention (see below).

## Usage

``` r
getPerc(
  intens_dat,
  num_marker,
  denom_marker,
  expand_num = FALSE,
  expand_denom = FALSE,
  keep_indicators = TRUE
)
```

## Arguments

- intens_dat:

  dataframe of gated data with indicator columns per marker of interest
  (specify in `num_marker` and `denom_marker`) with naming convention
  `marker_pos` per marker with values of 0 to indicate negative-, 1 to
  indicate positive-expressing

- num_marker:

  string for the marker(s) to specify the numerator for subpopulations
  of interest  
  See `expand_num` argument and examples for how to specify

- denom_marker:

  string for the marker(s) to specify the denominator for subpopulations
  of interest  
  See `expand_denom` argument and examples for how to specify.

- expand_num:

  logical, only accepts `TRUE` or `FALSE` with default of `FALSE`  
  if `expand_num=TRUE`, currently only considers up to pairs of markers
  specified in `num_marker` in the numerator of subpopulation
  calculations (e.g., CD4+ & CD8- of CD3+)  
  if `expand_num=FALSE`, only considers each marker specified in
  `num_marker` individually in the numerator of subpopulation
  calculations (e.g., CD4+ of CD3+)

- expand_denom:

  logical, only accepts `TRUE` or `FALSE` with default of `FALSE`  
  if `expand_denom=TRUE`, currently considers up to 1 marker from the
  `num_marker` and the unique combinations of `denom_marker` to generate
  list of subpopulations  
  e.g., if `denom_marker=c("CD8")`, `num_marker=c("LAG3", "KI67")`, and
  `expand_denom=TRUE`, the subpopulations will include:  
  1. LAG3+ of CD8+, LAG3- of CD8+, LAG3+ of CD8-, LAG3- of CD8-,  
  2. KI67+ of CD8+, KI67- of CD8+, KI67+ of CD8-, KI67- of CD8-,  
  3. KI67+ of (LAG3+ & CD8+), KI67- of (LAG3+ & CD8+), KI67+ of (LAG3+ &
  CD8-), KI67- of (LAG3+ & CD8-)...etc.,  
  4. LAG3+ of (KI67+ & CD8+), LAG3- of (KI67+ & CD8+), LAG3+ of (KI67+ &
  CD8-), LAG3- of (KI67+ & CD8-)...etc.,  
  if `expand_denom=FALSE`, only generates the list of subpopulations
  based on unique combinations of the `denom_marker` (e.g.,
  `denom_marker=c("CD4")` and `expand_denom=FALSE` only considers
  subpopulations with denominator CD4+ and CD4- whereas
  `denom_marker=c("CD4", "CD8"` and `expand_denom=FALSE` will consider
  subpopulations with denominators (CD4- & CD8-), (CD4+ & CD8-), (CD4- &
  CD8+) and (CD4+ & CD8+))

- keep_indicators:

  logical, only accepts `TRUE` or `FALSE` with default of `TRUE`  
  if `keep_indicators=TRUE`, will return indicator columns of 0/1 to
  specify which markers are considered in the numerator and denominators
  of the subpopulations.  
  Naming convention for the numerator cols are `<marker>_POS` and for
  denominator cols are `<marker>_POS_D`.  
  For both sets of columns, `0` indicates considered the negative cells,
  `1` indicates considered the positive cells and `NA_real_` indicates
  not in consideration for the subpopulation.  
  This is useful for matching to percentage data with potentially
  different naming conventions to avoid not having exact string matches
  for the same subpopulation  
  Take note that the order also matters when matching strings: "CD4+ &
  CD8- of CD3+" is different from "CD8- & CD4+ of CD3+"

## Value

tibble containing the percentage of cells where

- rows correspond to each subpopulation specified in the
  `subpopulation`,

- `n_num` indicates the number of cells that satisifies the numerator
  conditions,

- `n_denom` indicates the number of cells that satisifies the
  denominator conditions,

- `perc`=`n_num` divided by `n_denom` unless `n_denom`=0, then
  `perc=NA_real_`

## Details

The subpopulations are defined as (num marker(s)) out of (denom
marker(s)) where num denotes numerator, and denom denotes denominator
(these shorthands are used in the function arguments)

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

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

# Get percentage for CD4 based on gating
getPerc(intens_dat_2, num_marker=c("CD4"), denom_marker="CD3")
#> # A tibble: 4 × 6
#>   subpopulation      n_num n_denom  perc CD4_POS CD3_POS_D
#>   <chr>              <int>   <int> <dbl>   <dbl>     <dbl>
#> 1 CD4_NEG_OF_CD3_NEG     3      50     6       0         0
#> 2 CD4_POS_OF_CD3_NEG    47      50    94       1         0
#> 3 CD4_NEG_OF_CD3_POS    42      50    84       0         1
#> 4 CD4_POS_OF_CD3_POS     8      50    16       1         1
```
