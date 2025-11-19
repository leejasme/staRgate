# Internal function: Estimate derivatives for density of `marker` for each unique subset of `subset_col`

Internal function for `get_density_gates` For each unique value in
`subset_col`, estimate the derivatives for `marker` (intensity values)

## Usage

``` r
getDensityDerivs(dens)
```

## Arguments

- dens:

  `density` object from the
  [density](https://rdrr.io/r/stats/density.html)

## Value

list of dataframe with density estimation and corresponding 1st-4th
derivatives, indicators of local peaks, plateau_pre  
each element corresponds to each unique value of `subset_col`  
for each dataframe: rows correspond to each of the bins
