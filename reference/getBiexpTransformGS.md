# Applies Biexpeonential Transformation using specifications in csv file provided at `path_biexp_params`

The csv file at `path_biexp_params` should specify the channels to apply
the transformation to and the parameters (negative decades, width basis
and positive decades). The default is negative decades=0.5, width
basis=-30 and positive decades=4.5. The Transformation can be applied to
only a subset of the channels included in the GatingSet.

## Usage

``` r
getBiexpTransformGS(gs, path_biexp_params)
```

## Arguments

- gs:

  GatingSet to apply Biexponential Transformation to

- path_biexp_params:

  file path for .csv file that specifies the Biexponential
  Transformation

## Value

GatingSet with Biexponentially Transformed data

## Details

An example table is provided in the
extdata/biexp_transf_parameters_x50.csv

## Examples

``` r
# This example does not contain all the pre-processing steps required in
# getting the GatingSet (gs) ready for Biexp transformation.
# To see the steps that are required to creating the (gs),
# please see the vignette for a full tutorial

# To make this a runnable example, read in the FCS file to create gs and
# directly apply

# File path to the FCS file
path_fcs <- system.file("extdata",
                        "example_fcs.fcs",
                        package="staRgate",
                        mustWork=TRUE)
path_biexp_params <- system.file("extdata",
                                 "biexp_transf_parameters_x50.csv",
                                 package="staRgate",
                                 mustWork=TRUE)

# Create a cytoset then convert to gs
cs <- flowWorkspace::load_cytoset_from_fcs(path_fcs)
gs <- flowWorkspace::GatingSet(cs)

# gs must be a GatingSet object
gs <- getBiexpTransformGS(gs, path_biexp_params=path_biexp_params)

# To check the transformation parameters applied
flowWorkspace::gh_get_transformations(gs)
#> $`AF700-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dd5118>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`APC-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dd3b38>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`APC-f750-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dd6388>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BB515-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dd8f78>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BB660-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dd7998>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BB700-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dda1e8>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BB790-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dde958>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BUV395-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8de12d8>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BUV496-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8ddfcf8>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BUV563-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8de46d8>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BUV615-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8de6f28>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BUV661-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8de5948>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BUV737-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8de82c8>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BUV805-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8deab18>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BV421-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8de9538>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BV480-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8debd88>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BV510-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dee530>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BV570-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8df0f10>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BV605-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8def930>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BV650-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8df4840>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BV711-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8df7090>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BV750-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8df5ab0>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`BV786-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8df8300>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`PE-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dfab50>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`PE-CF594-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dfd3a0>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`PE-Cy5-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dfbdc0>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`PE-Cy5.5-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8dfe610>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
#> $`PE-Cy7-A`
#> function (x, deriv = 0) 
#> {
#>     deriv <- as.integer(deriv)
#>     if (deriv < 0 || deriv > 3) 
#>         stop("'deriv' must be between 0 and 3")
#>     if (deriv > 0) {
#>         z0 <- double(z$n)
#>         z[c("y", "b", "c")] <- switch(deriv, list(y = z$b, b = 2 * 
#>             z$c, c = 3 * z$d), list(y = 2 * z$c, b = 6 * z$d, 
#>             c = z0), list(y = 6 * z$d, b = z0, c = z0))
#>         z[["d"]] <- z0
#>     }
#>     res <- stats:::.splinefun(x, z)
#>     if (deriv > 0 && z$method == 2 && any(ind <- x <= z$x[1L])) 
#>         res[ind] <- ifelse(deriv == 1, z$y[1L], 0)
#>     res
#> }
#> <bytecode: 0x55e7d8419000>
#> <environment: 0x55e7d8e00e60>
#> attr(,"type")
#> [1] "biexp"
#> attr(,"parameters")
#> attr(,"parameters")$channelRange
#> [1] 4096
#> 
#> attr(,"parameters")$maxValue
#> [1] 262144
#> 
#> attr(,"parameters")$neg
#> [1] 0.5
#> 
#> attr(,"parameters")$pos
#> [1] 4.5
#> 
#> attr(,"parameters")$widthBasis
#> [1] -30
#> 
#> 
```
