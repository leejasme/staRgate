test_that("Returns an error if not passing in a density object", {
  set.seed(10)
  # Make a data.frame instead of a density object to pass in
  df <- data.frame(x=rnorm(100),
                  y=rnorm(100))


  testthat::expect_error(getDensityDerivs(df),
                         "is not a density object",
                         ignore.case=TRUE
  )
})

test_that("No error when passing in density object", {
  set.seed(10)
  # Make a data.frame
  df <- data.frame(x=rnorm(100),
                   y=rnorm(100))

  # Estimate the density of x
  density_x = stats::density(df$x)

  testthat::expect_no_error(getDensityDerivs(density_x))
})

test_that("Function returns a tibble", {
  set.seed(10)
  # Make a data.frame
  df <- data.frame(x=rnorm(100),
                   y=rnorm(100))

  # Estimate the density of x
  density_x <- stats::density(df$x)

  testthat::expect_s3_class(getDensityDerivs(density_x), c("data.frame", "tbl", "tbl_df"))
})

test_that("Function returns correct number of columns", {
  set.seed(10)
  # Make a data.frame
  df <- data.frame(x=rnorm(100),
                   y=rnorm(100))

  # Estimate the density of x
  density_x <- stats::density(df$x)

  # Calculate the dataframe
  calc_df <- getDensityDerivs(density_x)

  testthat::expect_equal(ncol(calc_df), 12)
})

test_that("Returns the first to fourth derivative columns and local peak", {
  set.seed(10)
  # Make a data.frame
  df <- data.frame(x=rnorm(100),
                   y=rnorm(100))

  # Estimate the density of x
  density_x <- stats::density(df$x)

  # Calculate the dataframe
  dens_derivs <- getDensityDerivs(density_x)

  chk_colnames = c("first_deriv", "second_deriv", "third_deriv", "fourth_deriv", "local_peak") %in% colnames(dens_derivs)

  testthat::expect_equal(all(chk_colnames), TRUE)
})

test_that("Returns a logical for local peak", {
  set.seed(10)
  # Make a data.frame
  df <- data.frame(x=rnorm(100),
                   y=rnorm(100))

  # Estimate the density of x
  density_x <- stats::density(df$x)

  # Calculate the dataframe
  dens_derivs <- getDensityDerivs(density_x)

  chk_type <- class(dens_derivs$local_peak)

  testthat::expect_equal(chk_type, "logical")
})
