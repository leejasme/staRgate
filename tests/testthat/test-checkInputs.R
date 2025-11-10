test_that("Returns error when matrix passed as input", {
  # Pass in a matrix instead of expected data.frame
  set.seed(100)
  intens_dat <- matrix(
    c(rep(c(0, 1), each=50), rnorm(100, 100, 10), rnorm(100, 100, 10)),
    ncol=3,
    byrow=FALSE,
    dimnames = list(NULL, c("CD3_pos", "CD4", "CD8"))
  )

  testthat::expect_error(checkInputs(intens_dat),
                         "must be of data.frame class",
                         ignore.case=TRUE
  )
})

test_that("Returns warning if not numeric values passed in", {
  # Pass in char for all expected numeric inputs
  testthat::expect_warning(checkInputs(bin_n="512"),
                           "Warning",
                           ignore.case=TRUE
  )
})

test_that("Returns warning if not numeric values passed in", {
  # Pass in char for all expected numeric inputs
  testthat::expect_warning(checkInputs(peak_detect_ratio=FALSE),
                           "Warning",
                           ignore.case=TRUE
  )
})

test_that("Returns warning if not numeric values passed in", {
  # Pass in char for all expected numeric inputs
  testthat::expect_warning(checkInputs(neg_intensity_threshold=TRUE),
                           "Warning",
                           ignore.case=TRUE
  )
})

test_that("Returns warning if not numeric values passed in", {
  # Pass in char for all expected numeric inputs
  testthat::expect_warning(checkInputs(neg_intensity_threshold=TRUE),
                           "Warning",
                           ignore.case=TRUE
  )
})

test_that("Returns error if column names specified do not exist in data",{
  # Pass in a matrix instead of expected data.frame
  set.seed(100)
  intens_dat <- tibble::tibble(
    CD3_pos = rep(c(0, 1), each = 50),
    CD4 = rnorm(100, 100, 10),
    CD8 = rnorm(100, 100, 10)
  )

  testthat::expect_error(checkInputs(intens_dat=intens_dat,
                                     marker="CD45RA"),
                         "must be string matching column name",
                         ignore.case=TRUE
  )
})

test_that("Returns error if column names are of incorrect class",{
  # Pass in a matrix instead of expected data.frame
  set.seed(100)
  intens_dat <- tibble::tibble(
    CD3_pos = rep(c(0, 1), each = 50),
    CD4 = rnorm(100, 100, 10),
    CD8 = rnorm(100, 100, 10)
  )

  testthat::expect_error(checkInputs(intens_dat=intens_dat,
                                     subset_col=450),
                         "must be string matching column name",
                         ignore.case=TRUE
  )
})

test_that("Returns error if column names are of incorrect class",{
  # Pass in a matrix instead of expected data.frame
  set.seed(100)
  intens_dat <- tibble::tibble(
    CD3_pos = rep(c(0, 1), each = 50),
    CD4 = rnorm(100, 100, 10),
    CD8 = rnorm(100, 100, 10)
  )

  testthat::expect_error(checkInputs(intens_dat=intens_dat,
                                     subset_col=450),
                         "must be string matching column name",
                         ignore.case=TRUE
  )
})

test_that("Returns warning for incorrect positive peak threshold data input",{
  ppt =
    tibble::tibble(
      marker = c("KI67"),
      threshold = 1800
    )

  testthat::expect_error(checkInputs(pos_peak_threshold=ppt),
                         "Column names for \\`pos_peak_threshold\\` should be \\`marker\\` and \\`pos_peak_threshold\\`",
                         ignore.case=TRUE
  )
})
