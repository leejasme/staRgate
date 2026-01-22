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

test_that("No errors when passing in a data.frame", {
  # Set up df
  set.seed(100)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  testthat::expect_no_error(checkInputs(intens_dat=intens_dat)
  )

})

test_that("Returns warning if not numeric values passed in for `bin_n`", {
  # Pass in char for all expected numeric inputs
  testthat::expect_warning(checkInputs(bin_n="512"),
                           "not specified as numeric",
                           ignore.case=TRUE
  )
})

test_that("Returns warning if not numeric values passed in for `peak_detect_ratio`", {
  # Pass in char for all expected numeric inputs
  testthat::expect_warning(checkInputs(peak_detect_ratio=FALSE),
                           "not specified as numeric",
                           ignore.case=TRUE
  )
})

test_that("Returns warning if not numeric values passed in for `neg_intensity_threshold`", {
  # Pass in logical for all expected numeric inputs
  testthat::expect_warning(checkInputs(neg_intensity_threshold=TRUE),
                           "must be a numeric",
                           ignore.case=TRUE
  )
})

test_that("Returns error if marker names specified do not exist in data",{
  # Set up df
  set.seed(100)
  intens_dat <- tibble::tibble(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  testthat::expect_error(checkInputs(intens_dat=intens_dat,
                                     marker="CD45RA"),
                         "must be string matching column name",
                         ignore.case=TRUE
  )
})

test_that("Returns error if column names are of incorrect class",{
  # Set up df
  set.seed(100)
  intens_dat <- tibble::tibble(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  testthat::expect_error(checkInputs(intens_dat=intens_dat,
                                     subset_col=450),
                         "must be string matching column name",
                         ignore.case=TRUE
  )
})

test_that("Returns warning for incorrect positive peak threshold data input",{
  ppt <-
    tibble::tibble(
      marker=c("KI67"),
      threshold=1800
    )

  testthat::expect_error(checkInputs(pos_peak_threshold=ppt),
                         "Column names for \\`pos_peak_threshold\\` should be \\`marker\\` and \\`pos_peak_threshold\\`",
                         ignore.case=TRUE
  )
})

test_that("All provided positive peak thresholds are within the input data",{
  ppt <-
    tibble::tibble(
      marker=c("KI67"),
      threshold=1800
    )

  # Set up df
  set.seed(100)
  intens_dat <- tibble::tibble(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10),
    EOMES=rnorm(100, 100, 10)
  )

  testthat::expect_error(checkInputs(intens_dat=intens_dat,
                                     marker="EOMES",
                                     pos_peak_threshold=ppt),
                         "Not all markers in \\`marker\\` argument have a \\`pos_peak_threshold\\` supplied",
                         ignore.case=TRUE
  )
})

test_that("Returns error if the markers of interest for calculating percentages are not gated.", {
  # Set up df
  set.seed(100)
  intens_dat <- tibble::tibble(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10),
    EOMES=rnorm(100, 100, 10)
  )

  testthat::expect_error(checkInputs(intens_dat=intens_dat,
                                     num_marker="EOMES",
                                     denom_marker="CD4"),
                         "\\`num_marker\\` and \\`denom_marker\\` must have indicator columns in `intens_dat`",
                         ignore.case=TRUE
  )
})

test_that("Returns error if the gated columns contain values other than 0 and/or 1.", {
  # Set up df
  set.seed(100)
  intens_dat <- tibble::tibble(
    CD3_POS=rep(c(0, 1), each=50),
    CD4_POS=rep(c(0, 1), each=50),
    CD8=rnorm(100, 100, 10),
    EOMES_POS=rep(c(0, 1, 2, 6, NA_real_), each=20),
  )

  testthat::expect_error(checkInputs(intens_dat=intens_dat,
                                     num_marker="EOMES",
                                     denom_marker="CD4"),
                         "Unique values in indicator columns corresponding to \\`num_marker\\` and \\`denom_marker\\` should only contain 0 and 1",
                         ignore.case=TRUE
  )
})

test_that("Returns error if the gated columns contain values of 0 or 1 but are characters instead of numeric.", {
  # Set up df
  set.seed(100)
  intens_dat <- tibble::tibble(
    CD3_POS=rep(c(0, 1), each=50),
    CD4_POS=as.character(rep(c(0, 1), each=50)),
    CD8=rnorm(100, 100, 10)
  )

  testthat::expect_error(checkInputs(intens_dat=intens_dat,
                                     num_marker="CD4",
                                     denom_marker="CD3"),
                         "Not all indicator columns corresponding to \\`num_marker\\` and \\`denom_marker\\` are numeric",
                         ignore.case=TRUE
  )
})

test_that("Returns warning if the `expand_num` or `expand_denom` are not logicals (TRUE/FALSE).",{

  testthat::expect_warning(checkInputs(expand_num="500",
                                       expand_denom=FALSE),
                         "not of class logical",
                         ignore.case=TRUE
  )

  testthat::expect_warning(checkInputs(expand_denom=TRUE,
                                       expand_num=100),
                         "not of class logical",
                         ignore.case=TRUE
  )

  testthat::expect_warning(checkInputs(expand_denom="TRUE",
                                       expand_num=FALSE),
                         "not of class logical",
                         ignore.case=TRUE
  )

  testthat::expect_warning(checkInputs(expand_denom="TRUE",
                                       expand_num="FALSE"),
                           "not of class logical",
                           ignore.case=TRUE
  )
})
