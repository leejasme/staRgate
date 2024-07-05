test_that("Returns correct class", {
  set.seed(100)
  intens_dat = tibble::tibble(
                 CD3_pos = rep(c(0, 1), each = 50),
                 CD4 = rnorm(100, 100, 10),
                 CD8 = rnorm(100, 100, 10)
  )

  # Run density gating, leaving other params at suggested defaults
  # number of bins suggested is 40 but default of the function is at `bin_n = 512`,
  # which is the default for the R base density() function
  gates = staRgate::getDensityGates(intens_dat, marker = "CD4", subset_col = "CD3_pos", bin_n = 40)

  testthat::expect_s3_class(gates, c("data.frame", "tbl", "tbl_df"))
})

test_that("Returns correct number of rows", {
  set.seed(100)
  intens_dat = tibble::tibble(
    CD3_pos = rep(c(0, 1), each = 50),
    CD4 = rnorm(100, 100, 10),
    CD8 = rnorm(100, 100, 10)
  )

  # Run density gating, leaving other params at suggested defaults
  # number of bins suggested is 40 but default of the function is at `bin_n = 512`,
  # which is the default for the R base density() function
  gates = staRgate::getDensityGates(intens_dat, marker = "CD4", subset_col = "CD3_pos", bin_n = 40)

  # The nrows should be equal to the number of unique values in the CD3_pos column
  testthat::expect_equal(nrow(gates), length(unique(intens_dat$CD3_pos)))
})


test_that("Returns correct number of columns", {
  set.seed(100)
  intens_dat = tibble::tibble(
    CD3_pos = rep(c(0, 1), each = 50),
    CD4 = rnorm(100, 100, 10),
    CD8 = rnorm(100, 100, 10)
  )

  # Run density gating, leaving other params at suggested defaults
  # number of bins suggested is 40 but default of the function is at `bin_n = 512`,
  # which is the default for the R base density() function
  gates = staRgate::getDensityGates(intens_dat, marker = "CD4", subset_col = "CD3_pos", bin_n = 40)

  # The expect ncols = 2 because we only specified 1 marker (CD4) and we will have 1 column for
  # Which CD3_pos unique value/subset
  testthat::expect_equal(ncol(gates), 2)
})
