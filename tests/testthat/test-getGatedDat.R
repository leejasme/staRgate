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

  # Tag on the 0/1 on intens_dat
  intens_dat_2 = staRgate::getGatedDat(intens_dat, cutoffs = gates, subset_col = "CD3_pos")

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

  # Tag on the 0/1 on intens_dat
  intens_dat_2 = staRgate::getGatedDat(intens_dat, cutoffs = gates, subset_col = "CD3_pos")

  # The nrows should be equal to that in intens_dat
  testthat::expect_equal(nrow(intens_dat_2), nrow(intens_dat))
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

  # Tag on the 0/1 on intens_dat
  intens_dat_2 = staRgate::getGatedDat(intens_dat, cutoffs = gates, subset_col = "CD3_pos")

  # The expect ncols = 1 additional column (for CD4 pos/neg)
  testthat::expect_equal(ncol(intens_dat_2), 1 + ncol(intens_dat))
})

test_that("The values in the new column are only 0 or 1s", {
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

  # Tag on the 0/1 on intens_dat
  intens_dat_2 = staRgate::getGatedDat(intens_dat, cutoffs = gates, subset_col = "CD3_pos")

  # The values in cd4_pos should only contain 0 or 1s
  testthat::expect_contains(unique(intens_dat_2$cd4_pos), c(0, 1))
})
