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

  # Get percentage for CD4 based on gating
  tbl_perc = staRgate::getPerc(intens_dat_2, num_marker = c("CD4"), denom_marker = "CD3")

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

  # Get percentage for CD4 based on gating
  tbl_perc = staRgate::getPerc(intens_dat_2, num_marker = c("CD4"), denom_marker = "CD3")

  # The nrows should be 4 (2 possible CD4 values and 2 possible CD3 values)
  testthat::expect_equal(nrow(tbl_perc), 4)
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

  # Get percentage for CD4 based on gating
  tbl_perc = staRgate::getPerc(intens_dat_2, num_marker = c("CD4"), denom_marker = "CD3")

  # The expect ncols = 6
  # subpop, n_num, n_denom, perc, the indicators for num CD4 and denom CD3
  testthat::expect_equal(ncol(tbl_perc), 6)
})

test_that("The values in the percent column are within 0-100%", {
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

  # Get percentage for CD4 based on gating
  tbl_perc = staRgate::getPerc(intens_dat_2, num_marker = c("CD4"), denom_marker = "CD3")

  # The values in perc should all be within 0-100
  testthat::expect_true(all(tbl_perc$perc <= 100 & tbl_perc$perc >= 0))
})
