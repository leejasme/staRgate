test_that("No errors when passing in a data.frame", {
  # Set up df
  set.seed(100)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  testthat::expect_no_error(getDensityMats(intens_dat=intens_dat, subset_col="CD3_pos", marker="CD4")
  )

})

test_that("Error when `subset_col` is not supplied.", {
  set.seed(100)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  testthat::expect_error(getDensityMats(intens_dat=intens_dat, marker="CD4"),
                         "is missing, with no default",
                         ignore.case=TRUE)

})

test_that("Error when `marker` is not supplied.", {
  set.seed(100)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  testthat::expect_error(getDensityMats(intens_dat=intens_dat),
                         "is missing, with no default",
                         ignore.case=TRUE)

})

test_that("Returns a data.frame/tibble", {
  set.seed(10)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  testthat::expect_s3_class(getDensityMats(intens_dat, subset_col="CD3_pos", marker="CD4"),
                            c("data.frame", "tbl", "tbl_df"))
})

test_that("Returned dataframe has nrows matching the n unique values from input data's `subset_col`",{
  set.seed(10)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  dens_mat = getDensityMats(intens_dat, subset_col="CD3_pos", marker="CD4")

  testthat::expect_equal(nrow(dens_mat), length(unique(intens_dat$CD3_pos)))
})

test_that("Returns `density` objects for each row of the `dens_obj` column",{
  set.seed(10)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  dens_mat = getDensityMats(intens_dat, subset_col="CD3_pos", marker="CD4")

  if("dens_obj" %in% colnames(dens_mat)){
    obj_type <- purrr::map(dens_mat$dens_obj, ~ class(.x))
  }else{
    obj_type <- FALSE
  }

  chk_obj_type <- sapply(obj_type, function(x){(x) == "density"})

  testthat::expect_equal(all(chk_obj_type), TRUE)
})

test_that("Dataframe is returned for `dens_peaks_flip` when cutting to the left of peak", {
  set.seed(100)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=500),
    # This creates a left tailed dist
    # Multiple by 2000 bc the returned random values range b/w 0-1
    # We need more variation to calculate the density!
    CD4=rbeta(1000, 5, 1)*2000,
    CD8=rnorm(1000, 100, 10)
  )

  dens_mats <- getDensityMats(intens_dat=intens_dat,
                              subset_col="CD3_pos",
                              marker="CD4",
                              bin_n=512,
                              # Set to lower ~ 1500 b/c the CD4 peak is at ~ 1900
                              pos_peak_threshold=1500)

  dens_peaks_flip_not_null <- sapply(dens_mats$dens_peaks_flip,
                                     function(x){!is.null(x)},
                                     simplify=TRUE)

  testthat::expect_equal(all(dens_peaks_flip_not_null), TRUE)

})
