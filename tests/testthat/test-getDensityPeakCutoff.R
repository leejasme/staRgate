test_that("Returns a data.frame/tibble", {
  set.seed(10)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  dens_mat <- getDensityMats(intens_dat, subset_col="CD3_pos", marker="CD4")

  dens_derivs <- dens_mat |>
    dplyr::mutate(dens_binned=
                    purrr::map(dens_obj,
                               function(d){getDensityDerivs(d)})
    )

  cutoffs <- purrr::map(dens_derivs$dens_binned,
                        ~ getDensityPeakCutoff(.x, marker="CD4", subset_col = "CD3_pos")
  )

  chk_has_cutoff_col <- sapply(cutoffs, function(x){"cutoff" %in% colnames(x)}, simplify=TRUE)

  testthat::expect_equal(all(chk_has_cutoff_col), TRUE)
})


test_that("Returned cutoff is to the right of the identified peak", {
  set.seed(10)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  dens_mat <- getDensityMats(intens_dat, subset_col="CD3_pos", marker="CD4")

  dens_derivs <- dens_mat |>
    dplyr::mutate(dens_binned=
                    purrr::map(dens_obj,
                               function(d){getDensityDerivs(d)})
    )

  cutoffs <- purrr::map(dens_derivs$dens_binned,
                        ~ getDensityPeakCutoff(.x, marker="CD4", subset_col = "CD3_pos")
  )

  location_of_peak <- sapply(cutoffs,
                       function(x){
                         dplyr::pull(dplyr::slice_min(dplyr::filter(x, peak == TRUE), original_row_num), original_row_num)
                       })

  location_of_cut <-
    sapply(cutoffs,
           function(x){
             dplyr::pull(dplyr::slice_min(dplyr::filter(x, cutoff == TRUE), original_row_num), original_row_num)
           })

  testthat::expect_equal(all(location_of_peak < location_of_cut), TRUE)
})


test_that("Returned cutoff is to the left of the identified peak for a flipped case", {

  set.seed(10)
  intens_dat <- data.frame(
    CD3_pos=rep(c(0, 1), each=50),
    CD4=rnorm(100, 100, 10),
    CD8=rnorm(100, 100, 10)
  )

  dens_mat = getDensityMats(intens_dat, subset_col="CD3_pos", marker="CD4")

  dens_derivs = dens_mat |>
    dplyr::mutate(dens_binned=
                    purrr::map(dens_obj,
                               function(d){getDensityDerivs(d)})
    )

  cutoffs = purrr::map(dens_derivs$dens_binned,
                       ~ getDensityPeakCutoff(.x, marker="CD4", subset_col="CD3_pos", dens_flip=TRUE)
  )

  location_of_peak <-
    sapply(cutoffs,
           function(x){
             dplyr::pull(dplyr::slice_min(dplyr::filter(x, peak == TRUE), original_row_num), original_row_num)
           })

  location_of_cut <-
    sapply(cutoffs,
           function(x){
             dplyr::pull(dplyr::slice_min(dplyr::filter(x, cutoff == TRUE), original_row_num), original_row_num)
           })

  testthat::expect_equal(all(location_of_peak > location_of_cut), TRUE)
})
