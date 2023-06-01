get_density_derivs = function(dens,
                              marker,
                              subset_col,
                              bin_n = 512,
                              peak_detect_ratio = 10,
                              pos_peak_threshold = 1600){

  #' Internal function: Estimate derivatives for density of `marker` for each unique subset of `subset_col`
  #'
  #' Internal function for `get_density_gates`
  #' For each unique value in `subset_col`, estimate the derivatives for
  #' `marker` (intensity values)
  #'
  #'
  #' @param dens `density` object from the \link[stats]{density}
  #' @param marker string for the marker to gate on
  #'              the name needs to match exactly the column name in `intens_dat`
  #' @param subset_col string for the column name to indicate the subsets to apply density gating on
  #'              will perform operation on subsets corresponding to each unique value in column
  #' @param bin_n numeric to be passed to `n` parameter of `density(n = bin_n)` for
  #'              number of equally spaced points at which the density is to be estimated
  #'              default is 512, which is the default of `density(n = 512)`
  #' @param peak_detect_ratio numeric threshold for eliminating small peaks where
  #'              a peak that is < than the highest peak by `peak_detect_ratio` times will be ignored
  #'              default = 10
  #' @param pos_peak_threshold numeric for threshold to identify a positive peak
  #'           '  default is 1600, which is on the biexponential scale
  #'
  #' @return list of dataframe with density estimation and corresponding 1st-4th derivatives,
  #'         indicators of local peaks, plateau_pre
  #'         each element corresponds to each unique value of `subset_col`
  #'         for each dataframe: rows correspond to each of the bins
  #' @export
  #'

  # Meant for internally calling within get_density_peaks and for debug/checking matrices/calculations
  # input data is diff- expecting the density obj for specific cd4/cd8 subset
  x = dens$x

  # x every 1
  x_binned =
    x %>%
    range() %>%
    {seq(.[[1]], .[[2]], by = 1)}

  # Find interval?
  new_d =
    # dplyr depends on tibble so should be ok to do
    tibble::tibble(
      x = x,
      y = dens$y,
      x_int = findInterval(x, x_binned)
    ) %>%
    # If not use interval of every 100 values
    # group_by(x) %>%
    dplyr::group_by(x_int) %>%
    dplyr::summarize(
      x_avg = mean(x),
      y_avg = mean(y)
    ) %>%
    dplyr::arrange(x_avg) %>%
    # 2022-09-06 add the original row num here
    # Want the original row num for merging back later
    # Based on the x continuous value might be problematic with rounding
    tibble::rownames_to_column(var = "original_row_num") %>%
    dplyr::mutate(original_row_num = as.numeric(original_row_num)) %>%
    # purrr::when(
    #   dens_flip == TRUE ~
    #     mutate(.,
    #            x_avg_flip = -1*x_avg),
    #   # When no flipping is needed, x_flip = x
    #   # Then we can arrange by x_flip before taking derivs
    #   ~
    #     mutate(.,
    #            x_avg_flip = x_avg)
    # ) %>%
    # arrange(x_avg_flip) %>%
  # 2022-03-25 may not need since we check ratio of peaks
  # these Fhould also be teased out
  # dplyr::filter(x_avg < aux_range[2], x_avg > aux_range[1]) %>%
  dplyr::mutate(
    # 2022-12-30: Changed to diff(y)/diff(x) instead of just diff(y) but
    # should not change shape of the derivs b/c x is equidistance
    # only changes the actual values
    first_deriv = c(0, diff(y_avg)/diff(x_avg)),
    first_deriv_sign = sign(first_deriv),
    second_deriv = c(0, diff(first_deriv)/diff(x_avg)),
    second_deriv_sign = sign(second_deriv),
    # 2022-06-09 add third deriv
    third_deriv = c(0, diff(second_deriv)/diff(x_avg)),
    third_deriv_sign = sign(third_deriv),
    fourth_deriv = c(0, diff(third_deriv)/diff(x_avg)),
    # Peak if change from + to - for first deriv (-1 - 1 = -2)
    # second derive < 0 = peak
    local_peak = (c(0, diff(first_deriv_sign)) == -2) & second_deriv_sign < 0
  )

  # Simpliest fix to shift local_peak 1 index "up"
  new_d$local_peak =
    c(new_d$local_peak[-1], FALSE)

  # 2022-06-15 split out the plateau_pre step bc need to fix the off-by-1 for local_peak identified?
  new_d =
    new_d %>%
    dplyr::mutate(
      # 2022-04-14 add condition for 2nd deriv positive
      # 2022-03-17 too early if 2nd deriv changes from neg to pos
      # plateau ~ slope changes from very neg to less neg
      # sign for 2nd deriv of less negative slope = + 1
      # second_deriv should be changing from neg to pos
      # plateau_pre = (c(0, diff(second_deriv)) >= 0) & second_deriv_sign > 0
      # 2022-03-17: want it to be farther out where the 2nd deriv "stabilizes"
      # which is the first peak in positive region after the negative peak in density identified
      # Identify peak of second_deriv with diff(sign(diff(second_deriv))) == -2 when slope
      # changes from pos to neg, need the extra 0 because of the outer diff() required
      # 2022-06-09 try 2nd deriv local max instead of the first dip in 2nd deriv with 2 consecutive points in the positive region of 2nd deriv
      # plateau_pre = ((c(0, diff(second_deriv_sign)) == 0) & (c(0, 0, diff(sign(diff(second_deriv)))) == -2) &
      #                  second_deriv > 0)
      # 2022-06-09: replaces above with identifing 2nd deriv local max in positive region of 2nd deriv
      # 2022-07-11: the special off case where pt follow peak in 2nd deriv is in second_deriv < 0
      # We should check
      # plateau_pre = ((c(0, diff(third_deriv_sign)) == -2) & sign(fourth_deriv) < 0), # & second_deriv_sign > 0),
      # plateau_pre_2 = sapply(row_number(),
      #                        function(r){if(r >= 2){(plateau_pre[r] == TRUE & second_deriv_sign[(r-1)] > 0)}else{FALSE}})

      # Due to how the diff and deriv are calculated, the point of interest is back tracked by 1 so check 2nd deriv of 1 point before
      plateau_pre = ((c(0, diff(third_deriv_sign)) == -2) & (sign(fourth_deriv)< 0)),
      plateau_pre_2 = sapply(dplyr::row_number(),
                             function(r){if(r >= 2){(plateau_pre[r] == TRUE & second_deriv_sign[(r-1)] > 0)}else{FALSE}})
    )

  return(new_d)

}
