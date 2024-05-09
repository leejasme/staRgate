getDensityDerivs = function(dens,
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
  #'
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


getDensityMats = function(intens_dat,
                          marker,
                          subset_col,
                          bin_n = 512,
                          peak_detect_ratio = 10,
                          pos_peak_threshold = 1800){
  #' Internal function: Matrix of calculations for density gating of intensity values in `marker` for each unique subset of `subset_col`
  #'
  #' Internal function for `getDensityGates`
  #' For each unique value in `subset_col`, there is a matrix for storing calculations for density gating
  #' contains: first to fourth derivatives of density,
  #' indicators for local peaks, "real peaks", plateau_pre and cutoff
  #'
  #'
  #' @param intens_dat dataframe of pre-gated (compensated, biexp. transf, gated CD4/CD8) intensity values where
  #'              cols = intensity value per marker,
  #'              rows = each sample
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
  #'           '  default is 1800, which is on the biexponential scale
  #'
  #' @return tibble of matrices for `marker` containing calculations for density gating
  #'         for each unique subset found in `subset_col`
  #'         rows correspond to unique values in `subset_col`, cols correspond to ?
  #'
  # Meant for using internally in getDensityGates but also for debug
  # to grab the full matrix of density, peak and cutoff
  # rows is each bin for the density estimation
  # cols for x (intens), y (dens) values, derivs (1-4th), sign of derivs,
  # indicators for local peak, plateau_pre (for calculating cutoff)
  # indicator of "real peaks" and cutoff
  # Same args as getDensityGates


  ## Calc density ---
  # Need to separate out per cd4/cd8 subsets
  dens =
    intens_dat %>%
    dplyr::group_by(!!dplyr::sym(subset_col)) %>%
    # each subset's data is separated out in `data`
    tidyr::nest() %>%
    # calculate density for each subset separately
    dplyr::mutate(
      dens_obj = purrr::map(data,
                            ~ stats::density(.x[[marker]],
                                             n = bin_n))
    ) %>%
    dplyr::ungroup()

  # Smooth the density by taking bins of `bin_size`,
  # calc 1st and 2nd deriv of each bin's avg y (density values) at avg x value
  # identify local peaks and plateuing
  dens_binned =
    dens %>%
    dplyr::mutate(
      dens_binned =
        purrr::map(dens_obj,
                   function(d){
                     getDensityDerivs(d,
                                      # dens_flip = FALSE,
                                      subset_col = subset_col,
                                      marker = marker,
                                      # bin_size = bin_size,
                                      bin_n = bin_n,
                                      peak_detect_ratio = peak_detect_ratio,
                                      pos_peak_threshold = pos_peak_threshold
                     )
                   }
        ),
      dens_peaks =
        purrr::map(
          dens_binned,
          function(d){
            getDensityPeakCutoff(d,
                                 marker = marker,
                                 subset_col = subset_col,
                                 dens_flip = FALSE,
                                 # bin_size = bin_size,
                                 bin_n = bin_n,
                                 peak_detect_ratio = peak_detect_ratio,
                                 pos_peak_threshold = pos_peak_threshold)
          }
        ),
      # 2022-11-01 check how many peaks
      # if > 2, reduce bin width before checking for desnity flip
      n_peaks =
        purrr::map(
          dens_peaks,
          function(d){
            sum(d$peak)
          }
        ) %>%
        unlist(),
      # bin_width = bin_size,
      bin_n = bin_n
    )

  # # bin_size_multipeaks = 50
  # # Keep the same bin_width = 1 for multipeaks
  # bin_size_multipeaks = 1
  #
  # for(i in which(dens_binned$n_peaks >= 2)){
  #   # 2022-12-05 use the bin_width variable and update bin width
  #   dens_binned$bin_width[[i]] = bin_size_multipeaks
  #
  #   # Update w smaller binsize
  #   dens_binned$dens_binned[[i]] =
  #     getDensityDerivs(dens_binned$dens_obj[[i]], dens_flip = FALSE,
  #                        marker = marker,
  #                        bin_size = bin_size_multipeaks,
  #                        peak_detect_ratio = peak_detect_ratio,
  #                        pos_peak_threshold = pos_peak_threshold
  #     )
  #
  #
  #   dens_binned$dens_peaks[[i]] =
  #     getDensityPeakCutoff(dens_binned$dens_binned[[i]], marker = marker, dens_flip = FALSE,
  #                             bin_size = bin_size_multipeaks,
  #                             peak_detect_ratio = peak_detect_ratio,
  #                             pos_peak_threshold = pos_peak_threshold)
  #
  #   #I dont think we want to update this for the more refined density?
  #   # dens_binned$n_peaks[[i]]
  # }

  # Check the peak value against threshold
  flag_peak =
    dens_binned %>%
    dplyr::select(!!dplyr::sym(subset_col), dens_peaks) %>%
    dplyr::mutate(
      peak_loc =
        purrr::map(dens_peaks,
                   function(x){
                     temp = x %>%
                       dplyr::filter(peak == TRUE) %>%
                       dplyr::slice_min(x_avg) %>%
                       dplyr::pull(x_avg)

                     # TO get around identifying 0 peaks
                     if(length(temp) == 0){
                       return(NA_real_)
                     }else{
                       return(temp)
                     }
                   }) %>%
        unlist(),
      flag_pos_peak =
        peak_loc > pos_peak_threshold
    ) %>%
    dplyr::select(-dens_peaks) %>%
    dplyr::filter(flag_pos_peak == TRUE)

  # Recalculate with flipped binned intensity
  # Chose to not use flipped raw intensity bc
  # worried about binning on flipped intensity will results
  # in diff bins and not a 1-1 match on flipped bins when
  # replacing cutoff from the first attempt
  # To prevent unexpected errors with empty rows
  if(nrow(flag_peak) > 0){
    dens_flipped =
      dens_binned %>%
      # Only filter to the subsets that needed to be recalculated
      # dplyr::filter(cd4_pos_cd8_pos %in% flag_peak$cd4_pos_cd8_pos) %>%
      # 2023-04-27 try to replace the explicit column name cd4_pos_cd8_pos w
      # subset_col
      dplyr::filter(
        # need a glue::glue to work
        !!rlang::parse_expr(glue::glue("{subset_col} %in% flag_peak[[subset_col]]"))
        # cd4_pos_cd8_pos %in% flag_peak$cd4_pos_cd8_pos
      ) %>%
      dplyr::mutate(
        # 2023-03-30 don't want to recalcualte density
        # # 2022-12-05 add the bin_width from flag_peak
        # # bin_width = flag_peak$bin_width,
        # dens_binned_flip =
        #   purrr::map2(dens_obj, bin_width,
        #      function(d, w){
        #        getDensityDerivs(d,
        #                           dens_flip = TRUE,
        #                           marker = marker,
        #                           bin_size = w,
        #                           bin_n = bin_n,
        #                           peak_detect_ratio = peak_detect_ratio,
        #                           pos_peak_threshold = pos_peak_threshold)
        #      }
        #      ),
        dens_peaks_flip =
          purrr::map(
            # 2023-03-30 dont use a recalculated density with --1*x_avg
            # Use original density but with dens_flip = TRUE arg to search on left vs. right
            dens_binned,
            # dens_binned_flip,
            # bin_width,
            function(d){
              getDensityPeakCutoff(d,
                                   marker = marker,
                                   subset_col = subset_col,
                                   dens_flip = TRUE,
                                   # bin_size = w,
                                   bin_n = bin_n,
                                   peak_detect_ratio = peak_detect_ratio,
                                   pos_peak_threshold = pos_peak_threshold)
            }
          )
      ) %>%
      dplyr::select(dplyr::all_of(subset_col),
                    dens_peaks_flip)
  }else if(nrow(flag_peak) == 0){
    dens_flipped =
      dens %>%
      dplyr::select(!!dplyr::sym(subset_col)) %>%
      dplyr::mutate(dens_peaks_flip =
                      rep(list(NULL), length(unique(dens[, subset_col]))))
  }

  # Replace cutoff col only for the subsets in dens_flipped
  # Need to merge the flipped data and use maps to replace?
  dens_final =
    dplyr::left_join(
      dens_binned,
      dens_flipped,
      by = subset_col
    ) %>%
    dplyr::mutate(
      dens_peaks_final =
        if(all(sapply(dens_peaks_flip, is.null, simplify = TRUE))){
          dens_peaks
        }else{
          purrr::pmap(
            list(
              dens_og = dens_peaks,
              dens_f = dens_peaks_flip
            ),
            function(dens_og, dens_f){
              if(!is.null(dens_f)){
                dplyr::left_join(
                  dens_og %>%
                    dplyr::select(-cutoff),
                  dens_f %>%
                    dplyr::select(original_row_num,
                                  cutoff),
                  by = "original_row_num"
                )
              }else{
                dens_og
              }
            }
          )
        }
    ) %>%
    # Try to remove unnecessary parts
    dplyr::select(-dens_binned)

  # return(dens_binned)
  return(dens_final)
}


getDensityPeakCutoff = function(dens_binned_dat,
                                marker,
                                subset_col,
                                bin_n = 512,
                                peak_detect_ratio = 10,
                                pos_peak_threshold = 1600,
                                dens_flip = FALSE){

  #' Internal function: Determine the "real peaks" and cutoff
  #'
  #' Internal function for `getDensityGates`
  #'
  #' @param dens_binned_dat list of dataframe output from the `getDensityDerivs`
  #' @param marker string for the marker to gate on
  #'              the name needs to match exactly the column name in `dens_binned_dat`
  #' @param subset_col string for the column name to indicate the subsets to apply density gating on
  #'              will perform operation on subsets corresponding to each unique value in column
  #' @param bin_n numeric to be passed to `n` parameter of `density(n = bin_n)` for
  #'              number of equally spaced points at which the density is to be estimated
  #'              default is 512, which is the default of `density(n = 512)`
  #' @param peak_detect_ratio numeric threshold for eliminating small peaks where
  #'              a peak that is < than the highest peak by `peak_detect_ratio` times will be ignored
  #'              default = 10
  #' @param pos_peak_threshold numeric for threshold to identify a positive peak
  #'              default is 1600, which is on the biexponential scale
  #' @param dens_flip logical for whether the gating should be applied "backwards" where the peak is
  #'              a positive peak and want to gate to the left of peak instead of right
  #'
  #' @return list of dataframe `dens_binned_dat` with additional columns added for
  #'         peak(s) identified and the cutoff
  #'         each element corresponds to each unique value of `subset_col`
  #'         for each dataframe: rows correspond to each of the bins
  #'
  #'

  # Identify "real peaks"
  new_d =
    dens_binned_dat %>%
    dplyr::filter(local_peak == TRUE) %>%
    dplyr::arrange(-y_avg) %>%
    # Get ratios relative to the highest peak
    dplyr::mutate(
      ratio = sapply(dplyr::row_number(),
                     function(r){y_avg[r]/y_avg[1]}
      )
    ) %>%
    # Remove any peaks that are "too small"
    # threshold is at ratio < 1/peak_detect_ratio is "too small"
    # If the highest peak is > peak_detect_ratio times the height of the peak in question
    # then the peak is too small
    dplyr::filter(ratio >= 1/peak_detect_ratio) %>%
    # denoting these as "true peaks"
    dplyr::mutate(peak = TRUE)  %>%
    dplyr::select(original_row_num, peak) %>%
    dplyr::left_join(
      dens_binned_dat,
      .,
      by = "original_row_num"
    ) %>%
    # Replace the NA in the peak col to avoid confusion
    dplyr::mutate(peak = #replace_na(peak, FALSE))
                    # Use ifelse instead to not rely on tidyr::replace_na
                    ifelse(is.na(peak), FALSE, peak)
    )

  # 2022-04-14 for all peaks need to correct off-by-one error where "peak" is identified at the
  #   # 1 index over from actual peak
  aux_ind_all_peaks =
    new_d %>%
    dplyr::filter(peak == TRUE) %>%
    # 2022-12-30 arrange by the height
    dplyr::arrange(-y_avg) %>%
    dplyr::pull(original_row_num)

  new_d =
    new_d %>%
    dplyr::mutate(peak = dplyr::case_when(
      # set the aux_ind_neg_peak as peak = TRUE to correct for off-by-one error
      # original_row_num %in% (aux_ind_all_peaks - 1) ~ TRUE,
      # 2022-06-15 fixed local_peak off-by-one so no need to -1
      original_row_num %in% (aux_ind_all_peaks) ~ TRUE,
      # the aux_ind_neg_peak is supposed to be FALSE then
      # original_row_num %in% (aux_ind_all_peaks) ~ FALSE,
      # all other values remain as the same
      TRUE ~ peak
    ))

  # Grab the index of the most neg peak
  # backtrack 1 to grab the actual peak
  aux_ind_neg_peak =
    min(aux_ind_all_peaks) # - 1

  # 2022-12-30 try to anchor the location for finding cutoff between the neg peak and the next closest?
  # Added arrange() above for aux_ind_all_peaks so we can grab the 2nd value for the second peak if it exists
  # 2023-03-23: Missing piece might be first checking if it's a flipped case -
  # For flipped case, only keep all indices to the left of the peak,
  # for regular case, keep all indices to the right of the peak
  # This method of subsetting should retain the order of the original vector
  if(dens_flip == TRUE){
    # Keep only the indices to the left of the peak (smaller)
    aux_ind_check_peaks =
      aux_ind_all_peaks[aux_ind_all_peaks <= aux_ind_neg_peak]
  }else{
    # if not flipped, keep only the indices to the right of the peak (larger)
    aux_ind_check_peaks =
      aux_ind_all_peaks[aux_ind_all_peaks >= aux_ind_neg_peak]
  }


  if(length(aux_ind_check_peaks) > 1){
    bound_cutoff = TRUE
    # aux_ind_second_peak = aux_ind_all_peaks[2]
    # If the most neg peak is = highest peak then take the next highest, otherwise take the highest peak
    if(aux_ind_check_peaks[1] == aux_ind_neg_peak){
      aux_ind_second_peak = aux_ind_check_peaks[2]
    }else{
      aux_ind_second_peak = aux_ind_check_peaks[1]
    }
  }else{
    bound_cutoff = FALSE
    if(dens_flip == TRUE){
      # Upper bound of the search when only 1 peak and dens flip is till the min row of the dataframe
      aux_ind_second_peak = -Inf

    }else{
      # Upper bound of the search when only 1 peak and no dens flip is till the max row of the dataframe
      aux_ind_second_peak = Inf
    }
  }

  # 2022-12-30: Rearrange the aux_ind_second_peak and aux_ind_neg_peak in case the ordering is flipped
  # due to the location of highest peak != neg peak
  aux_ind_bound =
    if(dens_flip == FALSE){
      range(c(aux_ind_neg_peak + 1, aux_ind_second_peak))
    }else{
      range(c(aux_ind_neg_peak - 1, aux_ind_second_peak))
    }

  # 2022-12-30: cutoff is max 2nd deriv/max curvature
  # 2022-04-14 the cutoff is at the start of the plateau after peak
  # plateau_pre == TRUE, min x_avg and then subtract 1 on the index
  aux_ind_cutoff_pre =
    new_d %>%
    # 2022-05-24 changed from >= to > b/c with the -1 (off by one error)
    # there is potential this would identify a cutoff before the most neg peak
    # 2022-09-07 when the density is flipped, we want the opposite, the left side which are
    #     row num less than the peak index on original row num
    # 2022-12-30 add conditions bounding the search on the filter()
    # Based on aux_ind_second_peak defined above
    # If 2nd peak exists, will bound between tallest and the next peak
    # If only 1 peak identified, then will bound search at min/max row depending on the dens_flip
    # purrr::when(
    #   dens_flip == FALSE  ~ dplyr::filter(., (original_row_num > (aux_ind_neg_peak + 1)) & (original_row_num < aux_ind_second_peak)),
    #   ~  dplyr::filter(., (original_row_num < (aux_ind_neg_peak - 1) & (original_row_num > aux_ind_second_peak)))
  #   # What if we put the indices of neg_peak and second_peak in order ?
  #   # Since it's already arranged, should not matter about whether dens_flip
  #   dens_flip == FALSE  ~ dplyr::filter(., (original_row_num > aux_ind_bound[1]) & (original_row_num <  aux_ind_bound[1])),
  # # 2023-03-14 Add a condition to check for right side of neg peak (left if flipped)
  dplyr::filter((original_row_num > aux_ind_bound[1]) & (original_row_num <= aux_ind_bound[2])) %>%
    #   ) %>%
    # # 2023-03-17 revert to 1st local max 2nd deriv
    # # 2022-12-30 What if we took max 2nd deriv as max curvature?
    # # Only between neg peak and the next peak to anchor the cutoff
    # dplyr::slice_max(second_deriv, n = 1, with_ties = FALSE) %>%
    # pull(original_row_num)
    # 2023-03-29: remove the condition of checking 2nd deriv >0 the point before point of interest
    # the problem of small blips in 2nd deriv might be solved with the optimized binning
    dplyr::slice_max(plateau_pre_2) %>%
    # 2022-09-07 when density flipped, want the x_avg that's max?
    purrr::when(
      dens_flip == FALSE ~ dplyr::slice_min(., x_avg),
      ~ dplyr::slice_max(., x_avg)
    ) %>%
    dplyr::pull(original_row_num)

  # 2023-03-30 Offset by 1 due to calculation of diff and derivs
  aux_ind_cutoff = aux_ind_cutoff_pre - 1

  # print(dim(new_d))
  # print(length(aux_ind_cutoff))

  # Creating a column for cutoff
  if(length(aux_ind_cutoff) == 0){
    new_d =
      new_d %>%
      dplyr::mutate(cutoff = FALSE)
  }else if (length(aux_ind_cutoff) > 0) {
    new_d =
      new_d %>%
      dplyr::mutate(
        cutoff = dplyr::case_when(
          original_row_num == aux_ind_cutoff ~ TRUE,
          TRUE ~ FALSE)

      )
  }

  return(new_d)
}
