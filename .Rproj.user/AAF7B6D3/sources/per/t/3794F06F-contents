
get_density_peak_cutoff = function(dens_binned_dat,
                                   marker,
                                   subset_col = subset_col,
                                   bin_n = 512,
                                   peak_detect_ratio = 10,
                                   pos_peak_threshold = 1600,
                                   dens_flip = FALSE){

  #' Internal function: Determine the "real peaks" and cutoff
  #'
  #' Internal function for `get_density_gates`
  #'
  #' @param dens_binned_dat list of dataframe output from the `get_density_derivs`
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
  #'              default is 1600, which is on the biexponential scale
  #' @param dens_flip logical for whether the gating should be applied "backwards" where the peak is
  #'              a positive peak and want to gate to the left of peak instead of right
  #'
  #' @return list of dataframe `dens_binned_dat` with additional columns added for
  #'         peak(s) identified and the cutoff
  #'         each element corresponds to each unique value of `subset_col`
  #'         for each dataframe: rows correspond to each of the bins
  #' @export
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
    dplyr::mutate(peak = case_when(
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
