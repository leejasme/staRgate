
get_density_mats = function(intens_dat,
                            marker,
                            subset_col,
                            bin_n = 512,
                            peak_detect_ratio = 10,
                            pos_peak_threshold = 1800){
  #' Internal function: Matrix of calculations for density gating of intensity values in `marker` for each unique subset of `subset_col`
  #'
  #' Internal function for `get_density_gates`
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
  #' @export
  # Meant for using internally in get_density_gates but also for debug
  # to grab the full matrix of density, peak and cutoff
  # rows is each bin for the density estimation
  # cols for x (intens), y (dens) values, derivs (1-4th), sign of derivs,
  # indicators for local peak, plateau_pre (for calculating cutoff)
  # indicator of "real peaks" and cutoff
  # Same args as get_density_gates


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
                     get_density_derivs(d,
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
            get_density_peak_cutoff(d,
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
  #     get_density_derivs(dens_binned$dens_obj[[i]], dens_flip = FALSE,
  #                        marker = marker,
  #                        bin_size = bin_size_multipeaks,
  #                        peak_detect_ratio = peak_detect_ratio,
  #                        pos_peak_threshold = pos_peak_threshold
  #     )
  #
  #
  #   dens_binned$dens_peaks[[i]] =
  #     get_density_peak_cutoff(dens_binned$dens_binned[[i]], marker = marker, dens_flip = FALSE,
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
        #        get_density_derivs(d,
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
              get_density_peak_cutoff(d,
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
