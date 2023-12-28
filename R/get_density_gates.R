get_density_gates = function(intens_dat,
                             marker,
                             subset_col,
                             bin_n = 512,
                             peak_detect_ratio = 10,
                             pos_peak_threshold = 1800,
                             neg_intensity_threshold = -1000){
  #' Density gating of intensity values in `marker` for each unique subset of `subset_col`
  #'
  #' For each unique value in `subset_col`, gate using density and estimated derivatives
  #' to identify cutoff at shoulder (i.e., point of tapering off) after peak for
  #' `marker` (intensity values)
  #'
  #' @param intens_dat dataframe of pre-gated (compensated, biexp. transf, openCyto steps) intensity values where
  #'              cols = intensity value per marker,
  #'              rows = each sample
  #' @param marker string for the marker(s) to gate on
  #'              the names need to match exactly the column name in `intens_dat`
  #' @param subset_col string for the column name to indicate the subsets to apply density gating on
  #'              will perform operation on subsets corresponding to each unique value in column
  #' @param bin_n numeric to be passed to `n` parameter of `density(n = bin_n)` for
  #'              number of equally spaced points at which the density is to be estimated
  #'              default is 512, which is the default of `density(n = 512)`
  #' @param peak_detect_ratio numeric threshold for eliminating small peaks where
  #'              a peak that is < than the highest peak by `peak_detect_ratio` times will be ignored
  #'              default = 10
  #' @param pos_peak_threshold either numeric for threshold to identify a positive peak for all or
  #'              a dataframe if supplying multiple `marker` to gate
  #'              Dataframe needs to be supplied with 2 columns named `marker` and `pos_peak_threshold`
  #'              and rows for the corresponding `marker` to gate
  #'              default is 1800 (note this is on the biexponential scale)
  #' @param neg_intensity_threshold numeric for threshold to filter out any "very negatively" expressed
  #'              cells in the density estimation to avoid over-compression and difficulty in distinguishing
  #'              peaks and the gates
  #'              This is only applied as a filter for the density estimation, the cells < `neg_intensity_threshold``
  #'              are retained in the intensity matrix for other steps
  #'              Expects the `neg_intensity_threshold` is on the same scale as the transformed data in `intens_dat`
  #'              Default is `NULL`: no filters applied and density estimation based on all cells in
  #'              corresponding subsets.
  #'              Suggested for biexp. transformed data is -1000 which corresponds to ~-3300 on the original intensity scale)
  #'
  #' @return tibble of gates/cutoffs for `marker` for each unique subset found in `subset_col`
  #'         rows correspond to unique values in `subset_col`, cols correspond to `marker`
  #' @export
  # Purpose: gate density of intensity value using 1st and 2nd derivatives
  #          to find the "shoulder point" after the negative peak as cutoff/gate
  #          this method does not require an isotype
  #
  # Inputs:
  # @intens_dat: df of pre-gated (compensated, biexp transf, gated CD4/CD8) intensity values where
  #              cols = intensity value per marker,
  #              rows = each sample
  #              assumes there are cols cd4_pos and cd8_pos of 0/1 for neg/pos CD4 and CD8
  # @marker: string for the marker to gate on
  #          the name needs to match exactly the column name in `intens_dat`
  # @bin_size: numeric for bin sizes when smoothing the density before taking derivs
  #            default = 120 (on the biexp transf scale)
  # @peak_detect_ratio: threshold for eliminating small peaks where
  #            a peak that is < than the highest peak by `peak_detect_ratio` times will be ignored
  #            default = 100
  # @pos_peak_threshold: numeric for threshold on identifying a positive peak
  #           current default is 1600 which still need to be checked/discussed as reasonable with MA
  #
  # Outputs:
  # tibble of gates for CD4+/CD8-, CD8+/CD4- and double pos
  #
  # Process: 1. Calc density of `marker` intensity values
  #          2. Smooth the density into bins of size `bin_size`,
  #             calc avg intensity and avg density in each bin
  #          3. calc first and second derivs, identify local peaks & shoulder pts after each peak
  #          4. identify highest peak & eliminate "small peaks" using `peak_detect_ratio`
  #          5. identify most negative peak and the cutoff as the 1st shoulder pt after the neg peak
  #          6. returns cutoffs to pass into get_gated_dat()

  ## Check inputs ---
  if(!(inherits(intens_dat, "data.frame"))){
    rlang::abort(message = "Error: `intens_dat` must be of data.frame class")
  }

  if(!(marker %in% colnames(intens_dat))){
    rlang::abort(message = "Error: `marker` must be string matching column name(s) of `intens_dat`")
  }


  if(!(inherits(bin_n, "numeric"))){
    rlang::warn(message = c("Warning: `bin_n` not specified as numeric.",
                            "i" = "Default value `bin_n = 512` will be used."))

    bin_n = 512
  }

  if(!(inherits(peak_detect_ratio, "numeric"))){
    rlang::warn(message = c("Warning: `peak_detect_ratio` not specified as numeric.",
                            "i" = "Default value `peak_detect_ratio = 10` will be used."))
    # Set back to default-- have tested and is corrected
    peak_detect_ratio = 10
  }

  # 2023-05-19 Add a check for `subset_col`
  if(!(subset_col %in% colnames(intens_dat))){
    rlang::abort(message = "Error: `subset_col` must be string matching column name of `intens_dat`")
  }

  # 2023-12-28 Add a check for `neg_intensity_threshold`
  if(!is.null(neg_intensity_threshold) & !(inherits(neg_intensity_threshold, "numeric"))){
    rlang::warn(message = c( "Warning: `neg_intensity_threshold` must be a numeric if supplied",
                            "i" = "Default value `neg_intensity_threshold = NULL` will be used (no filtering)."))

    neg_intensity_threshold = NULL
  }

  # 2023-12-28 add filtering to the function for very neg values
  # Only affects the density estimation and subsequent peak and gates identifications
  # This step is useful if there are a few outlier cells far off from the density
  # density estimation will lead to a compressed density squished in within a few bins
  # and the rest will be flat thus making it hard to detect peak and any features of the density
  if(!is.null(neg_intensity_threshold)){
    i_dat =
      intens_dat %>%
      dplyr::filter(!(if_any(all_of(markers_to_gate), ~.x < neg_intensity_thres)))

  }else{
    i_dat = intens_dat
  }

  # 2023-12-28 check pos_peak_threshold is structured correctly if it is not a numeric
  if(!(inherits(pos_peak_threshold, "numeric"))){
    # If not a single numeric, then check if its a df of correct format, otherwise warn
    if(!(iinherits(pos_peak_threshold, "data.frame"))){
      rlang::warn(message = c( "Warning: `pos_peak_threshold` must be a numeric or dataframe",
                               "i" = "Default value `pos_peak_threshold = 1800` will be used."))

      pos_peak_threshold = 1800
    }else{
      # If dataframe, then check if structured correctly
      # two columns with names marker and pos_peak_threshold
      # Change the colnames of pos_peak_thresholds to all caps if df
      pos_peak_thresholds =
        pos_peak_thresholds %>%
        janitor::clean_names(case = "all_caps")

      chk1 = all(colnames(pos_peak_thresholds) %in% c("MARKER", "POS_PEAK_THRESHOLD"))

      # Also check if all the names in marker column match those in the `marker` arg
      chk2 = all( marker %in% pos_peak_thresholds$marker)

      if(!(chk1 & chk2)){
        rlang::abort(message = c("Error: `pos_peak_thresholds` is not of correct format",
                                 "i" = "Column names for `pos_peak_thresholds` should be `marker` and `pos_peak_threshold`",
                                 "i" = "Not all markers in `marker` argument have a `pos_peak_threshold` supplied."))
      }else if(!chk1){
        rlang::abort(message = c("Error: `pos_peak_thresholds` is not of correct format",
                                 "i" = "Column names for `pos_peak_thresholds` should be `marker` and `pos_peak_threshold`"))
      }else if(!chk2){
        rlang::warn(message = c("Error: `pos_peak_thresholds` is not of correct format",
                                "i" = "Not all markers in `marker` argument have a `pos_peak_threshold` supplied.",
                                "i" = "Default value of `pos_peak_threshold = 1800` will be used for the markers without a value supplied."))

        # which names are not in the `marker` arg?
        nms_to_fill = marker[!(marker %in% pos_peak_thresholds$marker)]

        # dplyr::add_row should add a row per string in the nms_to_fill vector with pos_peak_threshold fixed
        pos_peak_thresholds =
          pos_peak_thresholds %>%
          dplyr::add_row(
            MARKER = nms_to_fill,
            POS_PEAK_THRESHOLD = 1800
          )
      }
    }
  }


  # TODO: Wrap around the `marker` supplied
  # make them lists
  dens_binned =
    lapply(marker,
           function(m){
             # Grab the pos peak threshold corresponding to the current marker to gate
             p_threshold =
               pos_peak_threshold %>%
               dplyr::filter(MARKER == m) %>%
               dplyr::pull(POS_PEAK_THRESHOLD)

            # Still expects just 1 marker at a time and the pos peak threshold to be a numeric
             get_density_mats(i_dat,
                              m,
                              subset_col,
                              bin_n = bin_n,
                              peak_detect_ratio = peak_detect_ratio,
                              pos_peak_threshold = p_threshold)
           }
           )

  # Name the list elements for easier grabbing?
  names(dens_binned) = marker

  # Grab the cutoff or gates and return as a tibble to match format in get_iso_ntil_gates()
  cutoffs =
    lapply(marker,
           function(m) {
             dens_binned[[m]] %>%
               dplyr::mutate("{m}" :=
                               purrr::map(dens_peaks_final,
                                          function(d) {
                                            c =
                                              d %>%
                                              dplyr::filter(cutoff == TRUE) %>%
                                              dplyr::select(x_avg)

                                            if (nrow(c) == 0) {
                                              return(NA_real_)
                                            } else{
                                              return(c)
                                            }
                                          }) %>%
                               unlist()) %>%
               # dplyr::filter(cd4_pos_cd8_pos != "cd4_neg_cd8_neg") %>%
               dplyr::select(dplyr::all_of(c(subset_col, marker))) #%>%
             # rename(subpop = cd4_pos_cd8_pos)
           }
    ) %>%
    # Reduce to a dataframe
    purrr::reduce(.,
                  dplyr::left_join,
                  # TO TEST- will it work if the subset_col is a string arg?
                  by = subset_col)

  return(cutoffs)
}
