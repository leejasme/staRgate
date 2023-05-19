get_density_gates = function(intens_dat,
                             marker,
                             subset_col,
                             bin_n = 512,
                             peak_detect_ratio = 10,
                             pos_peak_threshold = 1600){
  #' Density gating of intensity values in `marker` for each unique subset of `subset_col`
  #'
  #' For each unique value in `subset_col`, gate using density and estimated derivatives
  #' to identify cutoff at shoulder (i.e., point of tapering off) after peak for
  #' `marker` (intensity values)
  #'
  #' @param intens_dat dataframe of pre-gated (compensated, biexp. transf, openCyto steps) intensity values where
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
  #'           '  default is 1600, which is on the biexponential scale
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
    rlang::abort(message = "Error: `marker` must be string matching column name of `intens_dat`")
  }

  # if(!(inherits(bin_size, "numeric"))){
  #   rlang::warn(message = c("Warning: `bin_size` not specified as numeric.",
  #                           "i" = "Default value `bin_size = 1` will be used."))
  #   # Set back to default-- have tested and is corrected
  #   bin_size = 1
  # }

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

  # 2022-04-16 wrapped above to get_density_peaks for debug usage
  dens_binned =
    get_density_mats(intens_dat,
                     marker,
                     subset_col,
                     # bin_size = bin_size,
                     bin_n = bin_n,
                     peak_detect_ratio = peak_detect_ratio,
                     pos_peak_threshold = pos_peak_threshold)

  # Grab the cutoff or gates and return as a tibble to match format in get_iso_ntil_gates()
  cutoffs =
    dens_binned %>%
    dplyr::mutate("{marker}" :=
                    purrr::map(dens_peaks_final,
                               function(d){
                                 c =
                                   d %>%
                                   dplyr::filter(cutoff == TRUE) %>%
                                   dplyr::select(x_avg)

                                 if(nrow(c) == 0){
                                   return(NA_real_)
                                 }else{
                                   return(c)
                                 }
                               }) %>%
                    unlist()
    ) %>%
    # dplyr::filter(cd4_pos_cd8_pos != "cd4_neg_cd8_neg") %>%
    dplyr::select(dplyr::all_of(c(subset_col, marker))) #%>%
  # rename(subpop = cd4_pos_cd8_pos)

  return(cutoffs)
}
