get_marker_indicators = function(dat, col_nm = "marker", markers_nm){
  #' Get indicators of marker positive/negative in numerator and denominator for cell subpopulations
  #'
  #' Useful for matching subpopulation names ignoring order and naming of markers
  #'
  #' Inputs:
  #' @param dat dataframe to add indicators to
  #' expects rows to correspond to a subpopulation
  #' @param col_nm string indicating column in `dat` (i.e., `dat$col_nm`) to find the subpopulation names to apply indicators on
  #' Default is 'marker'
  #' Expectation that the naming convention is <num>_OF_<denom>
  #' @param markers_nm string vector for names of markers to generate indicators for
  #' No default currently
  #'
  #' @return dataframe `dat` returned with additional columns for indicators of positivity
  #' Each value of `markers_nm` will have 2 columns to indicate positivity in numerator or denominator
  #' Naming conventions are `tolower(marker_name)_pos` for positivity indicator in numerator
  #' and `tolower(marker_name)_pos_d` for positivity in denominator
  #' @export

  # Check inputs
  if(!(inherits(dat, "data.frame") & inherits(col_nm, "character"))){
    rlang::abort(message = "Error: `dat` must be of data.frame class, `col_nm` must be a character")
  }

  # split the string to numerator and denominator
  aux_dat =
    dat %>%
    tidyr::separate(col = !!dplyr::sym(col_nm),
                    into = c("num", "denom"),
                    sep = "_OF_",
                    remove = FALSE)

  # # iterate through every marker ?
  aux_dat =
    purrr::map(
      markers_nm,
      ~ dplyr::mutate(aux_dat,
                      # the created indicator cols will be all lower to follow convention in the iso_gated and
                      # intensity_dat
                      # but comparisons will be done on all upper to avoid any problems with caps
                      "{tolower(.x)}_pos" := dplyr::case_when(grepl(glue::glue("{toupper(.x)}_POS"), toupper(num)) ~ 1,
                                                              grepl(glue::glue("{toupper(.x)}_NEG"), toupper(num)) ~ 0,
                                                              TRUE ~ NA_real_),
                      "{tolower(.x)}_pos_d" := dplyr::case_when(grepl(glue::glue("{toupper(.x)}_POS"), toupper(denom)) ~ 1,
                                                                grepl(glue::glue("{toupper(.x)}_NEG"), toupper(denom)) ~ 0,
                                                                TRUE ~ NA_real_)) %>%
        # Only select the unique identifier and the indicator cols
        dplyr::select(!!dplyr::sym(col_nm), dplyr::ends_with("_pos"), dplyr::ends_with("_pos_d"))
    )


  # Join all indicator cols together iteratively
  # avoiding plyr::join_all()
  aux_dat_reduced = aux_dat[[1]]

  for(i in 2:length(aux_dat)){
    aux_dat_reduced =
      dplyr::left_join(
        aux_dat_reduced,
        aux_dat[[i]],
        by = col_nm
      )
  }

  return(aux_dat_reduced)
}
