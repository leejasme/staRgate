get_gated_dat = function(intens_dat = intensity_dat,
                         cutoffs,
                         subset_col){

  #' Attach indicator columns to `intens_dat` based on gates provided in `cutoffs`
  #'
  #' Adds an indicator column (0/1) to `intens_dat` for each marker in `cutoffs`
  #' as indicated by the columns in `cutoffs`
  #'
  #' The naming convention for the tagged on indicator columns will be
  #' `tolower(<marker_name>_pos)` where
  #' 0 indicates negativity or intensity < gate provided
  #' 1 indicates positivity or intensity > gate provided
  #'
  #' @param intens_dat dataframe of pre-gated (compensated, biexp. transf, openCyto steps) intensity values where
  #'        rows = each sample
  #' @param cutoffs tibble of gates/cutoffs for all markers to gate
  #'        Expects `cutoffs` to match format of output from [get_density_gates()] with
  #'        column corresponding to a marker, and rows to the subsets defined in the
  #'        `subset_col`
  #' @param subset_col string for the column name to indicate the subsets to apply density gating on
  #'        will perform operation on subsets corresponding to each unique value in column
  #' @return `intens_dat` with additional cols attached for each marker in `cutoffs`
  #' @export

  # Generalize getting 0/1 indicator for gated data based on cutoff
  # To replace get_iso_gated_mat


  ## Grab the markers in cutoffs
  mrks = (colnames(cutoffs) %>% .[!(.%in% c(subset_col, "subpop"))])

  ## Check inputs ---
  if(!(inherits(intens_dat, "data.frame"))){
    rlang::abort(message = "Error: `intens_dat` must be of data.frame class")
  }

  if(!(inherits(cutoffs, "data.frame"))){
    rlang::abort(message = "Error: `cutoffs` must be of data.frame class")
  }

  if(!all((mrks %in% colnames(intens_dat)))){
    rlang::abort(message = "Error: `cutoffs` must have marker names matching names in `intens_dat`")
  }

  # 2023-05-19 Add a check for `subset_col`
  if(!(subset_col %in% colnames(intens_dat))){
    rlang::abort(message = "Error: `subset_col` must be string matching column name of `intens_dat`")
  }


  # 2024-01-05 Add a check for `subset_col` on the cutoffs
  if(!(subset_col %in% colnames(cutoffs))){
    rlang::abort(message = "Error: `subset_col` must be string matching column name of `cutoffs`")
  }

  # 2023-05-19 should we build a capability for the subset_col to be passed in for
  # lazy eval (i.e., not string but just unquoted)


  ## use nest() to subset into separate dfs
  intens_dat %>%
    dplyr::group_by(.data[[subset_col]]) %>%
    tidyr::nest() %>%
    # Join the cutoffs
    dplyr::mutate(
      gated_data =
        purrr::map2(
          # 2023-05-19 does !!sym() in map work?
          !!rlang::sym(subset_col),
          data,
          function(s, d){
            # Filter to cutoffs for this subpop first
            c =
              cutoffs %>%
              # 2023-05-19 does !!sym() in filter work?
              dplyr::filter(!!rlang::sym(subset_col) == s)

            # If nrow(c) == 0 then add NAs for mrk_pos
            d %>%
              dplyr::mutate(
                dplyr::across(dplyr::all_of(mrks),
                       .fns = ~ {if(nrow(c) > 0){(.x >= c[[dplyr::cur_column()]])*1}else{NA_real_}},
                       # By specifying .names, we get a new col instead of overwrite
                       .names = "{tolower(.col)}_pos")
              )
          }
        )
    ) %>%
    # dont need the org data bc the indicators are added to the data as add'l cols
    dplyr::select(-data) %>%
    # use unnest() to combine back the data
    tidyr::unnest(cols = c(gated_data)) %>%
    # Ungroup data
    dplyr::ungroup()
}

