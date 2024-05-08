getPerc = function(intens_dat,
                   num_marker,
                   denom_marker,
                   expand_num = FALSE,
                   expand_denom = FALSE,
                   keep_indicators = TRUE){
  #' Calculate the percentage of positive cells for specific subpopulations
  #'
  #' Expects data input same as the output from `get_gated_dat` with indicator columns of specific
  #' naming convention (see below).
  #'
  #' The subpopulations are defined as (num marker(s)) out of (denom marker(s)) where
  #' num denotes numerator, and denom denotes denominator
  #' (these shorthands are used in the function arguments)
  #'
  #' @param intens_dat dataframe of gated data with indicator columns per marker of interest
  #'                   (specify in `num_marker` and `denom_marker`)
  #'                   with naming convention `marker_pos` per marker with
  #'                   values of 0 to indicate negative-, 1 to indicate positive-expressing
  #' @param num_marker string for the marker(s) to specify the numerator for subpopulations of interest
  #'                   See `expand_num` argument and examples for how to specify
  #' @param denom_marker string for the marker(s) to specify the denominator for subpopulations of interest
  #'                     See `expand_denom` argument and examples for how to specify.
  #' @param expand_num logical, only accepts `TRUE` or `FALSE` with default of `FALSE`
  #'                   if `expand_num = TRUE`, currently only considers up to pairs of markers
  #'                   specified in `num_marker` in the numerator of subpopulation calculations
  #'                   (e.g., CD4+ & CD8- of CD3+)
  #'                   if `expand_num = FALSE`, only considers each marker specified in `num_marker`
  #'                   individually in the numerator of subpopulation calculations
  #'                   (e.g., CD4+ of CD3+)
  #' @param expand_denom logical, only accepts `TRUE` or `FALSE` with default of `FALSE`
  #'                     if `expand_denom = TRUE`, currently considers up to 1 marker from the `num_marker` and
  #'                     the unique combinations of `denom_marker` to generate list of subpopulations
  #'                     e.g., if `denom_marker = c("CD8")`, `num_marker = c("LAG3", "KI67")`, and `expand_denom = TRUE`,
  #'                     the subpopulations will include:
  #'                     1. LAG3+ of CD8+, LAG3- of CD8+, LAG3+ of CD8-, LAG3- of CD8-,
  #'                     2. KI67+ of CD8+, KI67- of CD8+, KI67+ of CD8-, KI67- of CD8-,
  #'                     3. KI67+ of (LAG3+ & CD8+), KI67- of (LAG3+ & CD8+), KI67+ of (LAG3+ & CD8-), KI67- of (LAG3+ & CD8-)...etc.,
  #'                     4. LAG3+ of (KI67+ & CD8+), LAG3- of (KI67+ & CD8+), LAG3+ of (KI67+ & CD8-), LAG3- of (KI67+ & CD8-)...etc.,
  #'                     if `expand_denom = FALSE`, only generates the list of subpopulations based on
  #'                     unique combinations of the `denom_marker`
  #'                     (e.g., `denom_marker = c("CD4")` and `expand_denom = FALSE` only considers subpopulations with
  #'                     denominator CD4+ and CD4- whereas
  #'                     `denom_marker = c("CD4", "CD8"` and `expand_denom = FALSE` will consider subpopulations with denominators
  #'                     (CD4- & CD8-), (CD4+ & CD8-), (CD4- & CD8+) and (CD4+ & CD8+))
  #' @param keep_indicators logical, only accepts `TRUE` or `FALSE` with default of `TRUE`
  #'                        if `keep_indicators = TRUE`, will return indicator columns of 0/1 to specify which markers are considered in the
  #'                        numerator and denominators of the subpopulations. Naming convention for the numerator cols are `<marker>_POS`
  #'                        and for denominator cols are `<marker>_POS_D`. For both sets of columns, `0` indicates considered the negative cells,
  #'                        `1` indicates considered the positive cells and `NA_real_` indicates not in consideration for the subpopulation.
  #'                        This is useful for matching to percentage data
  #'                        with potentially different naming conventions to avoid not having exact string matches for the same
  #'                        subpopulation
  #'                        take note that the order also matters when matching strings: "CD4+ & CD8- of CD3+" is different from "CD8- & CD4+ of CD3+"
  #' @return tibble containing the percentage of cells where
  #'         rows correspond to each subpopulation specified in the `subpopulation`,
  #'         `n_num` indicates the number of cells that satisifies the numerator conditions,
  #'         `n_denom` indicates the number of cells that satisifies the denominator conditions,
  #'         `perc` = `n_num` divided by `n_denom` unless `n_denom` = 0, then `perc = NA_real_`
  #'
  #' @export
  #'




  ## Check inputs ---
  # Check names in num_marker and denom_marker are in the data
  # more specifically, it's the <marker>_pos column that matters
  col_nms =
    toupper(colnames(intens_dat))  %>%
    .[grepl("_POS$", .)]

  # Change all names to full caps to avoid errors?
  intens_dat =
    intens_dat %>%
    janitor::clean_names(case = "all_caps")

  # Create the colnames that we expect for num_marker and denom_marker with _POS tagged on
  col_nms_subset =
    paste(toupper(c(num_marker, denom_marker)), "POS", sep = "_")

  if (!(inherits(intens_dat, "data.frame"))) {
    rlang::abort(message = "Error: `intens_dat` must be of data.frame class")
  }

  if (!(inherits(col_nms_subset, "character"))) {
    rlang::abort(
      c(message = "Error: `num_marker` and/or `denom_marker` must be of character class.",
        "i" = "`num_marker` and `denom_marker` must be strings containing the names of the marker(s) of interest for calculating the subpopulations.")
    )
  }

  if (!all((col_nms_subset %in% col_nms))) {
    rlang::abort(
      message = c(
        "Error: `num_marker` and `denom_marker` must have indicator columns in `intens_dat`",
        "i" = "`intens_dat` should contain 0/1 columns to indicate negative/positive with column name <marker>_pos for each marker specified in `num_marker` and `denom_marker`"
      )
    )
  }

  # check the cols corresponding to the num and denom markers only contain 0 and 1
  # what to do about NAs?
  if (!all(unique(unlist(intens_dat[, col_nms_subset])) %in% c(0, 1))) {
    rlang::abort(
      message = c(
        "Error: Unique values in indicator columns corresponding to `num_marker` and `denom_marker` should only contain 0 and 1.",
        "i" = glue::glue(
          "Currently detected unique values: {paste(unique(unlist(intens_dat[, col_nms_subset])), collapse = ', ')}"
        )
      )
    )
  }

  # Perhaps also need to check if the 0/1 cols are numeric?
  # For now, make it an error but can consider for the future to convert it with readr::parse_number and print a warning
  # But not sure how badly that could break down if the parsing did not show the expected values?
  if (!all(sapply(intens_dat[, col_nms_subset], class, simplify = TRUE) == "numeric")) {
    rlang::abort(message = "Error: Not all indicator columns corresponding to `num_marker` and `denom_marker` are numeric")
  }

  # check expand_num and expand_denom as well
  if (!inherits(expand_num, "logical")) {
    rlang::warn(
      message = c(
        "Warning: `expand_num` not of class logical (either `TRUE` or `FALSE`)",
        "i" = "Default value of `FALSE` will be used"
      )
    )

    expand_num = FALSE
  }

  if (!inherits(expand_denom, "logical")) {
    rlang::warn(
      message = c(
        "Warning: `expand_denom` not of class logical (either `TRUE` or `FALSE`)",
        "i" = "Default value of `FALSE` will be used"
      )
    )

    expand_denom = FALSE
  }

  # Can only expand numerator if the supplied num_marker has 2 or more markers
  # Since we are considering pairs
  # Can only expand_num and expand_denom if num_marker has 3+ markers bc need to consider pairs in the numerator and 1 marker rotated
  # Into the denom
  if(expand_num == TRUE & expand_denom == TRUE & length(num_marker) < 3){

    rlang::warn(
      message = c(
        "Warning: Both `expand_num` and `expand_denom` can be `TRUE` only if >= 3 markers are supplied in `num_marker",
        "i" = stringr::str_glue("`expand_num` = {expand_num}, `num_marker` length = {length(num_marker)}"),
        "i" = "Default values of `expand_num = FALSE` and `expand_denom = FALSE` will be used"
      )
    )

    expand_num = FALSE
    expand_denom = FALSE
  }


  if(expand_num == TRUE & length(num_marker) < 2){
    rlang::warn(
      message = c(
        "Warning: `expand_num` only applies if multiple markers are supplied in `num_marker",
        "i" = stringr::str_glue("`expand_num` = {expand_num}, `num_marker` length = {length(num_marker)}"),
        "i" = "Default value of `expand_num = FALSE` will be used"
      )
    )

    expand_num = FALSE
  }

  # Also if expand_denom == TRUE
  if(expand_denom == TRUE & length(num_marker) < 2){
    rlang::warn(
      message = c(
        "Warning: `expand_denom` only applies if multiple markers are supplied in `num_marker",
        "i" = stringr::str_glue("`expand_denom` = {expand_denom}, `num_marker` length = {length(num_marker)}"),
        "i" = "Default value of `expand_denom = FALSE` will be used"
      )
    )

    expand_denom = FALSE
  }

  # List the subpopulations first -----
  # Tag on _POS for the col names
  denom_cols = paste(denom_marker, "POS", sep = "_")
  num_cols = paste(num_marker, "POS", sep = "_")

  # Numerator
  # If expand_num == FALSE
  # Only consider each marker individually (regardless of status of the other markers)
  num_pre =
    data.frame(num_filters = c(paste0(num_cols, " == 0"),
                               paste0(num_cols, " == 1")))


  # Tag on the indicators
  num =
    purrr::map(num_cols,
               function(c) {
                 num_pre %>%
                   dplyr::mutate(!!c := ifelse(grepl(paste0(c, " == 1"), num_filters),
                                               1,
                                               ifelse(
                                                 grepl(paste0(c, " == 0"), num_filters),
                                                 0,
                                                 NA_real_
                                               )))
               }) %>%
    purrr::reduce(dplyr::left_join,
                  by = "num_filters") %>%
    # for the numerator, only expand if expand_num = TRUE
    # If expand_num = TRUE, then tag on the expanded.
    {
      if (expand_num) {
        dplyr::bind_rows(
          .,
          # purrr::map_dfc(num_cols,
          #                function(n) {
          #                  # We expect the col to be an indicator 0/1 col
          #                  # SetNames will ensure named col in map outcome
          #                  stats::setNames(data.frame(c(
          #                    paste(n, "0", sep = " == "),
          #                    paste(n, "1", sep = " == ")
          #                  )),
          #                  n)
          #                }) %>%
          #   # This creates the combos of all length(num_marker) so not
          #   # going to create pairs when we supply a num_marker > length 2
          #   # use utils::combn instead
          #   # tidyr::expand(., !!!.) %>%
          #   # {{}else{.}} %>%
          #   # form a col that puts together the conditions
          #   tidyr::unite(
          #     .,
          #     col = num_filters,
          #     sep = " & ",
          #     remove = FALSE
          #   ) %>%

          # Using utils::combn to grab pairs of num markers
          utils::combn(num_pre$num_filters, 2) %>%
            # apply across cols to paste into 1 condition
            apply(MARGIN = 2, function(x){paste(x, collapse = " & ")}) %>%
            data.frame(num_filters = .)
            ) %>%
          # also parse back the 0/1 from the original cols to add indicators
          dplyr::mutate(dplyr::across(dplyr::ends_with("_POS"),
                                      ~
                                        ifelse(
                                          grepl(paste0(dplyr::cur_column(), " == 0"),
                                                num_filters),
                                          0,
                                          ifelse(
                                            grepl(paste0(dplyr::cur_column(), " == 1"),
                                                  num_filters),
                                            1,
                                            NA_real_
                                          )
                                        )

          )
          ) %>%
          # filter out: there are the same markers in the two combos -
          # such as ctla4_pos == 0 & ctla4_pos == 1
          dplyr::rowwise() %>% # Need rowwise to work with c_across
          dplyr::filter(
            !(stringr::str_detect(num_filters, "&") & (sum(1*is.na(dplyr::c_across(dplyr::ends_with("_POS")))) == (length(num_marker)-1)))
          ) %>%
          dplyr::ungroup()

      } else .
    }


  # Use dplyr filters with strings
  denom =
    purrr::map_dfc(denom_cols,
                   function(d) {
                     # We expect the col to be an indicator 0/1 col
                     # SetNames will ensure named col in map outcome
                     stats::setNames(data.frame(c(
                       paste(d, "0", sep = " == "),
                       paste(d, "1", sep = " == ")
                     )),
                     d)
                   }) %>%
    # Can use base expand.grid instead of tidyr::expand
    # But will do need tidyr for other functions
    # Equivalent as tidyr::expand(denom, denom[[1]], denom[[2]],...)
    # Except this should take care of all n cols
    # So far have tested denom with 1-3 markers
    tidyr::expand(.,!!!.) %>%
    # form a col that puts together the conditions
    tidyr::unite(.,
                 col = denom_filters,
                 sep = " & ",
                 remove = FALSE) %>%
    # also parse back the 0/1 from the original cols to add indicators
    dplyr::mutate(dplyr::across(dplyr::ends_with("_POS"),
                                ~ ifelse(grepl("== 0", .x), 0, 1))) %>%
    # to identify the indicator is denominator, tag on the _D to create _POS_D cols
    dplyr::rename_with(~ paste0(., "_D"), dplyr::ends_with("_POS")) %>%
    # If expand_denom = TRUE then need to tag on 1 additional 0/1 for each numerator marker
    {
      if (expand_denom == TRUE) {
        denom_pre_expand = .

        # Tag on this expanded set to the denom
        # expand.grid will create all combinations of the two vectors
        dplyr::bind_rows(
          denom_pre_expand,
          expand.grid(denom_pre_expand$denom_filters, num_pre$num_filters) %>%
            tidyr::unite(
              .,
              col = denom_filters,
              sep = " & ",
              remove = FALSE
            ) %>%
            # Use the original denom_filters to pull the 0/1 cols
            dplyr::left_join(.,
                             denom_pre_expand,
                             by = c("Var1" = "denom_filters")) %>%
            # same thing with the num 0/1s
            dplyr::left_join(.,
                             num,
                             by = c("Var2" = "num_filters")) %>%
            dplyr::select(-Var1,-Var2) %>%
            # The numerator cols pulled from num need to be renamed with the _D
            dplyr::rename_with(~ paste0(., "_D"), dplyr::ends_with("_POS"))
        )

      } else .
    }

  # List out all the subpops by expanding grid of num and denom filters
  # Need to filter out any pairs that have marker in the numerator and denom
  # This should only happen if expand_denom == TRUE
  tbl_subpop =
    expand.grid(num$num_filters, denom$denom_filters) %>%
    # This is a rowwise operation
    dplyr::rowwise() %>%
    dplyr::filter(!grepl(
      # Just need to detect the var name part
      stringr::str_remove(Var1, pattern = " == [[:digit:]]"),
      Var2)
      ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      num_filters = Var1,
      denom_filters = Var2
    ) %>%
    # Grab the indicators again
    dplyr::left_join(.,
                     num,
                     by = "num_filters") %>%
    # Grab the indicators again
    dplyr::left_join(.,
                     denom,
                     by = "denom_filters") %>%
    # Drop any that have the same marker in num and denom
    # sapply() is used to create the conditions of checking if not NA in both num and denom
      # for each marker in num_cols,
      # this means there are repeated markers in num and denom
    # then use a paste() after the sapply() to collapse using collapse = "|" for OR conditions
    # then finaly paste with a ! in front before parsing expr with rlang
    {
      if(expand_denom == TRUE){
        dplyr::filter(.,
                      !!rlang::parse_expr(sapply(num_cols,
                                                 function(x){
                                                   glue::glue("(!is.na({x}) & !is.na({x}_D))")
                                                 }) %>%
                                            paste(., collapse = "|") %>%
                                            paste0("!(", ., ")"))
        )
      }else{
        .
      }
    } %>%
    # # Get marker names
    dplyr::mutate(
      num_label =
        stringr::str_replace_all(num_filters, pattern = c("_POS == 0" = "_NEG",
                                                          "_POS == 1" = "_POS",
                                                          " \\& " = "_")),
      denom_label =
        stringr::str_replace_all(denom_filters, pattern = c("_POS == 0" = "_NEG",
                                                            "_POS == 1" = "_POS",
                                                            " \\& " = "_")),
      subpopulation =
        paste0(num_label, "_OF_", denom_label)

    )

  # Count the n and N,
  # Is there a way to achieve n and N as a list with 1 map statement
  tbl_counts =
    purrr::map2_dfr(
    tbl_subpop$num_filters,
    tbl_subpop$denom_filters,
    function(n_filter, d_filter){
      current_d =
        intens_dat %>%
        # First get our denom
        dplyr::filter(!!rlang::parse_expr(d_filter))

        N = nrow(current_d)
        # n cells is how many out of the denom satisfy the numerator expression
        n = nrow(current_d %>% dplyr::filter(!!rlang::parse_expr(n_filter)))

        data.frame(
          # Keep these to merge back onto tbl_subpop
          num_filters = n_filter,
          denom_filters = d_filter,
          n_num = n,
          n_denom = N,
          perc =
            # In anticipation of N = 0, divide by zero
            if(N > 0){
              (n/N)*100
            }else{
              NA_real_
            }

        )
    }
  )

  # Merge on for final table
  tbl_final =
    dplyr::left_join(
      tbl_subpop,
      tbl_counts,
      by = c("denom_filters", "num_filters")
    ) %>%
    # If not keep indicators, remove the POS cols
    {
      if(!keep_indicators){
        dplyr::select(., -dplyr::ends_with("_POS"), -dplyr::ends_with("_POS_D"))
      }else{
        .
      }
    } %>%
    # Move the marker col to front, and remove the num/denom label and filter columns
    dplyr::select(subpopulation, n_num, n_denom, perc, dplyr::everything(), -num_label, -denom_label, -num_filters, -denom_filters)

  return(tbl_final)

}


# Code for testing
# four scenarios are T/F for expand_num and expand_denom
library(tidyverse)

intens_dat = tibble::tibble(
  CD3_pos = c(0, 1, 0, 1, 1, 1, 1),
  CD4_pos = c(0, 1, 0, 0, 1, 0, 0),
  CD8_pos = c(0, 0, 0, 0, 1, 0, 1),
  CD45RA_pos = c(0, 1, 1, 0, 0, 0, 1),
  icos_pos = rep(0, times = 7),
  cd25_pos = rep(1, times = 7),
  tim3_pos = rep(0, times = 7),
  cd27_pos = rep(0, times = 7),
  Cd57_pos = rep(0, times = 7),
  Cxcr5_pos = rep(0, times = 7),
  CCR4_pos = rep(0, times = 7),
  ccr7_pos = rep(0, times = 7),
  hladr_pos = rep(0, times = 7),
  cd28_pos = rep(0, times = 7),
  pd1_pos = rep(0, times = 7),
  LAG3_pos = c(0, 1, 0, 1, 0, 0, 0),
  cd127_pos = rep(0, times = 7),
  cd38_pos = rep(0, times = 7),
  tigit_pos = rep(0, times = 7),
  eomes_pos = rep(0, times = 7),
  ctla4_pos = c(1, 1, 1, 0, 0, 0, 1),
  foxp3_pos = rep(0, times = 7),
  gitr_pos =  c(1, 1, 0, 0, 0, 0, 0),
  tbet_pos = rep(0, times = 7),
  ki67_pos = rep(0, times = 7),
  gzm_b_pos = rep(0, times = 7)


)

# Some examples
# A simpler combos
denom_marker = c("CD4", "CD8") # "CD3",
num_marker = c("LAG3", "PD1", "CTLa4", "ki67")

# for 29-marker panel
denom_marker = c("CD4", "CD8") # "CD3",
num_marker = c("CD45RA", "ICOS", "CD25", "TIM3", "CD27", "CD57",
                    "CXCR5", "CCR4", "CCR7", "HLADR", "CD28", "PD1", "LAG3",
                    "CD127", "CD38", "TIGIT", "EOMES", "CTLA4", "FOXP3",
                    "GITR", "TBET", "KI67", "GZM_B")

# For 11-color
denom_marker = c("CD4", "CD8") # "CD3",
num_marker = c("ICOS", "TIM3", "PD1", "LAG3",
               "CTLA4", "FOXP3", "KI67")

# test for all 4 combos of expand_num and expand_denom
test =
  getPerc(intens_dat,
           num_marker = num_marker,
           denom_marker = denom_marker,
           expand_num = FALSE,
           expand_denom = FALSE,
           keep_indicators = FALSE)

# View(test)

dim(test) # 184

test_2 =
  getPerc(intens_dat,
           num_marker = num_marker,
           denom_marker = denom_marker,
           expand_num = TRUE,
           expand_denom = FALSE,
           keep_indicators = TRUE)

dim(test_2) # 4232

test_3 =
  getPerc(intens_dat,
           num_marker = num_marker,
           denom_marker = denom_marker,
           expand_num = FALSE,
           expand_denom = TRUE,
           keep_indicators = FALSE)

dim(test_3) # 8280

test_4 =
  getPerc(intens_dat,
           num_marker = num_marker,
           denom_marker = denom_marker,
           expand_num = TRUE,
           expand_denom = TRUE,
           keep_indicators = FALSE)

dim(test_4) # 182344

# For each scenario, the n combinations are where n = number of markers supplied in num_marker
# n_d = number of markers for denom in denom_markers
# expand_num = FALSE, expand_num = FALSE: (2*(2^(n_d)))*n
# expand_num = TRUE, expand_num = FALSE: [2*(2^(n_d ) )]*n + 2^(n_d )*(n choose 2)*(2^2)
# expand_num = FALSE, expand_num = TRUE: [2*(2^(n_d ) )]*n + 4*(2^(n_d))*(n)(n-1)
# expand_num = TRUE, expand_num = TRUE: [2*(2^(n_d ) )]*n + 2^(n_d )*(n choose 2)*(2^2) + 4*(2^(n_d))*(n)(n-1) + (2^(n_d + 1))((2^2)*n*((n-1) choose 2))

# For a 23-marker vector, and CD4/CD8 4-combo subset
# expand_num = FALSE, expand_num = FALSE: 2*2*2*n = 184
# expand_num = TRUE, expand_num = FALSE: 4232
# expand_num = FALSE, expand_num = TRUE: 8280
# expand_num = TRUE, expand_num = TRUE: 182344

# For the 11-color panel, we have 7 markers
# expand_num = FALSE, expand_num = FALSE:  56
# expand_num = TRUE, expand_num = FALSE: 392
# expand_num = FALSE, expand_num = TRUE: 728
# expand_num = TRUE, expand_num = TRUE: 4424


