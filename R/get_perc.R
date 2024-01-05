get_perc = function(intens_dat,
                    num_marker,
                    denom_marker,
                    expand_num = FALSE,
                    expand_denom = FALSE) {
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
  #'                     (e.g., if `denom_marker = c("CD8")`, `num_marker = c("LAG3", "KI67")`, and `expand_denom = TRUE`,
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
  #' @param keep_indicators logical, only accepts `TRUE` or `FALSE` with default of `FALSE`
  #'                        if `keep_indicators = TRUE`, will return indicator columns of 0/1 to specify which markers are considered in the
  #'                        numerator and denominators of the subpopulations. This is useful for matching to percentage data
  #'                        with potentially different naming conventions to avoid not having exact string matches for the same
  #'                        subpopulation
  #'                        take note that the order also matters when matching strings: "CD4+ & CD8- of CD3+" is different from "CD8- & CD4+ of CD3+"
  #' @return tibble containing the percentage of cells where
  #'         rows correspond to each subpopulation
  #' @export

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

  if(!(inherits(intens_dat, "data.frame"))){
    rlang::abort(message = "Error: `intens_dat` must be of data.frame class")
  }

  if(!(inherits(col_nms_subset, "character"))){
    rlang::abort(c(message = "Error: `num_marker` and/or `denom_marker` must be of character class.",
                   "i" = "`num_marker` and `denom_marker` must be strings containing the names of the marker(s) of interest for calculating the subpopulations.")
    )
  }

  if (!(col_nms_subset %in% col_nms)) {
    rlang::abort(
      message = c(
        "Error: `num_marker` and `denom_marker` must have indicator columns in `intens_dat`",
        "i" = "`intens_dat` should contain 0/1 columns to indicate negative/positive with column name <marker>_pos for each marker specified in `num_marker` and `denom_marker`"
      )
    )
  }

  # check the cols corresponding to the num and denom markers only contain 0 and 1
  # what to do about NAs?
  if(!all(unique(unlist(intens_dat[, col_nms_subset])) %in% c(0, 1))){
    rlang::abort(
      message = c(
        "Error: Unique values in indicator columns corresponding to `num_marker` and `denom_marker` should only contain 0 and 1.",
        "i" = glue::glue("Currently detected unique values: {paste(unique(unlist(intens_dat[, col_nms_subset])), collapse = ', ')}")
      )
    )
  }

  # Perhaps also need to check if the 0/1 cols are numeric?
  # For now, make it an error but can consider for the future to convert it with readr::parse_number and print a warning
  # But not sure how badly that could break down if the parsing did not show the expected values?
  if(!all(sapply(intens_dat[, col_nms_subset], class, simplify = TRUE) == "numeric")){
    rlang::abort(message = "Error: Not all indicator columns corresponding to `num_marker` and `denom_marker` are numeric")
  }

  # check expand_num and expand_denom as well
  if(!inherits(expand_num, "logical")){
    rlang::warn(message = c("Warning: `expand_num` not of class logical (either `TRUE` or `FALSE`)",
                            "i" = "Default value of `FALSE` will be used"
                            )
    )

    expand_num = FALSE
  }

  if(!inherits(expand_denom, "logical")){
    rlang::warn(message = c("Warning: `expand_denom` not of class logical (either `TRUE` or `FALSE`)",
                            "i" = "Default value of `FALSE` will be used")
    )

    expand_denom = FALSE
  }

  # List the subpopulations first
  # Use dplyr filters with strings


}
