#' Calculate the percentage of positive cells for specific subpopulations
#'
#' Expects data input same as the output from `get_gated_dat` with indicator
#' columns of specific naming convention (see below).
#'
#' The subpopulations are defined as (num marker(s)) out of
#' (denom marker(s)) where num denotes numerator, and
#' denom denotes denominator
#' (these shorthands are used in the function arguments)
#'
#' @param intens_dat dataframe of gated data with indicator columns per
#'                   marker of interest
#'                   (specify in `num_marker` and `denom_marker`)
#'                   with naming convention `marker_pos` per marker with
#'                   values of 0 to indicate negative-, 1 to indicate
#'                   positive-expressing
#' @param num_marker string for the marker(s) to specify the numerator for
#'                   subpopulations of interest \cr
#'                   See `expand_num` argument and examples for how to specify
#' @param denom_marker string for the marker(s) to specify the denominator for
#'                   subpopulations of interest \cr
#'                   See `expand_denom` argument and examples for how to specify.
#' @param expand_num logical, only accepts `TRUE` or `FALSE` with default of `FALSE` \cr
#'                   if `expand_num=TRUE`, currently only considers up to pairs of markers
#'                   specified in `num_marker` in the numerator of subpopulation calculations
#'                   (e.g., CD4+ & CD8- of CD3+) \cr
#'                   if `expand_num=FALSE`, only considers each marker specified in `num_marker`
#'                   individually in the numerator of subpopulation calculations
#'                   (e.g., CD4+ of CD3+)
#' @param expand_denom logical, only accepts `TRUE` or `FALSE` with default of `FALSE` \cr
#'                     if `expand_denom=TRUE`, currently considers up to 1 marker from the `num_marker` and
#'                     the unique combinations of `denom_marker` to generate list of subpopulations \cr
#'                     e.g., if `denom_marker=c("CD8")`, `num_marker=c("LAG3", "KI67")`, and `expand_denom=TRUE`,
#'                     the subpopulations will include: \cr
#'                     1. LAG3+ of CD8+, LAG3- of CD8+, LAG3+ of CD8-, LAG3- of CD8-, \cr
#'                     2. KI67+ of CD8+, KI67- of CD8+, KI67+ of CD8-, KI67- of CD8-, \cr
#'                     3. KI67+ of (LAG3+ & CD8+), KI67- of (LAG3+ & CD8+), KI67+ of (LAG3+ & CD8-), KI67- of (LAG3+ & CD8-)...etc., \cr
#'                     4. LAG3+ of (KI67+ & CD8+), LAG3- of (KI67+ & CD8+), LAG3+ of (KI67+ & CD8-), LAG3- of (KI67+ & CD8-)...etc., \cr
#'                     if `expand_denom=FALSE`, only generates the list of subpopulations based on
#'                     unique combinations of the `denom_marker`
#'                     (e.g., `denom_marker=c("CD4")` and `expand_denom=FALSE` only considers subpopulations with
#'                     denominator CD4+ and CD4- whereas
#'                     `denom_marker=c("CD4", "CD8"` and `expand_denom=FALSE` will consider subpopulations with denominators
#'                     (CD4- & CD8-), (CD4+ & CD8-), (CD4- & CD8+) and (CD4+ & CD8+))
#' @param keep_indicators logical, only accepts `TRUE` or `FALSE` with default of `TRUE` \cr
#'                        if `keep_indicators=TRUE`, will return indicator columns of 0/1 to specify which markers are considered in the
#'                        numerator and denominators of the subpopulations.  \cr
#'                        Naming convention for the numerator cols are `<marker>_POS`
#'                        and for denominator cols are `<marker>_POS_D`. \cr
#'                        For both sets of columns, `0` indicates considered the negative cells,
#'                        `1` indicates considered the positive cells and `NA_real_`
#'                        indicates not in consideration for the subpopulation.  \cr
#'                        This is useful for matching to percentage data
#'                        with potentially different naming conventions to avoid not having exact string matches for the same
#'                        subpopulation \cr
#'                        Take note that the order also matters when matching strings: "CD4+ & CD8- of CD3+" is different from "CD8- & CD4+ of CD3+"
#' @return tibble containing the percentage of cells where
#' \itemize{
#'    \item rows correspond to each subpopulation specified in the `subpopulation`,
#'    \item `n_num` indicates the number of cells that satisifies the numerator conditions,
#'    \item `n_denom` indicates the number of cells that satisifies the denominator conditions,
#'    \item `perc`=`n_num` divided by `n_denom` unless `n_denom`=0, then `perc=NA_real_`
#' }
#'
#'
#' @export
#' @examples
#' library(dplyr)
#'
#' # Create a fake dataset
#' set.seed(100)
#' intens_dat <- tibble::tibble(
#'                CD3_pos=rep(c(0, 1), each=50),
#'                CD4=rnorm(100, 100, 10),
#'                CD8=rnorm(100, 100, 10)
#' )
#'
#' # Run getDensityGates to obtain the gates
#' gates <- getDensityGates(intens_dat, marker="CD4", subset_col="CD3_pos", bin_n=40)
#'
#' # Tag on the 0/1 on intens_dat
#' intens_dat_2 <- getGatedDat(intens_dat, cutoffs=gates, subset_col="CD3_pos")
#'
#' # Get percentage for CD4 based on gating
#' getPerc(intens_dat_2, num_marker=c("CD4"), denom_marker="CD3")

getPerc <- function(intens_dat,
                    num_marker,
                    denom_marker,
                    expand_num=FALSE,
                    expand_denom=FALSE,
                    keep_indicators=TRUE) {

  # Check names in num_marker and denom_marker are in the data
  # more specifically, it's the <marker>_pos column that matters
  col_nms <-
    toupper(colnames(intens_dat))[grepl("_POS$", (toupper(colnames(intens_dat))))]

  # Change all names to full caps to avoid errors?
  intens_dat <-
    intens_dat |>
    janitor::clean_names(case="all_caps")

  ## Check inputs ---
  checkInputs(intens_dat=intens_dat,
              num_marker=num_marker,
              denom_marker=denom_marker,
              expand_num=expand_num,
              expand_denom=expand_denom,
              keep_indicators=keep_indicators)

  # Create the colnames that we expect for num_marker and denom_marker with _POS tagged on
  col_nms_subset <-
    paste(toupper(c(num_marker, denom_marker)), "POS", sep="_")

  # use default values if incorrect
  if (!is.logical(expand_num)) {
    expand_num <- FALSE
  }

  if (!is.logical(expand_denom)) {
    expand_denom <- FALSE
  }

  if (!is.logical(keep_indicators)) {
    keepindicators <- TRUE
  }

  # Can only expand numerator if the supplied num_marker has 2 or more markers
  # Since we are considering pairs
  # Can only expand_num and expand_denom if num_marker has 3+ markers bc need to consider pairs in the numerator and 1 marker rotated
  # Into the denom
  if (expand_num == TRUE &
    expand_denom == TRUE & length(num_marker) < 3) {
    rlang::warn(
      message=c(
        "Both `expand_num` and `expand_denom` can be `TRUE` only if >= 3 markers are supplied in `num_marker",
        "i"=stringr::str_glue(
          "`expand_num`={expand_num}, `num_marker` length={length(num_marker)}"
        ),
        "i"="Default values of `expand_num=FALSE` and `expand_denom=FALSE` will be used"
      )
    )

    expand_num <- FALSE
    expand_denom <- FALSE
  }

  if (expand_num == TRUE & length(num_marker) < 2) {
    rlang::warn(
      message=c(
        "`expand_num` only applies if multiple markers are supplied in `num_marker",
        "i"=stringr::str_glue(
          "`expand_num`={expand_num}, `num_marker` length={length(num_marker)}"
        ),
        "i"="Default value of `expand_num=FALSE` will be used"
      )
    )

    expand_num <- FALSE
  }

  # Also if expand_denom == TRUE
  if (expand_denom == TRUE & length(num_marker) < 2) {
    rlang::warn(
      message=c(
        "`expand_denom` only applies if multiple markers are supplied in `num_marker",
        "i"=stringr::str_glue(
          "`expand_denom`={expand_denom}, `num_marker` length={length(num_marker)}"
        ),
        "i"="Default value of `expand_denom=FALSE` will be used"
      )
    )

    expand_denom <- FALSE
  }

  # List the subpopulations first -----
  # Tag on _POS for the col names
  denom_cols <- paste(denom_marker, "POS", sep="_")
  num_cols <- paste(num_marker, "POS", sep="_")

  # Numerator
  # If expand_num == FALSE
  # Only consider each marker individually (regardless of status of the other markers)
  num_pre <-
    data.frame(num_filters=c(paste0(num_cols, " == 0"), paste0(num_cols, " == 1")))

  # Tag on the indicators
  num <-
    purrr::map(num_cols, function(c) {
      num_pre |>
        dplyr::mutate(!!c :=
                        dplyr::case_when(grepl(paste0(c, " == 1"),
                                               .data$num_filters) ~ 1,
                                         grepl(paste0(c, " == 0"),
                                               .data$num_filters) ~ 0,
                                         .default=NA_real_)
        )
    }) |>
    purrr::reduce(dplyr::left_join, by="num_filters") |>
    # for the numerator, only expand if expand_num=TRUE
    # If expand_num=TRUE, then tag on the expanded.
    (function(df){
      if (expand_num) {
        dplyr::bind_rows(
          df,
          # Using utils::combn to grab pairs of num markers
          utils::combn(num_pre$num_filters, 2) |>
            # apply across cols to paste into 1 condition
            apply(MARGIN=2, function(x) {
              paste(x, collapse=" & ")
            }) |>
            data.frame(num_filters=_)
        ) |>
          # also parse back the 0/1 from the original cols to add indicators
          dplyr::mutate(dplyr::across(dplyr::ends_with("_POS"), ~
            # ifelse(
            #   grepl(paste0(dplyr::cur_column(), " == 0"), .data$num_filters),
            #   0,
            #   ifelse(grepl(
            #     paste0(dplyr::cur_column(), " == 1"), .data$num_filters
            #   ), 1, NA_real_)
            # )
            dplyr::case_when(
              grepl(paste0(dplyr::cur_column(), " == 0"), .data$num_filters) ~ 0,
              grepl(paste0(dplyr::cur_column(), " == 1"), .data$num_filters) ~ 1,
              .default=NA_real_
            )
            )
            ) |>
          # filter out: there are the same markers in the two combos -
          # such as ctla4_pos == 0 & ctla4_pos == 1
          dplyr::rowwise() |> # Need rowwise to work with c_across
          dplyr::filter(!(stringr::str_detect(.data$num_filters, "&") &
            (sum(
              1 * is.na(dplyr::c_across(dplyr::ends_with("_POS")))
            ) == (
              length(num_marker) - 1
            )))) |>
          dplyr::ungroup()
      } else {
        df
      }
    })()

  # Use dplyr filters with strings
  denom <-
    purrr::map_dfc(denom_cols, function(d) {
      # We expect the col to be an indicator 0/1 col
      # SetNames will ensure named col in map outcome
      stats::setNames(data.frame(c(
        paste(d, "0", sep=" == "), paste(d, "1", sep=" == ")
      )), d)
    }) |>
    # Can use base expand.grid instead of tidyr::expand
    # But will do need tidyr for other functions
    # Equivalent as tidyr::expand(denom, denom[[1]], denom[[2]],...)
    # Except this should take care of all n cols
    # So far have tested denom with 1-3 markers
    (function(x){tidyr::expand(x, !!!x)})() |>
    # form a col that puts together the conditions
    tidyr::unite(
      col=denom_filters,
      sep=" & ",
      remove=FALSE
    ) |>
    # also parse back the 0/1 from the original cols to add indicators
    dplyr::mutate(dplyr::across(dplyr::ends_with("_POS"),
                                ~ dplyr::case_when(grepl("== 0", .x) ~ 0,
                                                   .default=1)
                                )) |>
    # to identify the indicator is denominator, tag on the _D to create _POS_D cols
    # dplyr::rename_with(~ paste0("_D"), dplyr::ends_with("_POS")) |>
    dplyr::rename_with(function(nms){paste0(nms, "_D")}, dplyr::ends_with("_POS")) |>
    # If expand_denom=TRUE then need to tag on 1 additional 0/1 for each numerator marker
    (function(df){
      if (expand_denom == TRUE) {
        denom_pre_expand <- df

        # Tag on this expanded set to the denom
        # expand.grid will create all combinations of the two vectors
        dplyr::bind_rows(
          denom_pre_expand,
          expand.grid(denom_pre_expand$denom_filters, num_pre$num_filters) |>
            tidyr::unite(
              col=denom_filters,
              sep=" & ",
              remove=FALSE
            ) |>
            # Use the original denom_filters to pull the 0/1 cols
            dplyr::left_join(denom_pre_expand, by=c("Var1"="denom_filters")) |>
            # same thing with the num 0/1s
            dplyr::left_join(num, by=c("Var2"="num_filters")) |>
            dplyr::select(-"Var1", -"Var2") |>
            # The numerator cols pulled from num need to be renamed with the _D
            dplyr::rename_with(function(nms){paste0(nms, "_D")}, dplyr::ends_with("_POS"))
        )
      } else {
        df
      }
    })()

  # List out all the subpops by expanding grid of num and denom filters
  # Need to filter out any pairs that have marker in the numerator and denom
  # This should only happen if expand_denom == TRUE
  tbl_subpop <-
    expand.grid(num$num_filters, denom$denom_filters) |>
    # This is a rowwise operation
    dplyr::rowwise() |>
    dplyr::filter(!grepl(
      # Just need to detect the var name part
      stringr::str_remove(.data$Var1, pattern=" == [[:digit:]]"),
      .data$Var2
    )) |>
    dplyr::ungroup() |>
    dplyr::rename(
      num_filters="Var1",
      denom_filters="Var2"
    ) |>
    # Grab the indicators again
    dplyr::left_join(num, by="num_filters") |>
    # Grab the indicators again
    dplyr::left_join(denom, by="denom_filters") |>
    # Drop any that have the same marker in num and denom
    # vapply() is used to create the conditions of checking if not NA in both num and denom
    # for each marker in num_cols,
    # this means there are repeated markers in num and denom
    # then use a paste() after the vapply() to collapse using collapse="|" for OR conditions
    # then finally paste with a ! in front before parsing expr with rlang
    (function(df){
      if (expand_denom == TRUE) {
        dplyr::filter(df, !!rlang::parse_expr(
          vapply(num_cols, function(x) {
            glue::glue("(!is.na({x}) & !is.na({x}_D))")
          },
          FUN.VALUE=character(1)
          ) |>
            paste(collapse="|") |>
          (function(p){paste0("!(", p, ")")})()
        ))
      } else {
        df
      }
    })() |>
    # Get marker names
    dplyr::mutate(
      num_label =
        stringr::str_replace_all(
          .data$num_filters,
          pattern=c(
            "_POS == 0"="_NEG",
            "_POS == 1"="_POS",
            " \\& "="_"
          )
        ),
      denom_label =
        stringr::str_replace_all(
          .data$denom_filters,
          pattern=c(
            "_POS == 0"="_NEG",
            "_POS == 1"="_POS",
            " \\& "="_"
          )
        ),
      subpopulation =
        paste0(.data$num_label, "_OF_", .data$denom_label)
    )

  # Count the n and N,
  # Is there a way to achieve n and N as a list with 1 map statement
  tbl_counts <-
    purrr::map2_dfr(tbl_subpop$num_filters, tbl_subpop$denom_filters, function(n_filter, d_filter) {
      current_d <-
        intens_dat |>
        # First get our denom
        dplyr::filter(!!rlang::parse_expr(d_filter))

      N <- nrow(current_d)
      # n cells is how many out of the denom satisfy the numerator expression
      n <- nrow(current_d |> dplyr::filter(!!rlang::parse_expr(n_filter)))

      data.frame(
        # Keep these to merge back onto tbl_subpop
        num_filters=n_filter,
        denom_filters=d_filter,
        n_num=n,
        n_denom=N,
        perc =
        # In anticipation of N=0, divide by zero
          if (N > 0) {
            (n / N) * 100
          } else {
            NA_real_
          }
      )
    })

  # Merge on for final table
  tbl_final <-
    dplyr::left_join(tbl_subpop, tbl_counts, by=c("denom_filters", "num_filters")) |>
    # If not keep indicators, remove the POS cols
    (function(df){
      if (!keep_indicators) {
        dplyr::select(
          df,
          -dplyr::ends_with("_POS"),
          -dplyr::ends_with("_POS_D")
        )
      } else {
        df
      }
    })() |>
    # Move the marker col to front, and remove the num/denom label and filter columns
    dplyr::select(
      "subpopulation",
      "n_num",
      "n_denom",
      "perc",
      dplyr::everything(),
      -"num_label",
      -"denom_label",
      -"num_filters",
      -"denom_filters"
    )

  return(tbl_final)
}
