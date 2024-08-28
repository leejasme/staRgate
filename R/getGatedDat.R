getGatedDat <- function(intens_dat=intensity_dat, cutoffs, subset_col) {
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
  #'        rows=each cell and cols are the intensity values for each marker
  #' @param cutoffs tibble of gates/cutoffs for all markers to gate \cr
  #'        Expects `cutoffs` to match format of output from [getDensityGates()] with
  #'        column corresponding to a marker, and rows to the subsets defined in the
  #'        `subset_col`
  #' @param subset_col string for the column name to indicate the subsets to apply density gating on
  #'        will perform operation on subsets corresponding to each unique value in column
  #' @return `intens_dat` with additional columns attached for each marker in `cutoffs`
  #' @export
  #' @examples
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
  #' # intens_dat_2 now has the cd4_pos tagged on
  #' head(intens_dat_2)


  ## Grab the markers in cutoffs
  mrks <- (colnames(cutoffs) %>% .[!(. %in% c(subset_col, "subpop"))])

  ## Check inputs ---
  if (!(inherits(intens_dat, "data.frame"))) {
    rlang::abort(message="Error: `intens_dat` must be of data.frame class")
  }

  if (!(inherits(cutoffs, "data.frame"))) {
    rlang::abort(message="Error: `cutoffs` must be of data.frame class")
  }

  if (!all((mrks %in% colnames(intens_dat)))) {
    rlang::abort(message="Error: `cutoffs` must have marker names matching names in `intens_dat`")
  }

  if (!(subset_col %in% colnames(intens_dat))) {
    rlang::abort(message="Error: `subset_col` must be string matching column name of `intens_dat`")
  }

  if (!(subset_col %in% colnames(cutoffs))) {
    rlang::abort(message="Error: `subset_col` must be string matching column name of `cutoffs`")
  }

  ## use nest() to subset into separate dfs
  intens_dat %>%
    dplyr::group_by(.data[[subset_col]]) %>%
    tidyr::nest() %>%
    # Join the cutoffs
    dplyr::mutate(
      gated_data =
        purrr::map2(!!rlang::sym(subset_col), data, function(s, d) {
          # Filter to cutoffs for this subpop first
          c <-
            cutoffs %>%
            dplyr::filter(!!rlang::sym(subset_col) == s)

          # If nrow(c) == 0 then add NAs for mrk_pos
          d %>%
            dplyr::mutate(dplyr::across(dplyr::all_of(mrks),
              .fns=~ {
                if (nrow(c) > 0) {
                  (.x >= c[[dplyr::cur_column()]]) * 1
                } else {
                  NA_real_
                }
              }, # By specifying .names, we get a new col instead of overwrite
              .names="{tolower(.col)}_pos"
            ))
        })
    ) %>%
    # dont need the org data bc the indicators are added to the data as add'l cols
    dplyr::select(-.data$data) %>%
    # use unnest() to combine back the data
    tidyr::unnest(cols=c(gated_data)) %>%
    # Ungroup data
    dplyr::ungroup()
}
