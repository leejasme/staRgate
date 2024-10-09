getDensityGates <- function(intens_dat,
                            marker,
                            subset_col,
                            bin_n=512,
                            peak_detect_ratio=10,
                            pos_peak_threshold=1800,
                            neg_intensity_threshold=-1000) {
  #' Density gating of intensity values in `marker` for each unique subset of
  #' `subset_col`
  #'
  #' For each unique value in `subset_col`, gate using density and estimated
  #' derivatives to identify cutoff at shoulder (i.e., point of tapering off)
  #' relative to the peak for `marker` (intensity values).
  #' The strategy of cutting at the shoulder mimics the strategy to gate
  #' relative to a unimodal background negative subpopulation, which is capable
  #' of capturing dim subpopulations.
  #'
  #' @param intens_dat dataframe of pre-gated (compensated, biexp. transf,
  #'                   openCyto steps) intensity values where
  #'                   cols=intensity value per marker,
  #'                   rows=each sample
  #' @param marker string for the marker(s) to gate on the names need to match
  #'               exactly the column name in `intens_dat`
  #' @param subset_col string for the column name to indicate the subsets to
  #'                   apply density gating on will perform operation on subsets
  #'                   corresponding to each unique value in column
  #' @param bin_n numeric to be passed to `n` parameter of `density(n=bin_n)`
  #'              for number of equally spaced points at which the density is to
  #'              be estimated \cr
  #'              Default is 512, which is the default of
  #'              `density(n=512)`
  #' @param peak_detect_ratio numeric threshold for eliminating small peaks where
  #'              a peak that is < than the highest peak by `peak_detect_ratio`
  #'              times will be ignored \cr
  #'              Default=10
  #' @param pos_peak_threshold either:
  #' \itemize{
  #'     \item numeric for threshold to identify a positive peak for all or
  #'     \item a dataframe if supplying multiple `marker` to gate. The dataframe
  #'           needs to be supplied with 2 columns named `marker` and
  #'           `pos_peak_threshold` and rows for the `marker` to gate
  #' }
  #' Default is 1800 (note this is on the biexponential scale) for all `marker`
  #' @param neg_intensity_threshold numeric for threshold to filter out any
  #'              "very negatively" expressed cells in the density estimation to
  #'              avoid over-compression and difficulty in distinguishing peaks
  #'              and the gates \cr
  #'              This is only applied as a filter for the density estimation,
  #'              the cells < `neg_intensity_threshold` are retained in the
  #'              intensity matrix for other steps \cr
  #'              Expects the `neg_intensity_threshold` is on the same scale as
  #'              the transformed data in `intens_dat` \cr
  #'              Default is `NULL`: no filters applied and density estimation
  #'              based on all cells in corresponding subsets.\cr
  #'              Suggested for biexp. transformed data is -1000 which
  #'              corresponds to ~-3300 on the original intensity scale)
  #'
  #' @return tibble of gates/cutoffs for `marker` for each unique subset found
  #'         in `subset_col` where
  #' \itemize{
  #'     \item rows correspond to unique values in `subset_col`
  #'     \item , columns correspond to`marker`
  #' }
  #' @export
  #'
  #' @examples
  #' # Create a fake dataset
  #' set.seed(100)
  #' intens_dat<-tibble::tibble(
  #'                CD3_pos=rep(c(0, 1), each=50),
  #'                CD4=rnorm(100, 100, 10),
  #'                CD8=rnorm(100, 100, 10)
  #' )
  #'
  #' # Run density gating, leaving other params at suggested defaults
  #' # number of bins suggested is 40 but default is at `bin_n=512`,
  #' # which is the default for the R base density() function
  #' getDensityGates(intens_dat, marker="CD4", subset_col="CD3_pos", bin_n=40)

  # getDensityGates() calls a few internal functions but they are not exported
  # These functions are still included in the internal.R file that users can see
  # In case there are any debugging needs, etc.
  # To keep the package functions clean and less confusing about which functions
  # to use
  # Order of calls/operation:
  # (1) getDensityGates() (exported function) calls:
  # (2) getDensityMats (internal) calls (in the following order)
  #     (3) getDensityDerivs() (internal)
  #     (4) getDensityPeakCutoff() (internal)

  ## Check inputs ---
  if (!(inherits(intens_dat, "data.frame"))) {
    rlang::abort(message="Error: `intens_dat` must be of data.frame class")
  }

  if (!all(marker %in% colnames(intens_dat))) {
    rlang::abort(message="Error: `marker` must be string matching column name(s) of `intens_dat`")
  }


  if (!(inherits(bin_n, "numeric"))) {
    rlang::warn(
      message=c("Warning: `bin_n` not specified as numeric.",
                "i"="Default value `bin_n=512` will be used.")
    )

    bin_n <- 512
  }

  if (!(inherits(peak_detect_ratio, "numeric"))) {
    rlang::warn(
      message=c(
        "Warning: `peak_detect_ratio` not specified as numeric.",
        "i"="Default value `peak_detect_ratio=10` will be used."
      )
    )
    # Set back to default-- have tested and is corrected
    peak_detect_ratio <- 10
  }

  if (!(subset_col %in% colnames(intens_dat))) {
    rlang::abort(message="Error: `subset_col` must be string matching column name of `intens_dat`")
  }

  if (!is.null(neg_intensity_threshold) &
    !(inherits(neg_intensity_threshold, "numeric"))) {
    rlang::warn(
      message=c(
        "Warning: `neg_intensity_threshold` must be a numeric if supplied",
        "i"="Default value `neg_intensity_threshold=NULL` will be used (no filtering)."
      )
    )

    neg_intensity_threshold <- NULL
  }

  if (!is.null(neg_intensity_threshold)) {
    i_dat <-
      intens_dat |>
      dplyr::filter(!(
        dplyr::if_any(dplyr::all_of(marker),
                      ~ .x < neg_intensity_threshold)
      ))
  } else {
    i_dat <- intens_dat
  }

  if (!(inherits(pos_peak_threshold, "numeric"))) {
    # If not a single numeric, then check if its a df of correct format, otherwise warn
    if (!(inherits(pos_peak_threshold, "data.frame"))) {
      rlang::warn(
        message=c(
          "Warning: `pos_peak_threshold` must be a numeric or dataframe",
          "i"="Default value `pos_peak_threshold=1800` will be used."
        )
      )

      pos_peak_threshold <-
        tibble::tibble(MARKER=marker, POS_PEAK_THRESHOLD=1800)
    } else {
      # If dataframe, then check if structured correctly
      # two columns with names marker and pos_peak_threshold
      # Change the colnames of pos_peak_thresholds to all caps if df
      pos_peak_threshold <-
        pos_peak_threshold |>
        janitor::clean_names(case="all_caps")

      chk1 <- all(colnames(pos_peak_threshold) %in% c("MARKER", "POS_PEAK_THRESHOLD"))

      # Also check if all the names in marker column match those in the `marker` arg
      chk2 <- all(marker %in% pos_peak_threshold$MARKER)

      # If the col names of pos_peak_threshold are correct but not all markers have a corresponding threshold
      if (!chk2 & chk1) {
        rlang::warn(
          message=c(
            "Error: Not all markers in `marker` argument have a `pos_peak_threshold` supplied",
            "i"="Default value of `pos_peak_threshold=1800` will be used for the markers without a value supplied."
          )
        )

        # which names are not in the `marker` arg?
        nms_to_fill <- marker[!(marker %in% pos_peak_threshold$MARKER)]

        # dplyr::add_row should add a row per string in the nms_to_fill vector with pos_peak_threshold fixed
        pos_peak_threshold <-
          pos_peak_threshold |>
          dplyr::add_row(MARKER=nms_to_fill, POS_PEAK_THRESHOLD=1800)
      } else if (!(chk1 & chk2)) {
        rlang::abort(
          message=c(
            "Error: `pos_peak_threshold` is not of correct format",
            "i"="Column names for `pos_peak_threshold` should be `marker` and `pos_peak_threshold`",
            "i"="Not all markers in `marker` argument have a `pos_peak_threshold` supplied."
          )
        )
      } else if (!chk1) {
        rlang::abort(
          message=c(
            "Error: `pos_peak_threshold` is not of correct format",
            "i"="Column names for `pos_peak_threshold` should be `marker` and `pos_peak_threshold`"
          )
        )
      }
    }
  } else {
    pos_peak_threshold <-
      tibble::tibble(MARKER=marker, POS_PEAK_THRESHOLD=pos_peak_threshold)
  }


  # Wrap around the `marker` supplied
  # make them lists
  dens_binned <-
    lapply(marker, function(m) {
      # Grab the pos peak threshold corresponding to the current marker to gate
      p_threshold <-
        pos_peak_threshold |>
        dplyr::filter(MARKER == m) |>
        dplyr::pull("POS_PEAK_THRESHOLD")

      # Still expects just 1 marker at a time and the pos peak threshold to be a numeric
      getDensityMats(
        i_dat,
        marker=m,
        subset_col,
        bin_n=bin_n,
        peak_detect_ratio=peak_detect_ratio,
        pos_peak_threshold=p_threshold
      )
    })

  # Name the list elements for easier grabbing?
  names(dens_binned) <- marker

  # Grab the cutoff or gates and return as a tibble
  cutoffs <-
    lapply(marker, function(m) {
      dens_binned[[m]] |>
        dplyr::mutate("{m}" :=
          purrr::map(.data$dens_peaks_final, function(d) {
            c <-
              d |>
              dplyr::filter(cutoff == TRUE) |>
              dplyr::select(x_avg)

            if (nrow(c) == 0) {
              return(NA_real_)
            } else {
              return(c)
            }
          }) |>
          unlist()) %>%
        dplyr::select(dplyr::all_of(c(subset_col, m)))
    }) %>%
    # Reduce to a dataframe
    purrr::reduce(., dplyr::left_join, by=subset_col)

  return(cutoffs)
}
