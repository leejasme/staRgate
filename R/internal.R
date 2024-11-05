# Consolidate checks on argument inputs to functions
#' Internal function: Check the argument inputs for functions
#'
#' @keywords internal

check_inputs <- function(fun_name,
                         arg_list){

  # Create a list of args that need to be set to default values
  return_list <- list()

  # Based on the different function names, the checks are diff

  if(fun_name == "getDensityGates"){
    if (!(inherits(arg_list[["intens_dat"]], "data.frame"))) {
      rlang::abort(message="Error: `intens_dat` must be of data.frame class")
    }

    if (!all(arg_list[["marker"]] %in% colnames(arg_list[["intens_dat"]]))) {
      rlang::abort(message="Error: `marker` must be string matching column name(s) of `intens_dat`")
    }


    if (!(inherits(arg_list[["bin_n"]], "numeric"))) {
      rlang::warn(
        message=c("Warning: `bin_n` not specified as numeric.",
                  "i"="Default value `bin_n = 512` will be used.")
      )

      return_list <- append(return_list, list("bin_n"=512))
    }

    if (!(inherits(peak_detect_ratio, "numeric"))) {
      rlang::warn(
        message=c(
          "Warning: `peak_detect_ratio` not specified as numeric.",
          "i"="Default value `peak_detect_ratio = 10` will be used."
        )
      )
      # Set back to default
      return_list <- append(return_list, list("peak_detect_ratio"=10))
    }

    if (!(subset_col %in% colnames(intens_dat))) {
      rlang::abort(message="Error: `subset_col` must be string matching column name of `intens_dat`")
    }

    if (!is.null(neg_intensity_threshold) &
        !(inherits(neg_intensity_threshold, "numeric"))) {
      rlang::warn(
        message=c(
          "Warning: `neg_intensity_threshold` must be a numeric if supplied",
          "i"="Default value `neg_intensity_threshold = NULL` will be used (no filtering)."
        )
      )

      return_list <- append(return_list, list("neg_intensity_threshold"=NULL))
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


  }

  return(return_list)
}

# getDensityGates() calls a few internal functions but they are not exported
# To keep the package functions clean and less confusing about which functions to use
# Order of calls/operation:
# (1) getDensityGates() (exported function) calls:
# (2) getDensityMats (internal) calls (in the following order)
#     (3) getDensityDerivs() (internal)
#     (4) getDensityPeakCutoff() (internal)

#' Internal function: Estimate derivatives for density of `marker` for each unique subset of `subset_col`
#'
#' Internal function for `get_density_gates`
#' For each unique value in `subset_col`, estimate the derivatives for
#' `marker` (intensity values)
#'
#'
#' @param dens `density` object from the \link[stats]{density}
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
#'              default is 1800, which is on the biexponential scale
#'
#' @return list of dataframe with density estimation and corresponding 1st-4th derivatives,
#'         indicators of local peaks, plateau_pre \cr
#'         each element corresponds to each unique value of `subset_col` \cr
#'         for each dataframe: rows correspond to each of the bins
#'
#' @keywords internal

getDensityDerivs <- function(dens,
                             marker,
                             subset_col,
                             bin_n = 512,
                             peak_detect_ratio = 10,
                             pos_peak_threshold = 1800) {

  # Meant for internally calling within get_density_peaks and for debug/checking matrices/calculations
  # input data is diff- expecting the density obj for general subset_col
  x <- dens$x

  # x every 1
  x_binned <-
    range(x) |>
    (function(r){
      seq(r[[1]], r[[2]], by = 1)
    })()

  # Find interval
  new_d <-
    tibble::tibble(
      x = x,
      y = dens$y,
      x_int = findInterval(x, x_binned)
    ) |>
    dplyr::group_by(.data$x_int) |>
    dplyr::summarize(
      x_avg = mean(x),
      y_avg = mean(y)
    ) |>
    dplyr::arrange(x_avg) |>
    # 2022-09-06 add the original row num here
    # Want the original row num for merging back later
    # Based on the x continuous value might be problematic with rounding
    tibble::rownames_to_column(var = "original_row_num") |>
    dplyr::mutate(original_row_num = as.numeric(original_row_num)) |>
    dplyr::mutate(
      first_deriv = c(0, diff(y_avg) / diff(x_avg)),
      first_deriv_sign = sign(first_deriv),
      second_deriv = c(0, diff(first_deriv) / diff(x_avg)),
      second_deriv_sign = sign(second_deriv),
      third_deriv = c(0, diff(second_deriv) / diff(x_avg)),
      third_deriv_sign = sign(third_deriv),
      fourth_deriv = c(0, diff(third_deriv) / diff(x_avg)),
      # Peak if change from + to - for first deriv (-1 - 1 = -2)
      # second derive < 0 = peak
      local_peak = (c(0, diff(first_deriv_sign)) == -2) & second_deriv_sign < 0
    )

  # Simpliest fix to shift local_peak 1 index "up"
  new_d$local_peak <-
    c(new_d$local_peak[-1], FALSE)

  # # 2022-06-15 split out the plateau_pre step bc need to fix the off-by-1 for local_peak identified
  # new_d <-
  #   new_d |>
  #   dplyr::mutate(
  #     # Due to how the diff and deriv are calculated, the point of interest is back tracked by 1 so check 2nd deriv of 1 point before
  #     plateau_pre = ((c(0, diff(third_deriv_sign)) == -2) & (sign(fourth_deriv) < 0)),
  #     plateau_pre_2 = vapply(
  #       dplyr::row_number(),
  #       function(r) {
  #         if (r >= 2) {
  #           (plateau_pre[r] == TRUE & second_deriv_sign[(r - 1)] > 0)
  #         } else {
  #           FALSE
  #         }
  #       },
  #       logical(1)
  #     )
  #   )

  return(new_d)
}


#' Internal function: Matrix of calculations for density gating of intensity values in `marker` for each unique subset of `subset_col`
#'
#' Internal function for `getDensityGates`
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
#'         for each unique subset found in `subset_col` \cr
#'         rows correspond to unique values in `subset_col`, \cr
#'         cols correspond to the information for density gating
#' @keywords internal

getDensityMats <- function(intens_dat,
                           marker,
                           subset_col,
                           bin_n = 512,
                           peak_detect_ratio = 10,
                           pos_peak_threshold = 1800) {

  # Meant for using internally in getDensityGates but also for debug
  # to grab the full matrix of density, peak and cutoff
  # rows is each bin for the density estimation
  # cols for x (intens), y (dens) values, derivs (1-4th), sign of derivs,
  # indicators for local peak, plateau_pre (for calculating cutoff)
  # indicator of "real peaks" and cutoff
  # Same args as getDensityGates


  ## Calc density ---
  # Need to separate out per subset
  dens <-
    intens_dat |>
    dplyr::group_by(!!dplyr::sym(subset_col)) |>
    # each subset's data is separated out in `data`
    tidyr::nest() |>
    # calculate density for each subset separately
    dplyr::mutate(
      dens_obj = purrr::map(
        data,
        ~ stats::density(.x[[marker]],
          n = bin_n
        )
      )
    ) |>
    dplyr::ungroup()

  # Smooth the density by taking n bins = bin_n
  # calc 1st and 2nd deriv of each bin's avg y (density values) at avg x value
  # identify local peaks and plateuing
  dens_binned <-
    dens |>
    dplyr::mutate(
      dens_binned =
        purrr::map(
          dens_obj,
          function(d) {
            getDensityDerivs(d,
              subset_col = subset_col,
              marker = marker,
              bin_n = bin_n,
              peak_detect_ratio = peak_detect_ratio,
              pos_peak_threshold = pos_peak_threshold
            )
          }
        ),
      dens_peaks =
        purrr::map(
          dens_binned,
          function(d) {
            getDensityPeakCutoff(d,
              marker = marker,
              subset_col = subset_col,
              dens_flip = FALSE,
              bin_n = bin_n,
              peak_detect_ratio = peak_detect_ratio,
              pos_peak_threshold = pos_peak_threshold
            )
          }
        ),
      # 2022-11-01 check how many peaks
      # if > 2, reduce bin width before checking for desnity flip
      n_peaks =
        purrr::map(
          dens_peaks,
          function(d) {
            sum(d$peak)
          }
        ) |>
          unlist(),
      # bin_width = bin_size,
      bin_n = bin_n
    )

  # Check the peak value against threshold
  flag_peak <-
    dens_binned |>
    dplyr::select(!!dplyr::sym(subset_col), dens_peaks) |>
    dplyr::mutate(
      peak_loc =
        purrr::map(
          dens_peaks,
          function(x) {
            temp <- x |>
              dplyr::filter(peak == TRUE) |>
              dplyr::slice_min(x_avg) |>
              dplyr::pull(x_avg)

            # TO get around identifying 0 peaks
            if (length(temp) == 0) {
              return(NA_real_)
            } else {
              return(temp)
            }
          }
        ) |>
          unlist(),
      flag_pos_peak =
        peak_loc > pos_peak_threshold
    ) |>
    dplyr::select(-dens_peaks) |>
    dplyr::filter(flag_pos_peak == TRUE)

  # Recalculate with flipped binned intensity
  # Chose to not use flipped raw intensity bc
  # worried about binning on flipped intensity will result
  # in diff bins and not a 1-1 match on flipped bins when
  # replacing cutoff from the first attempt
  # To prevent unexpected errors with empty rows
  if (nrow(flag_peak) > 0) {
    dens_flipped <-
      dens_binned |>
      # Only filter to the subsets that needed to be recalculated
      dplyr::filter(
        # need a glue::glue to work
        !!rlang::parse_expr(glue::glue("{subset_col} %in% flag_peak[[subset_col]]"))
      ) |>
      dplyr::mutate(
        dens_peaks_flip =
          purrr::map(
            # 2023-03-30 dont use a recalculated density with --1*x_avg
            # Use original density but with dens_flip = TRUE arg to search on left vs. right
            dens_binned,
            function(d) {
              getDensityPeakCutoff(d,
                marker = marker,
                subset_col = subset_col,
                dens_flip = TRUE,
                bin_n = bin_n,
                peak_detect_ratio = peak_detect_ratio,
                pos_peak_threshold = pos_peak_threshold
              )
            }
          )
      ) |>
      dplyr::select(
        dplyr::all_of(subset_col),
        dens_peaks_flip
      )
  } else if (nrow(flag_peak) == 0) {
    dens_flipped <-
      dens |>
      dplyr::select(!!dplyr::sym(subset_col)) |>
      dplyr::mutate(
        dens_peaks_flip =
          rep(list(NULL), length(unique(dens[, subset_col])))
      )
  }

  # Replace cutoff col only for the subsets in dens_flipped
  dens_final <-
    dplyr::left_join(
      dens_binned,
      dens_flipped,
      by = subset_col
    ) |>
    dplyr::mutate(
      dens_peaks_final =
        if (all(vapply(dens_peaks_flip, is.null, FUN.VALUE = logical(1)))) {
          dens_peaks
        } else {
          purrr::pmap(
            list(
              dens_og = dens_peaks,
              dens_f = dens_peaks_flip
            ),
            function(dens_og, dens_f) {
              if (!is.null(dens_f)) {
                dplyr::left_join(
                  dens_og |>
                    dplyr::select(-cutoff),
                  dens_f |>
                    dplyr::select(
                      original_row_num,
                      cutoff
                    ),
                  by = "original_row_num"
                )
              } else {
                dens_og
              }
            }
          )
        }
    ) |>
    # Try to remove unnecessary parts
    dplyr::select(-dens_binned)

  # return(dens_binned)
  return(dens_final)
}

#' Internal function: Determine the "real peaks" and cutoff based on the density estimation and its derivs
#'
#' Internal function for `getDensityGates`
#'
#' @param dens_binned_dat list of dataframe output from the `getDensityDerivs`
#' @param marker string for the marker to gate on
#'              the name needs to match exactly the column name in `dens_binned_dat`
#' @param subset_col string for the column name to indicate the subsets to apply density gating on
#'              will perform operation on subsets corresponding to each unique value in column
#' @param bin_n numeric to be passed to `n` parameter of `density(n = bin_n)` for
#'              number of equally spaced points at which the density is to be estimated
#'              default is 512, which is the default of `density(n = 512)`
#' @param peak_detect_ratio numeric threshold for eliminating small peaks where
#'              a peak that is < than the highest peak by `peak_detect_ratio` times will be ignored
#'              default = 10
#' @param pos_peak_threshold numeric for threshold to identify a positive peak
#'              default is 1800, which is on the biexponential scale
#' @param dens_flip logical for whether the gating should be applied "backwards" where the peak is
#'              a positive peak and want to gate to the left of peak instead of right
#'
#' @return list of dataframe `dens_binned_dat` with additional columns added for
#'         peak(s) identified and the cutoff
#'         each element corresponds to each unique value of `subset_col`
#'         for each dataframe: rows correspond to each of the bins
#' @keywords internal

getDensityPeakCutoff <- function(dens_binned_dat,
                                 marker,
                                 subset_col,
                                 bin_n = 512,
                                 peak_detect_ratio = 10,
                                 pos_peak_threshold = 1800,
                                 dens_flip = FALSE) {

  # Identify "real peaks"
  new_d <-
    dens_binned_dat |>
    dplyr::filter(local_peak == TRUE) |>
    dplyr::arrange(-y_avg) |>
    # Get ratios relative to the highest peak
    dplyr::mutate(
      ratio = vapply(
        dplyr::row_number(),
        function(r) {
          y_avg[r] / y_avg[1]
        },
        numeric(1)
      )
    ) |>
    # Remove any peaks that are "too small"
    # threshold is at ratio < 1/peak_detect_ratio is "too small"
    # If the highest peak is > peak_detect_ratio times the height of the peak in question
    # then the peak is too small
    dplyr::filter(ratio >= 1 / peak_detect_ratio) |>
    # denoting these as "true peaks"
    dplyr::mutate(peak = TRUE) |>
    dplyr::select(original_row_num, peak) |>
    dplyr::left_join(
      dens_binned_dat,
      y = _,
      by = "original_row_num"
    ) |>
    # Replace the NA in the peak col to avoid confusion
    dplyr::mutate(
      peak = dplyr::case_when(is.na(peak) ~ FALSE,
                              .default=peak)
    )

  # 2022-04-14 for all peaks need to correct off-by-one error where "peak" is identified at the
  #   # 1 index over from actual peak
  aux_ind_all_peaks <-
    new_d |>
    dplyr::filter(peak == TRUE) |>
    # 2022-12-30 arrange by the height
    dplyr::arrange(-y_avg) |>
    dplyr::pull(original_row_num)

  new_d <-
    new_d |>
    dplyr::mutate(peak = dplyr::case_when(
      # set the aux_ind_neg_peak as peak = TRUE to correct for off-by-one error
      # 2022-06-15 fixed local_peak off-by-one so no need to -1
      original_row_num %in% (aux_ind_all_peaks) ~ TRUE,
      # the aux_ind_neg_peak is supposed to be FALSE then
      # all other values remain as the same
      TRUE ~ peak
    ))

  # Grab the index of the most neg peak
  aux_ind_neg_peak <-
    min(aux_ind_all_peaks)

  # Grab other peaks locations for anchoring the search for cutoff
  # For flipped case, only keep all indices to the left of the peak,
  # for regular case, keep all indices to the right of the peak
  if (dens_flip == TRUE) {
    # Keep only the indices to the left of the peak (smaller)
    aux_ind_check_peaks <-
      aux_ind_all_peaks[aux_ind_all_peaks <= aux_ind_neg_peak]
  } else {
    # if not flipped, keep only the indices to the right of the peak (larger)
    aux_ind_check_peaks <-
      aux_ind_all_peaks[aux_ind_all_peaks >= aux_ind_neg_peak]
  }

  # Anchor search points
  if (length(aux_ind_check_peaks) > 1) {
    bound_cutoff <- TRUE
    # If the most neg peak is = highest peak then take the next highest,
    # otherwise take the highest peak
    if (aux_ind_check_peaks[1] == aux_ind_neg_peak) {
      aux_ind_second_peak <- aux_ind_check_peaks[2]
    } else {
      aux_ind_second_peak <- aux_ind_check_peaks[1]
    }
  } else {
    bound_cutoff <- FALSE
    if (dens_flip == TRUE) {
      # Upper bound of the search when only 1 peak and dens flip is till
      # the min row of the dataframe
      aux_ind_second_peak <- -Inf
    } else {
      # Upper bound of the search when only 1 peak and no dens flip is till
      # the max row of the dataframe
      aux_ind_second_peak <- Inf
    }
  }

  # 2022-12-30: Rearrange the aux_ind_second_peak and aux_ind_neg_peak in
  # case the ordering is flipped
  # due to the location of highest peak != neg peak
  aux_ind_bound <-
    if (dens_flip == FALSE) {
      range(c(aux_ind_neg_peak + 1, aux_ind_second_peak))
    } else {
      range(c(aux_ind_neg_peak - 1, aux_ind_second_peak))
    }

  # IF `dens_flip` == FALSE
  # Place cutoff at max 2nd deriv where 1st deriv is still < 0
  # There are some instances where the 1st deriv > 0 so we need to filter
  # more specifically before finding max(2nd deriv)
  if(dens_flip == FALSE){
    aux_ind_cutoff <-
      new_d |>
      # In case there are multiple, would not want to simple subset to 1st deriv < 0
      # instead subset to the first set of rows where 1st deriv < 0
      dplyr::filter((original_row_num > aux_ind_bound[1]) & (original_row_num <= aux_ind_bound[2])) |>
      dplyr::arrange(original_row_num) |>
      dplyr::mutate(first_deriv_change = c(0, diff(first_deriv_sign)),
                    new_anchor = (first_deriv_change != 0) & (lag(first_deriv_sign) < 0)) |>
      (function(df){
        # If 1st deriv <0 for all of search region, no need to additional filter
        if(all(df$new_anchor == FALSE)){
          df
        }else{
          # Else, need to filter to only 1st deriv <0 before searching
          dplyr::filter(df, original_row_num < min(original_row_num[new_anchor == TRUE]))
        }
      })() |>
      # cutoff at max of 2nd deriv
      dplyr::slice_max(second_deriv) |>
      # There is no need to offset by 1 here because we are not calculating
      # diff in the derivs
      dplyr::pull(original_row_num)
  }else{
    # if `dens_flip` == TRUE, place cutoff at max 2nd deriv where 1st deriv >0
    aux_ind_cutoff <-
      new_d |>
      # In case there are multiple, would not want to simple subset to 1st deriv > 0
      # instead subset to the first set of rows where 1st deriv > 0
      dplyr::filter((original_row_num > aux_ind_bound[1]) & (original_row_num <= aux_ind_bound[2])) |>
      dplyr::arrange(-original_row_num) |>
      dplyr::mutate(first_deriv_change = c(0, diff(first_deriv_sign)),
                    new_anchor = (first_deriv_change != 0) & (stats::lag(first_deriv_sign) > 0)) |>
      (function(df){
        # If 1st deriv <0 for all of search region, no need to additional filter
        if(all(df$new_anchor == FALSE)){
          df
        }else{
          # Else, need to filter to only 1st deriv >0 before searching
          dplyr::filter(df, original_row_num > max(original_row_num[new_anchor == TRUE]))
        }
      })() |>
      # cutoff at max of 2nd deriv
      dplyr::slice_max(second_deriv) |>
      # There is no need to offset by 1 here because we are not calculating
      # diff in the derivs
      dplyr::pull(original_row_num)
  }


  # aux_ind_cutoff_pre <-
  #   new_d |>
  #   # # 2023-03-14 Add a condition to check for right side of neg peak (left if flipped)
  #   dplyr::filter((original_row_num > aux_ind_bound[1]) & (original_row_num <= aux_ind_bound[2])) |>
  #   dplyr::slice_max(plateau_pre_2) |>
  #   (function(df){
  #     if(dens_flip == FALSE){
  #       dplyr::slice_min(df, x_avg)
  #     }else{
  #       dplyr::slice_max(df, x_avg)
  #     }
  #   })()|>
  #   dplyr::pull(original_row_num)
  #
  # # 2023-03-30 Offset by 1 due to calculation of diff and derivs
  # aux_ind_cutoff <- aux_ind_cutoff_pre - 1

  # Creating a column for cutoff
  if (length(aux_ind_cutoff) == 0) {
    new_d <-
      new_d |>
      dplyr::mutate(cutoff = FALSE)
  } else if (length(aux_ind_cutoff) > 0) {
    new_d <-
      new_d |>
      dplyr::mutate(
        cutoff = dplyr::case_when(
          original_row_num == aux_ind_cutoff ~ TRUE,
          TRUE ~ FALSE
        )
      )
  }

  return(new_d)
}
