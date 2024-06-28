getCompGS = function(gs, path_comp_mat){
#' Applies Compensation using specifications in csv file provided at `path_comp_mat`
#'
#' The csv file at `path_comp_mat` should specify the channels to apply the compensation to.
#' The format is a matrix where the col and row names correspond to the channel names
#'
#' An example matrix is provided in the extdata/comp_mat_example_fcs.csv
#'
#' @param gs GatingSet to apply Biexponential Transformation to
#' @param path_comp_mat file path for .csv file that specifies the Compensation Matrix
#'
#' @return GatingSet with compensated data
#' @examples
#' \dontrun{
#' path_comp_mat = system.file("extdata", "comp_mat_example_fcs.csv",
#'                              package = "staRgate", mustWork = TRUE)
#'
#' # gs is a GatingSet object
#' gs = getCompGS(gs, path_comp_mat = path_comp_mat)
#' }
#'
#' @export

  # Compensation
  ## Import comp. mat. csv exported from flowJo
  comp.mat = utils::read.csv(path_comp_mat,
                             header = TRUE,
                             skip = 0) %>%
    ## Can remove the X col because that's the row names
    tibble::column_to_rownames(var = "X")

  # Apply comp mat only works if the chnl names in cs match the names in comp.mat
  # Can check two diff scenarios, if neither checks out then give error to ask user to correct

  # This should work generally for both cases because
  # Condition 1 is if the rownames are just the channel names w/out markers
  # There will be NAs for marker
  # And we can grab column chnl
  # Condition 2 is if we clean up the rownames(comp.mat) with
  # names <channel> :: <marker>
  marker_chnl_names =
    tibble::tibble(nms = rownames(comp.mat)) %>%
    tidyr::separate_wider_delim(cols = .data$nms,
                                names = c("chnl", "marker"),
                                # delim = "\\.\\.\\.\\.",
                                delim = " :: ",
                                too_few = "align_start",
                                too_many = "merge",
                                cols_remove = FALSE)

  # Replace the colnames of comp.mat with the chnl col
  # B/c it's a square matrix, they should always match up
  colnames(comp.mat) = marker_chnl_names$chnl

  ## Create the `compensation` object with `flowCore::compensation()`
  comp = flowCore::compensation(comp.mat)

  ## Apply compensation to the `GatingSet` with `flowCore::compensate()`
  ### Must have the colnames in `comp.mat` match the channel names in the `cytoset`
  # IF they don't match, there will be an error stating which comp parameter (chnl) is not
  # in the gs
  gs = flowCore::compensate(gs, comp)

  return(gs)
}
