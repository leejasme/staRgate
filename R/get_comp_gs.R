get_comp_gs = function(gs, path_comp_mat){
  #' @export

  # Compensation
  ## Import comp. mat. csv exported from flowJo
  comp.mat = read.csv(path_comp_mat,
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
    # tibble::tibble(colnms = colnames(comp.mat)) %>%
    # tidyr::separate_wider_delim(cols = colnms,
    #                            names = c("chnl", "marker"),
    #                            # delim = "\\.\\.\\.\\.",
    #                            delim = "....",
    #                            too_few = "align_start",
    #                            too_many = "merge",
    #                            cols_remove = FALSE)
    # switch to using rownames b/c the names should be read in exactly rather than cols
    # That can changed when reading in such as the space converted to a .
    tibble::tibble(nms = rownames(comp.mat)) %>%
    tidyr::separate_wider_delim(cols = nms,
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
