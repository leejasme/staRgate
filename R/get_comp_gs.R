get_comp_gs = function(gs, path_comp_mat){

  # ## change "Viability" or "L/D" to = "L_D" for consistency with the comp mat label
  # {if(any(c("Viability", "L/D", "LD") %in% flowWorkspace::markernames(cs))){
  #   aux_log = flowWorkspace::markernames(cs) == "Viability" | flowWorkspace::markernames(cs) == "L/D" | flowWorkspace::markernames(cs) == "LD"
  #
  #   markernames(cs) = replace(markernames(cs), aux_log, "L_D")}
  # }
  #

  # If we take the GS as the "source of truth" for the marker to channel names mapping
  # we can clean up the comp mat according to the mapping instead of the "exceptions" listed above and below


  # Compensation
  ## Import comp. mat. csv exported from flowJo
  comp.mat = read.csv(path_comp_mat,
                      header = TRUE,
                      skip = 0) %>%
    ## Can remove the X col because that's the row names
    tibble::column_to_rownames(var = "X")

  ## Optional to check comp mat
  ## comp.mat %>%
  ##   head

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

  #
  # ## further clean up some of the col nanmes to match channel names in `gs`
  # comp.mat =
  #   comp.mat %>%
  #   dplyr::rename(stats::setNames(marker_chnl_names$colnms, gsub("\\.", "-", marker_chnl_names$chnl))) %>%
  #   {if("PE.Cy5.5.A" %in% marker_chnl_names$chnl) dplyr::rename(., "PE-Cy5.5-A" = "PE-Cy5-5-A")
  #     else dplyr::rename(.)} %>%
  #   {if("Alexa.Fluor.700.A" %in% marker_chnl_names$chnl) dplyr::rename(., "Alexa Fluor 700-A" = "Alexa-Fluor-700-A")
  #     else dplyr::rename(.)} %>%
  #   {if("Horizon.V450.A" %in% marker_chnl_names$chnl) dplyr::rename(., "Horizon V450-A" = "Horizon-V450-A")
  #     else dplyr::rename(.)} %>%
  #   {if("Pacific.Orange.A" %in% marker_chnl_names$chnl) dplyr::rename(., "Pacific Orange-A" = "Pacific-Orange-A")
  #     else dplyr::rename(.)} %>%
  #   {if("Qdot.605.A" %in% marker_chnl_names$chnl) dplyr::rename(., "Qdot 605-A" = "Qdot-605-A")
  #     else dplyr::rename(.)} %>%
  #   {if("Qdot.655.A" %in% marker_chnl_names$chnl) dplyr::rename(., "Qdot 655-A" = "Qdot-655-A")
  #     else dplyr::rename(.)}
  #

  # ## Check col names are equal
  # colnames(comp.mat) == (colnames(cs) %>% .[!(grepl("FSC|SSC|Time", .))])

  ## Create the `compensation` object with `flowCore::compensation()`
  comp = flowCore::compensation(comp.mat)

  ## Apply compensation to the `GatingSet` with `flowCore::compensate()`
  ### Must have the colnames in `comp.mat` match the channel names in the `cytoset`
  # IF they don't match, there will be an error stating which comp parameter (chnl) is not
  # in the gs
  gs = flowCore::compensate(gs, comp)

  return(gs)
}
