run_pipelin = function(){
  # String together the steps
  dtTemplate = data.table::fread(gtFile)

  ### Load
  gt_tcell = openCyto::gatingTemplate(gtFile)

  cs  <- flowWorkspace::load_cytoset_from_fcs(glue::glue("{path_fcs}"))

  ## Create a GatingSet of 1 sample
  gs = flowWorkspace::GatingSet(cs)


  # Apply comp -----
    # # Check no comp applied
    # flowWorkspace::gh_get_compensations(gs)
    gs = get_comp_gs(gs, path_comp_mat = path_comp_mat)

    # # Can check that the comp was applied
    # flowWorkspace::gh_get_compensations(gs)

  # Transformation ----
    # Check no transformation before
    gh_get_transformations(gs)

    # Apply biexp trans
    gs = get_biexpTransform_gs(gs, path_biexp_params = path_biexp_params)

    # Check transformation applied
    gh_get_transformations(gs)


}
