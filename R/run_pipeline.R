run_pipeline = function(path_gate_template,
                        path_fcs,
                        path_comp_mat,
                        path_biexp_params){
  # String together the steps
  dtTemplate = data.table::fread(path_gate_template)

  ### Load
  gt_tcell = openCyto::gatingTemplate(path_gate_template)

  cs  <- flowWorkspace::load_cytoset_from_fcs(path_fcs)

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
