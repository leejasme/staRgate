run_pipeline = function(path_out,
                        path_gate_template,
                        path_fcs,
                        path_comp_mat,
                        path_biexp_params,
                        flowAI_yn = TRUE){
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

# optional flowAI step
  if(flowAI_yn){
    # Auto creates a flowAI folder within the path_out folder specified
    ## If dir doesn't exist, create one
    if(!dir.exists(here::here(glue("{path_out}/flowAI_results")))){
      dir.create(here::here(glue("{path_out}/flowAI_results")))
    }

    qc = flowAI::flow_auto_qc(flowWorkspace::gh_pop_get_data(gs),
                              folder_results = here::here(glue("{path_out}/flowAI_results")))

    ## Convert to a flowSet in order to convert back to GatingSet
    qc_fs =
      qc %>%
      # First convert to a flowSet
      flowCore::flowSet()

    ## Rename the sample
    flowCore::sampleNames(qc_fs) = pt_samp_nm

    ## convert the flowFrame obj returned in qc_fs to a GatingSet to pass to openCyto
    gs_qc =
      qc_fs %>%
      flowWorkspace::GatingSet()

    ## Remote the qc_fs object as it's no longer needed
    rm(qc_fs)

  }

  # Pre-Gate
  set.seed(glue::glue({format(Sys.Date(), format = "%Y%m%d")}))

  openCyto::gt_gating(gt_tcell, gs)

}
