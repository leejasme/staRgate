get_biexpTransform_gs = function(gs, path_biexp_params){

  # Transformation
  ## Read in the table of parameters
  tbl_biexp_params =
    # readxl::read_xlsx(path_biexp_params) %>%
    read.csv(path_biexp_params) %>%
    janitor::clean_names(case = "all_caps")

  ## loop through all the channels and create a biexpTrans object per channel
  biexpTrans =
    purrr::pmap(
      list(p = tbl_biexp_params$POSITIVE_DEC,
           n = tbl_biexp_params$EXT_NEG_DEC,
           wb = tbl_biexp_params$WIDTH_BASIS),
      function(p, n, wb, ...){
        flowWorkspace::flowjo_biexp_trans(pos = p,
                                          neg = n,
                                          widthBasis = wb)

      }
    )


  ## Give names
  names(biexpTrans) = tbl_biexp_params$FULL_NAME


  ## Optional - Checks that all channels have a matching biexpTrans object
  # lapply(biexpTrans,
  #        FUN = function(x){x[["transform"]] %>% attr(., "parameters")}) %>%
  #   dplyr::bind_rows(.id = "channel")
  #
  # ## Do channel names match the names of biexpTrans?
  # all(names(biexpTrans) %in% parameters(comp))
  # all(names(biexpTrans) %in% names(markernames(gs)))

  ## Create transformerList in order to apply to gs
  trans = flowWorkspace::transformerList(names(biexpTrans), biexpTrans)

  ## Optional - Check the values before transformation
  # gs %>%
  #   gh_pop_get_data() %>%
  #   summary()

  ## Apply to gs
  gs = flowWorkspace::transform(gs, trans)

  return(gs)
}
