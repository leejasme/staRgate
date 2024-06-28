getBiexpTransformGS = function(gs, path_biexp_params){
#' Applies Biexpeonential Transformation using specifications in csv file provided at `path_biexp_params`
#'
#' The csv file at `path_biexp_params` should specify the channels to apply the transformation to
#' and the parameters (negative decades, width basis and positive decades).
#'
#' The default is negative decades = 0.5, width basis = -30 and positive decades = 4.5.
#'
#' The Transformation can be applied to only a subset of the channels included in the
#' GatingSet.
#'
#' An example table is provided in the extdata/biexp_transf_parameters_x50.csv
#'
#' @param gs GatingSet to apply Biexponential Transformation to
#' @param path_biexp_params file path for .csv file that specifies the Biexponential Transformation
#'
#' @return GatingSet with Biexponentially Transformed data
#'
#' @examples
#' \dontrun{
#' path_biexp_params = system.file("extdata", "biexp_transf_parameters_x50.csv",
#'                                  package = "staRgate", mustWork = TRUE)
#'
#' # gs must be a GatingSet object
#' gs = getBiexpTransformGS(gs, path_biexp_params = path_biexp_params)
#' }
#'
#' @export

  # Transformation
  ## Read in the table of parameters
  tbl_biexp_params =
    # readxl::read_xlsx(path_biexp_params) %>%
    utils::read.csv(path_biexp_params) %>%
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

  ## Create transformerList in order to apply to gs
  trans = flowWorkspace::transformerList(names(biexpTrans), biexpTrans)

  ## Apply to gs
  gs = flowWorkspace::transform(gs, trans)

  return(gs)
}
