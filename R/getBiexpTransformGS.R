#' Applies Biexpeonential Transformation using specifications in csv file
#' provided at `path_biexp_params`
#'
#' The csv file at `path_biexp_params` should specify the channels to apply
#' the transformation to and the parameters
#' (negative decades, width basis and positive decades).
#' The default is negative decades=0.5, width basis=-30 and
#' positive decades=4.5.
#' The Transformation can be applied to only a subset of the channels
#' included in the GatingSet.
#'
#' An example table is provided in the extdata/biexp_transf_parameters_x50.csv
#'
#' @param gs GatingSet to apply Biexponential Transformation to
#' @param path_biexp_params file path for .csv file that specifies the
#'        Biexponential Transformation
#'
#' @return GatingSet with Biexponentially Transformed data
#'
#' @examples
#' # This example does not contain all the pre-processing steps required in
#' # getting the GatingSet (gs) ready for Biexp transformation.
#' # To see the steps that are required to creating the (gs),
#' # please see the vignette for a full tutorial
#'
#' # To make this a runnable example, read in the FCS file to create gs and
#' # directly apply
#'
#' # File path to the FCS file
#' path_fcs <- system.file("extdata",
#'                         "example_fcs.fcs",
#'                         package="staRgate",
#'                         mustWork=TRUE)
#' path_biexp_params <- system.file("extdata",
#'                                  "biexp_transf_parameters_x50.csv",
#'                                  package="staRgate",
#'                                  mustWork=TRUE)
#'
#' # Create a cytoset then convert to gs
#' cs <- flowWorkspace::load_cytoset_from_fcs(path_fcs)
#' gs <- flowWorkspace::GatingSet(cs)
#'
#' # gs must be a GatingSet object
#' gs <- getBiexpTransformGS(gs, path_biexp_params=path_biexp_params)
#'
#' # To check the transformation parameters applied
#' flowWorkspace::gh_get_transformations(gs)
#' @export

getBiexpTransformGS <- function(gs, path_biexp_params) {
  ## Read in the table of parameters
  tbl_biexp_params <-
    utils::read.csv(path_biexp_params) |>
    janitor::clean_names(case="all_caps")

  ## loop through all the channels and create a biexpTrans object per channel
  biexpTrans <-
    purrr::pmap(list(
      p=tbl_biexp_params$POSITIVE_DEC,
      n=tbl_biexp_params$EXT_NEG_DEC,
      wb=tbl_biexp_params$WIDTH_BASIS
    ), function(p, n, wb, ...) {
      flowWorkspace::flowjo_biexp_trans(
        pos=p,
        neg=n,
        widthBasis=wb
      )
    })


  ## Give names
  names(biexpTrans) <- tbl_biexp_params$FULL_NAME

  ## Create transformerList in order to apply to gs
  trans <- flowWorkspace::transformerList(names(biexpTrans), biexpTrans)

  ## Apply to gs
  gs <- flowWorkspace::transform(gs, trans)

  return(gs)
}
