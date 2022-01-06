#' See the list of probes needed for each Dunedin Pace of Aging Methylation Score Model
#'
#' \code{getRequiredProbes} returns a list containing the IDs of the probes needed for each Dunedin Pace of Aging Methylation Score Model.
#' @param backgroundList (default: FALSE).  If set to true, it will return the full set of 20,000 probes in the background set used for normalization
#' @return A list of string vectors.  Each element in the list is a different model, and the character vector contains the IDs needed for that particular model.
#' @examples
#' getRequiredProbes()
#'
getRequiredProbes = function(backgroundList=F) {
  if( backgroundList == FALSE ) {
    mPOA_Models$model_probes
  } else {
    mPOA_Models$gold_standard_probes
  }
}






