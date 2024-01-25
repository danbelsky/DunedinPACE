#' Generate the Dunedin Methylation Pace of Aging Scores!
#'
#' The \code{PoAmProjector} function has been deprecated! Please use \code{PACEProjector} instead.
#'
#' @param betas A numeric matrix containing the percent-methylation for each probe.  Missing data should be 'NA's.  The rows should be probes, with the probe ID as the row name, and the columns should be samples, with sample names as the column name.
#' @param proportionOfProbesRequired (default: 0.8).  This value specificies the threshold for missing data (see description for more details on how missing data is handled)
#' @return A list of mPoA values.  There will be one element in the list for each mPoA model.  Each element will consist of a numeric vector with mPoA values.  The names of the values in the vector will be the sample names from the 'betas' matrix.
#' @details This function returns the Dunedin Methylation Pace of Aging scores for methylation data generated from either the Illumina 450K array or the Illumina EPIC array.  The Age45 score is one that has been trained on data based on 3 waves of collection (26, 38, and 45).  The manuscript is currently in preparation, but has been shown to be more accurate than the Age38 score.
#' Missing data handled in two different ways (and the threshold for both is set by the 'proportionOfProbesRequired' parameter).  First, if a sample is missing data for more probes than the threshold, the sample will get an NA back for a score.  If a particular probe is missing fewer samples than the threshold, then missing data is set to the mean in the provided 'betas' matrix.  If a probe is missing more samples than the threshold, then all samples in the 'betas' matrix have their value replaced with the mean of the training data for that particular model.
#' Because of how we handle missing data, it is reccomended that entire cohorts be run at once as a large 'betas' matrix.
#' @examples
#' PoAmProjector(betas)

PoAmProjector = function(betas, proportionOfProbesRequired=0.8 ) {
  .Deprecated("PACEProjector()")
}

