\name{PACEProjector}
\alias{PACEProjector}
\title{
%%   ~~ data name/kind ... ~~
Calculate DunedinPACE
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This function returns the Dunedin Methylation Pace of Aging scores for methylation data generated from either the Illumina 450K array or the Illumina EPIC array.

Missing data handled in two different ways (and the threshold for both is set by the 'proportionOfProbesRequired' parameter). First, if a sample is missing data for more probes than the threshold, the sample will get an NA back for a score.  If a particular probe is missing fewer samples than the threshold, then missing data is set to the mean in the provided 'betas' matrix.  If a probe is missing more samples than the threshold, then all samples in the 'betas' matrix have their value replaced with the mean of the training data for that particular model.

}
\usage{PACEProjector(betas, proportionOfProbesRequired = 0.8)}

\arguments{
  \item{betas}{A numeric matrix containing the percent-methylation for each probe.  Missing data should be 'NA's.  The rows should be probes, with the probe ID as the row name, and the columns should be samples, with sample names as the column name.}
  \item{proportionOfProbesRequired}{Default = 0.8.  This value specificies the threshold for missing data (see description for more details on how missing data is handled). Note: that if the function detects EPICv2 data based on rownames that contain the replicate suffix, the proportion of probes required will automatically be lowered to 0.7.}
}

\details{
%%  ~~ If necessary, more details than the __description__ above ~~
A list of mPACE values.  There will be one element in the list for each mPACE model.  Each element will consist of a numeric vector with mPACE values.  The names of the values in the vector will be the sample names from the 'betas' matrix. The output is equivalent to the pace of aging (in years) expected over a 1-year period.

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Belsky DW, Caspi A, Corcoran DL, Sugden K, Poulton R, Arseneault L, Baccarelli A, Chamarti K, Gao X, Hannon E, et al. 2022. DunedinPACE, a DNA methylation biomarker of the pace of aging. Deelen J, editor. eLife. 11:e73420. doi:10.7554/eLife.73420.

}
\examples{
\dontrun{PACEProjector(betas)}
}
\keyword{file}
