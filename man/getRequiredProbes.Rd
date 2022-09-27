\name{getRequiredProbes}
\alias{getRequiredProbes}
\title{
%%   ~~ data name/kind ... ~~
View key probes involved in calculating DunedinPACE
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Returns a list containing the IDs of the probes needed for each Dunedin Pace of Aging Methylation Score Model. By default does not include all the required probes (20k) used in the normalization and DunedinPACE calculation.
}
\usage{getRequiredProbes(backgroundList=F)}
\arguments{
  \item{backgroundList}{Default = FALSE. If set to true, it will return the full set of 20,000 probes in the background set used for normalization}
}

\details{
%%  ~~ If necessary, more details than the __description__ above ~~

The probes derived from this function are NOT sufficient to calculate DunedinPACE. If you subset using this list, you will get `NA` from the `PoAmProjector()` function.

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Belsky DW, Caspi A, Corcoran DL, Sugden K, Poulton R, Arseneault L, Baccarelli A, Chamarti K, Gao X, Hannon E, et al. 2022. DunedinPACE, a DNA methylation biomarker of the pace of aging. Deelen J, editor. eLife. 11:e73420. doi:10.7554/eLife.73420.

}
\examples{
getRequiredProbes()
}
\keyword{file}