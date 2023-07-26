# News

#### DunedinPACE now works on EPICv2! 

Follow instructions and formatting as for EPICv2 or 450k (below). If you have performed probe QC prior to DunedinPACE, you may need to lower the `proportionOfProbesRequired = 0.8` argument in the function, as EPICv2 is already missing some probes that were originally in EPICv1. 

# Introduction to DunedinPACE
DunedinPACE is a novel DNA-methylation based blood biomarker of the pace of aging for gerontology and geroscience. It shows high test-retest reliability, is associated with morbidity, disability, and mortality, and indicated faster aging in young adults with childhood adversity. DunedinPACE effect-sizes are similar to GrimAge Clock effect-sizes. In some analyses of incident morbidity, disability, and mortality, DunedinPACE added incremental prediction beyond GrimAge.

Citation for the original paper describing DunedinPACE in more detail:

`Belsky DW, Caspi A, Corcoran DL, Sugden K, Poulton R, Arseneault L, Baccarelli A, Chamarti K, Gao X, Hannon E, et al. 2022. DunedinPACE, a DNA methylation biomarker of the pace of aging. Deelen J, editor. eLife. 11:e73420. doi:10.7554/eLife.73420.`

# The DunedinPACE R package
The DunedinPACE R package allows users to calculate DunedinPACE in their own DNA methylation dataset. Given a set of methylation beta values, where rows are probes and columns are subjects, this tool will produce individual DunedinPACE.


#### Installation (via devtools):
To install the DunedinPACE and associated vignette, please use the following code:
```r
devtools::install_github("danbelsky/DunedinPACE", build_vignettes = TRUE)
```

#### Requirements (preprocesscore Bioconductor package):
DundeinPACE utilizes an internal normalization process that is built on the `preprocessCore` package. 
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore")
```

#### Usage:
Requires a matrix of beta DNA methylation values. Rows must be probes (CpGs) and columns must be subjects. Rownames should be Illumina probe names (i.e. cg#########). 
```r
# load the DunedinPACE package
library("DunedinPACE")

# Apply function to matrix of beta values (here called 'betas'). 
PACEProjector(betas)
```

#### To see list of probes necessary for each model:
```r
getRequiredProbes()
```

`getRequiredProbes(backgroundList = FALSE)` returns 173 probes used for calculating DunedinPACE directly. 
`getRequiredProbes(backgroundList = TRUE)` returns 173 probes used for calculating DunedinPACE directly, as well as 19827 probes used in the normalization process. 

We do not recommend excluding the 19827 probes used for normalization and calculating DunedinPACE using only the 173 DunedinPACE associated probes, as this will affect DunedinPACE estimates.

## Input:

####  betas:
    Matrix or data.frame of beta values where rownames are probe ids and column names should correspond to sample names.
    Ensure beta values are numeric and that missing values should be coded as 'NA'

####  proportionOfProbesRequired:
    This is the proportion of probes required to have a non-missing value for both the sample to have
    DunedinPACE calculated, as well as to determine if we can impute the mean from the current cohort.
    By default, this is set to 0.8

## Output:
   A vector of DunedinPACE estimates.      

