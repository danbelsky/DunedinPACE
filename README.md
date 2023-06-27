# DunedinPACE
Pace of Age calculator for Illumina methyl-array data (DunedinPACE)

PACEProjector.R -- 20220105

Given a set of methylation beta values, this tool will calculate the Dunedin Methylation Pace of Aging Methylation Score (DunedinPACE)

#### Installation (via devtools):
```r
devtools::install_github("danbelsky/DunedinPACE")
```

#### Requirements (preprocesscore Bioconductor package):
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore")
```

#### Usage:
```r
library("DunedinPACE")
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

