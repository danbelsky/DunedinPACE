# DunedinPACE
Pace of Age calculator for Illumina methyl-array data (DunedinPACE)

PoAmProjector.R -- 20220105

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
PoAmProjector(betas)
```

#### To see list of probes necessary for each model:
```r
getRequiredProbes()
```
Note: In order to calculate `DunedinPACE`, you will need many other probes used in the pre-processing steps prior to `DunedinPACE` selection. Do __not__ filter your beta matrix using `getRequiredProbes()` prior to using the `PoAmProjector()` function. 


## Input:
####  betas:
    Matrix or data.frame of beta values where rownames are probe ids and column names should correspond to sample names.
    Ensure beta values are numeric and that missing values should be coded as 'NA'

####  proportionOfProbesRequired:
    This is the proportion of probes to have a non-missing value for both the sample to have
    DunedinPACE calculated, as well as to determine if we can impute the mean from the current cohort
    By default, this is set to 0.8

## Output:
   A list containing the mPoAs for each model

