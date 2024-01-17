#' Generate the Dunedin Methylation Pace of Aging Scores!
#'
#' \code{PACEProjector} returns the Dunedin Pace of Aging Methylation Scores
#'
#' @param betas A numeric matrix containing the percent-methylation for each probe.  Missing data should be 'NA's.  The rows should be probes, with the probe ID as the row name, and the columns should be samples, with sample names as the column name.
#' @param proportionOfProbesRequired (default: 0.7).  This value specificies the threshold for missing data (see description for more details on how missing data is handled)
#' @return A list of mPACE values.  There will be one element in the list for each mPACE model.  Each element will consist of a numeric vector with mPACE values.  The names of the values in the vector will be the sample names from the 'betas' matrix. The output is equivalent to the pace of aging (in years) expected over a 1-year period.
#' @details This function returns the Dunedin Methylation Pace of Aging scores for methylation data generated from either the Illumina 450K array or the Illumina EPIC array.  The Age45 score is one that has been trained on data based on 3 waves of collection (26, 38, and 45).  The manuscript is currently in preparation, but has been shown to be more accurate than the Age38 score.
#' Missing data handled in two different ways (and the threshold for both is set by the 'proportionOfProbesRequired' parameter).  First, if a sample is missing data for more probes than the threshold, the sample will get an NA back for a score.  If a particular probe is missing fewer samples than the threshold, then missing data is set to the mean in the provided 'betas' matrix.  If a probe is missing more samples than the threshold, then all samples in the 'betas' matrix have their value replaced with the mean of the training data for that particular model.
#' Because of how we handle missing data, it is recommended that entire cohorts be run at once as a large 'betas' matrix.
#' @examples
#' PACEProjector(betas)


PACEProjector = function( betas, proportionOfProbesRequired=0.7 ) {
  requireNamespace("preprocessCore")

  # Check if row names have more than 10 characters and rename if necessary
  if (any(nchar(rownames(betas)) > 10)) {
    rownames(betas) <- substr(rownames(betas), 1, 10)
  }
  
 if( any( grepl( "TC11", rownames(betas)))) {
      print("This looks like EPICv2 array data. If EPICv2, DunedinPACE will lower the proportion of probes required to from 0.8 to 0.7 and proceed with some missing probes.")
  }
  
  # loop through models
  model_results <- lapply(mPACE_Models$model_names, function(model_name) {
    # make sure it has been converted to a matrix
    if( !is.numeric(as.matrix(betas)) ) { stop("betas matrix/data.frame is not numeric!") }
    probeOverlap <- length(which(rownames(betas) %in% mPACE_Models$model_probes[[model_name]])) / length(mPACE_Models$model_probes[[model_name]])
    probeOverlap_background <- length(which(rownames(betas) %in% mPACE_Models$gold_standard_probes[[model_name]])) / length(mPACE_Models$gold_standard_probes[[model_name]])
    # make sure enough of the probes are present in the data file
    if( probeOverlap < proportionOfProbesRequired | probeOverlap_background < proportionOfProbesRequired ) {
      result <- rep(NA, ncol(betas))
      names(result) <- colnames(betas)
      result
    } else {
      # Work with a numeric matrix of betas
      betas.mat <- as.matrix(betas[which(rownames(betas) %in% mPACE_Models$gold_standard_probes[[model_name]]),])
      # If probes don't exist, we'll add them as rows of values based on their mean in the gold standard dataset
      probesNotInMatrix <- mPACE_Models$gold_standard_probes[[model_name]][which(mPACE_Models$gold_standard_probes[[model_name]] %in% rownames(betas.mat) == F)]
      if( length(probesNotInMatrix) > 0 ) {
        for( probe in probesNotInMatrix ) {
          tmp.mat <- matrix(0, nrow=1, ncol=ncol(betas.mat))
          rownames(tmp.mat) <- probe
          colnames(tmp.mat) <- colnames(betas.mat)
          tmp.mat[probe,] <- rep(mPACE_Models$gold_standard_means[[model_name]][probe], ncol(tmp.mat))
          betas.mat <- rbind(betas.mat, tmp.mat)
        }
      }

      # Identify samples with too many missing probes and remove them from the matrix
      samplesToRemove <- colnames(betas.mat)[which(apply(betas.mat, 2, function(x) { 1 - ( length(which(is.na(x))) / length(x) ) < proportionOfProbesRequired}))]
      if( length(samplesToRemove) > 0 ) {
        betas.mat <- betas.mat[,-which(colnames(betas.mat) %in% samplesToRemove)]
      }
      if(ncol(betas.mat) > 0) {
        # Identify missingness on a probe level
        pctValuesPresent <- apply( betas.mat, 1, function(x) { 1 - (length(which(is.na(x))) / length(x)) } )
        # If they're missing values, but less than the proportion required, we impute to the cohort mean
        probesToAdjust <- which(pctValuesPresent < 1 & pctValuesPresent >= proportionOfProbesRequired)
        if( length(probesToAdjust) > 0 ) {
          if( length(probesToAdjust) > 1 ) {
            betas.mat[probesToAdjust,] <- t(apply( betas.mat[probesToAdjust,], 1 , function(x) {
              x[is.na(x)] = mean( x, na.rm = TRUE )
              x
            }))
          } else {
            betas.mat[probesToAdjust,which(is.na(betas.mat[probesToAdjust,]))] <- mean(betas.mat[probesToAdjust,], na.rm=T)
          }
        }
        # If they're missing too many values, everyones value gets replaced with the mean from the Dunedin cohort
        if( length(which(pctValuesPresent < proportionOfProbesRequired)) > 0 ) {
          probesToReplaceWithMean <- rownames(betas.mat)[which(pctValuesPresent < proportionOfProbesRequired)]
          for( probe in probesToReplaceWithMean ) {
            betas.mat[probe,] <- rep(mPACE_Models$model_means[[model_name]][probe], ncol(betas.mat))
          }
        }

        # Normalize the matrix to the gold standard dataset
        betas.norm <- preprocessCore::normalize.quantiles.use.target(betas.mat, target=mPACE_Models$gold_standard_means[[model_name]])
        rownames(betas.norm) <- rownames(betas.mat)
        colnames(betas.norm) <- colnames(betas.mat)
        # Calculate score:
        score = mPACE_Models$model_intercept[[model_name]] + rowSums(t(betas.norm[mPACE_Models$model_probes[[model_name]],]) %*% diag(mPACE_Models$model_weights[[model_name]]))
        names(score) <- colnames(betas.norm)
        if( length(samplesToRemove) > 0 ) {
          score.tmp <- rep(NA, length(samplesToRemove))
          names(score.tmp) <- samplesToRemove
          score <- c(score, score.tmp)
        }
        score <- score[colnames(betas)]
        score
      } else {
        result <- rep(NA, ncol(betas.mat))
        names(result) <- colnames(betas.mat)
        result
      }
    }
  })
  names(model_results) <- mPACE_Models$model_names
  model_results
}



# Note: This is a function to print an error for PoAmProjector as a result of the name change of the function.

PoAmProjector <- function(beta_matrix) {
  .Defunct(msg = "Function `PoAmProjector()` has been replaced by `PACEProjector()`...\nUse `PACEProjector()` instead")
}
