
##
## Get SNR as mean(signal) / sigma over all measurements
##
##
dwiSNR <- function(data, bvals, rois, sigma) {

  labelIndices = sort( unique( rois[rois > 0] ) )

  numLabels = length(labelIndices)

  numMeas = length(bvals)

  output = matrix(nrow = numLabels, ncol = numMeas)
  
  for (i in 1:numLabels) {

    roiData = data[rois == labelIndices[i]]

    dim(roiData) = c(length(roiData) / numMeas, numMeas)

    ## Here we are summing over voxels in a single DWI volume, not the same voxel over
    ## multiple b=0, so we mean cols not rows
    meanSignal = colMeans(roiData)

    output[i,] = meanSignal / sigma[i]
    
  }

  return(output)
  
  
}


##
## Get CNR as sd(signal) / sigma over all measurements
##
##
dwiCNR <- function(data, bvals, rois, sigma) {

  labelIndices = sort( unique( rois[rois > 0] ) )

  numLabels = length(labelIndices)

  numMeas = length(bvals)

  output = matrix(nrow = numLabels, ncol = numMeas)
  
  for (i in 1:numLabels) {

    roiData = data[rois == labelIndices[i]]

    dim(roiData) = c(length(roiData) / numMeas, numMeas)

    sdSignal = apply(sd, 2, roiData)

    output[i,] = sdSignal / sigma[i]
    
  }

  return(output)
  
  
}





## Estimate the SNR in b=0 in some ROIs, assuming the presence of multiple b0
##
## Uses the unweighted data (b <= unweightedB).
##
## See Dietrich et al JMRI 26:376-385 (2007)
## 
##
estimateB0SNR <- function(data, bvals, rois, unweightedB = 5) {

  b0Ind = which(bvals <= unweightedB)

  numB0 = length(b0Ind)

  if (numB0 < 2) {
    stop("Cannot compute sigma with fewer than 2 b=0 measurements")
  }
  
  labelIndices = sort( unique( rois[rois > 0] ) )

  numLabels = length(labelIndices)

  ## results are sigma Diff and sigma Mult for each ROI
  ##
  ## Sigma mult is NA if there are 2 b=0
  ##
  sigmaDiff = vector("numeric", length = numLabels)
  snrDiff = vector("numeric", length = numLabels)
  sigmaMult = vector("numeric", length = numLabels)
  snrMult = vector("numeric", length = numLabels)
  
  for (i in 1:numLabels) {

    labelVoxels = which(as.array(rois) == labelIndices[i], arr.ind = T)
    
    numLabelVoxels = nrow(labelVoxels)

    print(paste("Processing", numLabelVoxels, "voxels in label", labelIndices[i]))
    
    b0Data = matrix(nrow = numLabelVoxels, ncol = numB0)
    
    for (v in 1:numLabelVoxels) {
      b0Data[v,] = as.numeric(data[labelVoxels[v,1], labelVoxels[v,2], labelVoxels[v,3], b0Ind])
    }

    diff = snrDiff(b0Data)
    mult = snrMult(b0Data)

    snrDiff[i] = diff$SNR
    sigmaDiff[i] = diff$Sigma

    snrMult[i] = mult$SNR
    sigmaMult[i] = mult$Sigma
    
  }

  return(data.frame(list(LabelID = labelIndices, SNRDiff = snrDiff, SigmaDiff = sigmaDiff, SNRMult = snrMult, SigmaMult = sigmaMult)))
  
}



## snrMult(dataMat) where dataMat is a N x M matrix, for M b=0 measurements over N voxels
##
## Returns the SNR and sigma
## 
## Only defined for M > 2
##
snrMult <- function(dataMat) {

  result = list(SNR = NA, Sigma = NA) 
  
  numVox = nrow(dataMat)
  
  numB0 = ncol(dataMat)

  if (numB0 < 3) {
    return(result)
  }

  sumB0 = 0

  sumB0 = sum(colSums(dataMat))

  meanSD = mean(apply(dataMat, 1, sd))

  result$SNR = sumB0 / (numB0 * numVox * meanSD)
  result$Sigma = meanSD
  
  return(result)
  
}

## Computes the SNR using only the first two b=0 images
snrDiff <- function(dataMat) {

  numVox = nrow(dataMat)
  
  numB0 = ncol(dataMat)

  result = list(SNR = NA, Sigma = NA) 

  if (numB0 != 2) {
    if (numB0 < 2) {
      return(result)
    }
    else {
      dataMat = dataMat[,1:2]
    }
  }
  
  
  meanB0 = sum(colSums(dataMat)) / (2 * numVox)
  
  diffB0 = dataMat[,1] - dataMat[,2]

  sigmaDiff = ( 1 / sqrt(2) ) * sd(diffB0)
  
  result$SNR = meanB0 / sigmaDiff
  result$Sigma = sigmaDiff

  return(result)
}
