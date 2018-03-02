
##
## Checks if a set of bvecs and bvals is the same as a reference set, to within some margin of error.
##
## The number of measurements must be the same between the schemes. The bvecs should be unit vectors
##
## bvecs - bvectors to be tested
## bvals - bvalues to be tested
## refBvecs - reference scheme bvecs
## refBvals - reference scheme bvals
## shells - optionally specify shells, and the acceptable varation in b-values for each. For example,
## a variation of 50 s / mm^2 might be acceptable at b=3000, but not at b=1000.
## epsB - acceptable variation in b-values, either a constant or one value per shell
## epsOrient - acceptable variation in bvec orientation, in degrees
##
## Returns a list containing, for each measurement:
##
## bvecAngle - angle in degrees between input and reference bvec
##
## bvalDiff - difference between input and reference bvalue
##
## schemesConsistent - TRUE if none of the measurements differ in bvalue or orientation
##                     by more than the allowed amount
## 
compareDiffusionSchemes <- function(bvals, bvecs, refBvals, refBvecs, shells = NA, epsB = 50, epsOrient = 0.5) {

  numMeas = nrow(refBvecs)
  
  shellIDs = 0
  
  if (! is.na(shells[1]) ) {
    shellIDs = bvalsToShells(bvals, shells)
  }
  else {
    shells = 1
    shellIDs = rep(1,numMeas)
  }
  
  if (length(epsB) != length(shells)) {
    epsB = rep(epsB, length(shells))
  }
  
  result = list(bvecAngle = vector("numeric", numMeas), bvalDiff = vector("numeric", numMeas), schemesConsistent = F)
  
  if (length(refBvals) != numMeas) {
    stop("Reference bvals and bvecs are mismatched");
  }
  if (nrow(bvecs) != numMeas) {
    stop("Input scheme does not have the same number of measurements as reference scheme");
  }
  if (length(bvals) != numMeas) {
    stop("Input bvals and bvecs are mismatched");
  }
  
  result$schemesConsistent = T
  
  ## We require the bvecs to be unit vectors. This should be the case with bvecs from dcm2nii
  for (i in 1:numMeas) {

    scalarProd = dot(bvecs[i,], refBvecs[i,])

    if (abs(scalarProd) - 1 > 1e-5) {
      stop("Invalid dot product, are gradient directions normalized?")
    }
    if (scalarProd > 1) {
      scalarProd = 1
    }
    if (scalarProd < -1) {
      scalarProd = -1
    }
    if (scalarProd == 0 && sum(c(bvecs[i,], refBvecs[i,]) == 0)) {
      scalarProd = 1
    }
    
    result$bvecAngle[i] = 180 * acos(scalarProd) / pi

    if (result$bvecAngle[i] > epsOrient) {
      result$schemesConsistent = F
    }

    result$bvalDiff[i] = bvals[i] - refBvals[i]
    
    if (result$bvalDiff[i] > epsB[shellIDs[i]]) {
      result$schemesConsistent = F
    }
    
  }

  return(result)
  
}

