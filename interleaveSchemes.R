##
## Interleave schemes. Useful to merge shells or insert b=0 measurements at regular intervals
##
## Pass the longest scheme first, and merge the smaller one. If length(bvals) > length(refBvals), the two
## are switched
##
interleaveSchemes <- function(refBvals, refBvecs, bvals, bvecs) {

  ## Usual sanity checks
  if ( nrow(bvecs) != length(bvals) ) {
    stop("Number of bvals and bvecs do not match")
  }
  if ( nrow(refBvecs) != length(refBvals) ) {
    stop("Number of refBvecs and refBvals do not match")
  }

  ## merge the shorter scheme into the longer one
  if (length(bvals) > length(refBvals)) {
    tmpBvals = bvals
    bvals = refBvals
    refBvals = tmpBvals

    tmpBvecs = bvecs
    bvecs = refBvecs
    refBvecs = tmpBvecs

  }

  numMeasRef = length(refBvals)

  numMeasToInsert = length(bvals)

  mergedNumMeas = numMeasRef + numMeasToInsert
  
  interval = mergedNumMeas / numMeasToInsert
  
  ## Where to insert measurements relative to the refBvecs / refBvals
  insertionIndices = round(seq(interval,mergedNumMeas, interval))

  ## Sanity checks
  if (length(insertionIndices) != numMeasToInsert) {
    stop(paste("insertion indices has length ", length(insertionIndices), " but number of measurements to interleave is ", numMeasToInsert, sep = ""))
  }

  if (length(unique(insertionIndices)) != numMeasToInsert) {
     stop("Duplicate values in insertion indices")
  }

  mergedBvecs = matrix(nrow = mergedNumMeas, ncol = 3)

  mergedBvecs[insertionIndices,] = bvecs
  mergedBvecs[-insertionIndices,] = refBvecs

  mergedBvals = vector("numeric", length = mergedNumMeas)

  mergedBvals[insertionIndices] = bvals

  mergedBvals[-insertionIndices] = refBvals

  ## Final sanity checks
  if (! all(mergedBvecs[insertionIndices,] == bvecs) ) {
    stop("ERROR bvecs not correct in output")
  }

  if (! all(mergedBvals[insertionIndices] == bvals) ) {
    stop("ERROR bvals not correct in output")
  }

  if (! all(mergedBvecs[-insertionIndices,] == refBvecs) ) {
    stop("ERROR refBvecs not correct in output")
  }

  if (! all(mergedBvals[-insertionIndices] == refBvals) ) {
    stop("ERROR refBvals not correct in output")
  }

  
  return(list(bvecs = mergedBvecs, bvals = mergedBvals))
  
}
