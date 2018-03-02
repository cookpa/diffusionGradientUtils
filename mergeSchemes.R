##
## Merge schemes by interleaving measurements. Useful to merge different DWI shells or
## insert b=0 measurements at regular intervals.
##
## Pass the longest scheme first, the function will merge the smaller one into it.
## If length(bvals) > length(refBvals), the two are switched.
##
## If the option offsetFromEnd is TRUE, the final measurement in the merged scheme is the
## final measurement in bvecs / bvals. Otherwise, the first measurement in the merged scheme
## is the first measuremeng in bvecs / bvals. 
##
mergeSchemes <- function(refBvals, refBvecs, bvals, bvecs, offsetFromEnd = T) {

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
  if (offsetFromEnd) {
    insertionIndices = round(seq(interval,mergedNumMeas, interval))
  }
  else {
    insertionIndices = round(seq(1, mergedNumMeas - interval + 1, interval))
  }

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
