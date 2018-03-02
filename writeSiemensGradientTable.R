
##
## Writes a Siemens style gradient information file 
## 
##  scaledGradDirs = writeSiemensGradientTable(bvals, bvecd, fileName)
##
writeSiemensGradientTable <- function(bvals, bvecs, fileName) {

  numMeas = nrow(bvecs)

  if (nrow(bvecs) != length(bvals)) {
    stop("bvecs and bvals are mismatched")
  }

  maxB = max(bvals)

  shellB = sort(unique(bvals))

  numShells = length(shellB)

  shellIDs = bvalsToShells(bvals, shellB)

  numMeasPerShell = vector("numeric", numShells);

  for (i in 1:numShells) {
    numMeasPerShell[i] = length(which(shellIDs == i)) 
  }

  bFrac = shellB / maxB

  out = file(fileName, open = "w")

  ## Write some info for the table file
  writeLines("## Multi-shell gradient directions with vectors scaled by sqrt(b / bMax)", out)
  writeLines(c("##", "## Shell Num Meas\tb/bMax"), out, sep = "\n")

  for (i in 1:numShells) {
    writeLines(sprintf("## %d\t\t%.3f", numMeasPerShell[i],bFrac[i]), out)
  }

  writeLines(c("##", sprintf("[directions=%d]", numMeas), "CoordinateSystem = xyz", "Normalisation = none"), out, sep = "\n")

  
  ## write format Vector[0] = ( 0.000, 0.000, 0.000 )

  scaledVecs = matrix(nrow = numMeas, ncol = 3)
  
  for (i in 1:numMeas) {
    scaledVecs[i,] = bvecs[i,] * sqrt(bvals[i] / maxB)
    vecLine = sprintf("Vector[%d] = ( %.6f,  %.6f,  %.6f )", i-1, scaledVecs[i,1], scaledVecs[i,2], scaledVecs[i,3])
    writeLines(vecLine, out)
  }
 
  close(out)

  ## Return some information including the set of scaled directions so they can be checked
  scaledGradMod = apply(scaledVecs, 1, function(x) { sqrt(sum(x^2)) } )
  
  return(list(gradientDirections = scaledVecs, gradientMagnitude = scaledGradMod))
  
}
