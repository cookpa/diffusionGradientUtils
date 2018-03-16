## Bunch of general use functions for dealing with schemes


#
# Maps b-values to shell IDs, given shells, eg shells = c(0,1000,2000)
#
# It's expected that the actual b-values from the scanner will differ slightly from the nominal values
# but assumed that the closest one will still be correct
#
bvalsToShells <- function(bvals, shells) {
    
    bvalShell = rep(0, length(bvals))
    
    for (i in 1:length(bvals)) {
        
        dist = abs(bvals[i] - shells)
        
        bvalShell[i] = which(dist == min(dist))
        
    }
    
    return(bvalShell)
    
}


## Of all shells, excluding b=0
## eg electrostaticEnergy(bvecs, bvalsToShells(bvals, c(0,1000)))
##
## The first shell is assumed to be zero or close to zero, and is ignored
##
## If no shell IDs are provided, all non-zero gradient directions are included
##
electrostaticEnergy <- function(bvecs, shellIDs = NA) {

    numMeas = nrow(bvecs)

    if (missing(shellIDs)) {
      shellIDs = vector("numeric", numMeas)
      
      bvecMod = apply(bvecs, 1, function(x) { sqrt(sum(x^2)) })

      shellIDs = 1 + (bvecMod > 0)
    }
    
    energy = matrix(nrow = numMeas, ncol = numMeas)

    # Do lower triangle, then fill in upper triangle for speed  
    for (i in 1:numMeas) {
        for (j in 1:i) {
            if (i == j || shellIDs[i] == 1 || shellIDs[j] == 1) {
                energy[i,j] = 0
            }
            else {

                dot = sum(bvecs[i,] * bvecs[j,])
                
                # Actual energy would be proportional to this for equal charge pairs
                energy[i,j] = 2 * sqrt(2) * ( 1 / sqrt(1 - dot) + 1 / sqrt(1 + dot) )
            }
        }
    }
    
    # Do lower triangle, then fill in upper triangle for speed  
    for (i in 1:numMeas) {
        for (j in i:numMeas) {
            energy[i,j] = energy[j,i]
        }
    }
    
    return(energy)
}


#
# Read a bvals file without annoyances of warnings or extra fields from trailing whitespace
#
readBvals <- function(bvalFile) {
  
  valsString = readLines(bvalFile, warn = F)[[1]]

  valsString = trimws(valsString)

  tokens = unlist(strsplit(valsString, " "))

  return(as.numeric(tokens))
}


#
# Read a bvecs file and transposes it into the more readable Nx3 format expected by the other functions
#
readBvecs <- function(bvecFile) {
  
  bString = readLines(bvecFile, warn = F)

  bString = trimws(bString)

  tokens = lapply(bString, strsplit, " ")

  bvecs = c(tokens[[1]][[1]], tokens[[2]][[1]], tokens[[3]][[1]])

  bvecs = as.numeric(bvecs)
  
  dim(bvecs) = c(length(tokens[[1]][[1]]), 3)

  return(bvecs)
  
}

## Vector dot product
dot <- function(vec1, vec2) {

  sum = 0

  dim = length(vec1)

  if (length(vec2) != dim) {
    stop("Vectors are not of the same dimension")
  }
  
  for (i in 1:dim) {
    sum = sum + vec1[i] * vec2[i]
  }

  return(sum)
}

## 3D Vector cross product
cross <- function(v1,v2) {
  
  if (length(v1) != length(v2) || length(v1) != 3) {
    stop("Vectors must be 3D")
  }

  result = matrix(c(v1[2] * v2[3] - v2[2] * v1[3], v2[1] * v1[3] - v1[1] * v2[3], v1[1] * v2[2] - v2[1] * v1[2]), ncol = 3)

  return(result)
  
}

normalizePoints <- function(vecs) {

  uvecs = t(apply(vecs, 1, function(x) { x / sqrt(sum(x^2)) }))

  return(uvecs)
}

##
## Returns the number of measurements on each shell for a given bvals vector 
## 
##
shellMeasurements <- function(bvals, shells) {

  numShells = length(shells)
  
  shellVals = bvalsToShells(bvals, shells)

  meas = vector("numeric", numShells)

  for (i in 1:numShells) {
    meas[i] = length(which(shellVals == i))
  }

  return(meas)
  
}


readCaminoElecPoints <- function(numPoints, caminoDir = "/data/jet/pcook/code/camino") {

  numPointsString = sprintf("%03d", numPoints)

  file = file.path(caminoDir, "PointSets", paste("Elec", numPointsString, ".txt", sep = ""))

  if (!file.exists(file)) {
    stop(paste("No point set available for", numPoints, "points"))
  }

  points = read.table(file, header = FALSE, sep = "\n")

  numPointsInFile = points[1,1]

  if (numPointsInFile != numPoints) {
    stop("Point set on disk contains the wrong number of points")
  }

  points = points[2:nrow(points),1]
  
  dirs = matrix(points, ncol = 3, byrow = T)

  if (nrow(dirs) != numPoints) {
    stop("Point set on disk contains the wrong number of points")
  }
  
  return(dirs)
  
}
