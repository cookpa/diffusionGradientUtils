## Functions to import INRIA multi-shell gradient dirs from
## http://www.emmanuelcaruyer.com/q-space-sampling.php
##
## Input format is
## shell,u_x,u_y,u_z


##
## inriaDirs - INRIA format directions with shell IDs, with columns shellID,u_x,u_y,u_z for vector u
##
## shells - vector of b-values corresponding to the shell IDs, eg c(0, 1000, 2000), where shells[i] is the b-value
## for a measurement with shellID == i
##
##
## Returns a list containing bvals and bvecs
##
inriaShellsToBvecsBvals <- function(inriaDirs, shells) {

  numMeas = nrow(inriaDirs)

  if (length(unique(inriaDirs[,1])) != length(shells)) {
    stop("Number of shells in input does not match length of shell b-value arg. Cannot translate shell ID to b-value")
  }

  bvecs = matrix(ncol = 3, nrow = numMeas)
  bvals = vector("numeric", length = numMeas)
  mod = vector("numeric", length = numMeas)

  for (i in 1:numMeas) {
    gradDir = as.numeric(inriaDirs[i,2:4])

    mod[i] = sqrt(gradDir[1]^2 + gradDir[2]^2 + gradDir[3]^2)
    
    if (mod[i] > 0) {
      bvecs[i,] = gradDir / mod[i]
    }
    else {
      bvecs[i,] = gradDir
    }
 
    bvals[i] = shells[inriaDirs[i,1]]
    
  }
  
  fslScheme = list(bvals = bvals, bvecs = bvecs, originalVectorMagnitude = mod)
 
  return(fslScheme)
  
}

