##
## Compares the energy of a point set to an optimized point set with the same N, perturbations of an optimized
## point set, and random point sets. Any reasonable point set should have much lower energy than a random set.
##
##
evaluatePointSetCoverage <- function(points, referencePoints, plot = F, iterations = 1000) {

    numPoints = nrow(points)

    energy = sum(electrostaticEnergy(points))
    
    refEnergy = sum(electrostaticEnergy(referencePoints))

    perturbedEnergy = rep(0,iterations)

    ## Mean and SD of perturbation, degrees
    perturbationM = 0

    ## Test a few sigmas
    perturbationSigma = c(1,2,5)

    zEnergyUnderPerturbation = rep(0, length(perturbationSigma))
    pEnergyUnderPerturbation = rep(0, length(perturbationSigma))
    perturbedEnergy = matrix(nrow = iterations, ncol = length(perturbationSigma))
    
    ## Compare to perturbed elec set
    for (n in 1:length(perturbationSigma)) {
        
        for (i in 1:iterations) {
            perturbedEnergy[i,n] = sum(electrostaticEnergy(perturbPoints(referencePoints, thetaSD = perturbationSigma[n], thetaMean = perturbationM)))
        }
        
        meanPE = mean(perturbedEnergy[,n])
        sdPE = sd(perturbedEnergy[,n])
        
        # p of seeing an energy less than the point set energy by perturbation of the elec point set
        pEnergyUnderPerturbation[n] = pnorm(energy, mean = meanPE, sd = sdPE, lower.tail = T)
        
        zEnergyUnderPerturbation[n] = (energy - meanPE) / sdPE
        
    }

    ## Other relevant info: min angle between pairs of points, median angle between pairs of points
    
    return(list(pointSetEnergy = energy, referenceEnergy = refEnergy, perturbedEnergy = perturbedEnergy, perturbationZ = zEnergyUnderPerturbation, perturbationPVal = pEnergyUnderPerturbation, perturbationSigma = perturbationSigma))
   

    ## Compare to random point set

    ## Any reasonable point set should do much better than random
    

    ## Return some results
    
    
}


##
## Compute a random point set
##
## These point sets are NOT optimal and should definitely not be used for data acquisition
##
## minAngle - reject samples that are closer than this angle (in degrees) to each other. This
##            reduces outlier energy terms, but don't make this too large or it will take a long time
##            to compute acceptable samples
##
getRandomPointSet <- function(numPoints, minAngle = 0) {
  

    vecs = matrix(rnorm(3 * numPoints,0,1), nrow = numPoints, ncol = 3, byrow = T)
    
    vecs = t(apply(vecs, 1, function(x) { x / sqrt(sum(x^2)) }))

    if (minAngle > 0) {

        maxDot = cos(minAngle * pi / 180)

        repeat {
            
            absDot = matrix(rep(0,numPoints*numPoints), ncol = numPoints, nrow = numPoints)

            ## Upper triangle
            for (i in 1:(numPoints-1)) {
                for (j in (i+1):numPoints) {
                    absDot[i,j] = abs(dot(vecs[i,],vecs[j,]))
                }
            }
            
            closePoints = which(absDot > maxDot, arr.ind = T)
            
            if(length(closePoints) == 0) {
                break
            }

            replace = unique(sort(closePoints))

            numReplacements = length(replace)

            replacementVecs = matrix(rnorm(3 * numReplacements,0,1), nrow = numReplacements, ncol = 3, byrow = T)

            replacementVecs = t(apply(replacementVecs, 1, function(x) { x / sqrt(sum(x^2)) }))

            vecs[replace,] = replacementVecs
            
        }

    }
    
    return(vecs)

}


## Get a rotation matrix about some axis
## Axis is a unit vector, theta is in radians
getRotationMatrix <- function(axis, theta) {

    ct = cos(theta);
    st = sin(theta);
    oMCT = 1 - ct
    
    axis = axis / sqrt(sum(axis^2))
    
    rx = axis[1]
    ry = axis[2]
    rz = axis[3]
  
    rotation = matrix(nrow = 3, ncol = 3)
    
    rotation[1,1] = rx*rx*oMCT + ct
    rotation[2,1] = ry*rx*oMCT + rz*st
    rotation[3,1] = rz*rx*oMCT - ry*st
    rotation[1,2] = rx*ry*oMCT - rz*st
    rotation[2,2] = ry*ry*oMCT + ct
    rotation[3,2] = rz*ry*oMCT + rx*st
    rotation[1,3] = rx*rz*oMCT + ry*st
    rotation[2,3] = ry*rz*oMCT - rx*st
    rotation[3,3] = rz*rz*oMCT + ct
    
    return(rotation)
    
}


## perturb an (Nx3) point set by applying random rotations to each point
##
## The rotation axes are random vectors in the plane perpendicular to each point
##
## The rotation angle is drawn from a normal distribution with user defined parameters
##
perturbPoints <- function(points, thetaSD = 1, thetaMean = 0) {

  numPoints = nrow(points)
  
  perturbedPoints = points

  for (i in 1:numPoints) {

    point = matrix(points[i,], ncol = 3)

    ## Cross this vector with the z axis, unless the point is too close to z
    crossAxis = matrix(c(0,0,1), ncol = 3)
    
    if (abs(dot(point, crossAxis)) > 0.98) {
      crossAxis = matrix(c(1,0,0), ncol = 3)
    }

    ## This now perpendicular to the point and to the crossAxis
    norm = cross(point, crossAxis)

    ## Now randomly rotate this axis about crossAxis, to get our final rotation axis
    axisRotMat = getRotationMatrix(point, dunif(1, 0, 2 * pi))

    perturbationAxis = axisRotMat %*% t(norm)

    rotation = getRotationMatrix(perturbationAxis, rnorm(1,thetaMean,thetaSD) * pi / 180)

    perturbedPoints[i,] = rotation %*% t(point)
    
  }

  return(perturbedPoints)
  
}

