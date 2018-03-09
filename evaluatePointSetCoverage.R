##
## Compares the energy of a point set to an optimized point set with the same N and perturbations of those points.
## point set.
##
##
evaluatePointSetCoverage <- function(points, referencePoints, plot = F) {

    numPoints = nrow(points)

    if (nrow(referencePoints) != numPoints) {
        stop("Point sets have a different number of points")
    }

    energy = sum(electrostaticEnergy(points))
    
    refEnergy = sum(electrostaticEnergy(referencePoints))

    print(paste("energy =", energy))

    ## Mean and SD of perturbation, degrees
    perturbationM = 0

    ## Test a few sigmas
    perturbationSigma = c(2,4,8)
  
    ## Compare to perturbed elec set at various sigmas
    ## Want to find the smallest perturbation that is has higher average energy than the points being evaluated
    finalPerturbationSigma = -1
    
    counter = 1

    iterations = 100
    
    perturbedEnergy = rep(0,iterations)

    meanPE = 0
    sdPE = 0

    ## Do a quick search for the perturbation that produces a similar mean energy
    while(counter <= length(perturbationSigma) && finalPerturbationSigma < 0) {
        
        for (i in 1:iterations) {
            perturbedEnergy[i] = sum(electrostaticEnergy(perturbPoints(referencePoints, thetaSD = perturbationSigma[counter],
                                                                         thetaMean = perturbationM)))
        }
        
        meanPE = mean(perturbedEnergy)
        sdPE = sd(perturbedEnergy)

        print(paste("  Sigma =", perturbationSigma[counter], " |  meanPE =", meanPE, " |  sdPE =", sdPE))

        if (meanPE > energy || counter == length(perturbationSigma)) {
            finalPerturbationSigma = perturbationSigma[counter]
        }

        counter = counter + 1
    }

    ## Now do some more iterations for a nicer histogram
    iterations = 800
    perturbedEnergy = rep(0,iterations)

    
    for (i in 1:iterations) {
        perturbedEnergy[i] = sum(electrostaticEnergy(perturbPoints(referencePoints, thetaSD = finalPerturbationSigma,
                                                                   thetaMean = perturbationM)))
    }
    
    meanPE = mean(perturbedEnergy)
    sdPE = sd(perturbedEnergy)

    ## winsorize for display
    peQuantiles = quantile(perturbedEnergy, prob = seq(0,1,0.005))

    peWinsor = perturbedEnergy

    clip = peQuantiles[200]

    peWinsor[peWinsor > clip] = clip
    
        
    # p of seeing an energy less than the point set energy by perturbation of the elec point set
    pEnergyUnderPerturbation = pnorm(energy, mean = meanPE, sd = sdPE, lower.tail = T)
    
    zEnergyUnderPerturbation = (energy - meanPE) / sdPE
        
    ## Other relevant info: min angle between pairs of points
    minAngle = getClosestPointAngles(points)
    
    minRefAngle = getClosestPointAngles(referencePoints)

    if(plot) {
        xMin = min(refEnergy, energy)

        h = hist(peWinsor, n = 16, plot = F)
        
        xMax = max(energy, h$breaks[length(h$breaks)])

        hist(peWinsor, n = 16, xlim = c(xMin, xMax), xlab = "Electrostatic energy", main = paste("Energy compared to perturbed ref points (sigma = ", finalPerturbationSigma, ")", sep = ""))

        abline(v = energy, lty = 2)

        abline(v = refEnergy, lty = 3)

        legend("topright", legend = c("points", "ref points"), lty = c(2,3))

        plot.new()
        dev.new()

        plot(sort(minAngle, decreasing = T), type = "l", lty = 1, main = "Minimum angle between points", xlab = "Index", ylab = "angle (deg)")

        points(sort(minRefAngle, decreasing = T), type = "l", lty = 2)

        legend("bottomleft", legend = c("points", "ref points"), lty = c(1,2))
        
    }
    
    return(list(pointSetEnergy = energy, refEnergy = refEnergy, minAnglePoints = minAngle, minAngleRefPoints = minRefAngle,
                perturbedEnergy = perturbedEnergy, perturbationZ = zEnergyUnderPerturbation,
                perturbationPVal = pEnergyUnderPerturbation, perturbationSigma = finalPerturbationSigma))
    
    
}


## Get the minimum angle between each point and the closest neighboring point
getClosestPointAngles <- function(points) {

    numPoints = nrow(points)

    absDot = matrix(rep(0,numPoints*numPoints), ncol = numPoints, nrow = numPoints)
    
    ## Upper triangle
    for (i in 1:(numPoints-1)) {
        for (j in (i+1):numPoints) {
            absDot[i,j] = abs(dot(points[i,],points[j,]))
        }
    }

    minAngle = rep(0, numPoints)
    
    for (i in 1:numPoints) {
        minAngle[i] = acos( max(c(absDot[i,], absDot[,i])) ) * 180 / pi
    }

    return(minAngle)
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
