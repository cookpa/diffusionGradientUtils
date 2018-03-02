
##
## bvecsFullSphere = distributeShellsOverSphere(bvecs, bvals, shells, iterations = 250)
##
## Projects input bvecs to the upper hemisphere (z >= 0), then flips half of the
## vector so they point to the lower hemisphere (z < 0).
##
## Procedure for each shell:
##
## 1. Randomly choose half the points on the shell and put into lower hemisphere (vec <- -vec)
## 2. Compute electrostatic energy of points on upper hemisphere
## 3. Iterate over all points and swap any pair of points that reduces energy 
## 4. Repeat 2-4 until no swaps can be made
## 
## shells[1] is assumed to be unweighted and is not optimized
##
distributeShellsOverSphere <- function(bvals, bvecs, shells, iterations = 250) {
    
    shellIDs = bvalsToShells(bvals, shells)
    
    bvecsProj = bvecs
    
    numMeas = nrow(bvecs)
    
    ## Make unit vectors projected into hemisphere, bvecs may not have unit magnitude
    ## in multi-shell schemes
    uvecsProj = matrix(nrow = numMeas, ncol = 3)
    
    for (r in 1:numMeas) {
        
        if (shellIDs[r] > 1 ) {
            
            uvecsProj[r,] = as.numeric(bvecs[r,] / sqrt(sum(bvecs[r,]^2)))
            
            if (uvecsProj[r,3] < 0) {
                uvecsProj[r,] = -uvecsProj[r,]
            }
            
            if (bvecsProj[r,3] < 0) {
                bvecsProj[r,] = -bvecsProj[r,]
            }
            
        }
        else {
            uvecsProj[r,] = as.numeric(bvecs[r,])

            if (bvecsProj[r,3] < 0) {
                bvecsProj[r,] = -bvecsProj[r,]
            }

        }
        
    }

    ## Some of these will get negated as we do the optimization
    bvecsFinal = bvecsProj

    energy = electrostaticEnergy(uvecsProj,shellIDs)

    ## Energy of shell is the energy between all points in the shell

    ## Split each shell such that half the directions are negated, thus distributing over the full sphere
    
    
    for (sh in 2:length(shells)) {

        bestEnergy = Inf

        bestLower = c()

        print(paste("Optimizing shell ", sh))
        
        ## minimize energy of directions in the upper hemisphere (z > 0)
        for (it in 1:iterations) {
            
            ##      print(paste("Iteration ", it, " Shell", sh))
            
            allMeasOnShell = which(shellIDs == sh)
            
            ranOrder = sample(allMeasOnShell, length(allMeasOnShell))
            
            halfWay = ceiling( length(allMeasOnShell) / 2 )
            
            upperSubset = ranOrder[1:halfWay]
            lowerSubset = ranOrder[(halfWay+1):length(allMeasOnShell)]
            
            swapped = T
            
            upperSubsetEnergy = sum ( colSums(energy[upperSubset, upperSubset]) )
            
            ##    print(paste("  Initial energy ==" , upperSubsetEnergy))
            
            counter = 1
            
            while (swapped == TRUE) {
                
                swapped = FALSE
                
                ##  print(paste("    Iteration", counter))
                
                for (u in 1:length(upperSubset)) {
                    for (v in 1:length(lowerSubset)) {
                        tmpUpper = upperSubset
                        
                        ## Shall we swap this point into the upper subset?
                        candidate = lowerSubset[v]
                        
                        tmpLower = lowerSubset
                        
                        tmpLower[v] = tmpUpper[u]
                        
                        tmpUpper[u] = candidate
                        
                        tmpEnergy = sum ( colSums(energy[tmpUpper, tmpUpper]) )
                        
                        if (tmpEnergy < upperSubsetEnergy) {
                            upperSubset = tmpUpper
                            lowerSubset = tmpLower
                            upperSubsetEnergy = tmpEnergy
                            swapped = TRUE
                        }
                        
                    }
                }
                
                counter = counter + 1
                
            }
            
            ##     print(paste("  Final energy ==" , upperSubsetEnergy))
            
            if (upperSubsetEnergy < bestEnergy) {
                bestEnergy = upperSubsetEnergy
                bestLower = lowerSubset
            }
            
        }

        print(paste("Best energy for shell ", sh, " == " , bestEnergy))
        
        bvecsFinal[bestLower,] = -bvecsFinal[bestLower,] 
        
    }
    
    
    return(bvecsFinal)


}
