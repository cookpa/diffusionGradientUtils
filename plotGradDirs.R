library(rgl)
                                        ##
## plotGradDirs(bvecs, bvals, shells, projectToHemi = F, opacity = 0.5)
##
## Plots gradient directions on a unit sphere. Requires RGL
##
## Required args:
##
## bvals   -  N x 1 b-values
## bvecs   -  N x 3 gradient dirs
## shells  -  list of shells, eg c(0,300,800,2000). Individual b-values will be assigned to the closest shell
##
##
## Options:
##
## projectToUpper - if TRUE, project all data into the upper hemisphere (z >= 0)
##
## opacity - controls the opacity of the unit sphere, adjust as needed to assist depth perception.
##
## whichShells - list of shells to plot, eg c(2,3,4). Defaults to 2:length(shells). The first shell is assumed
## to be zero or close to zero, and is not plotted by default.
## 
## whichHemi - either "upper", "lower", or "both". Used to plot only directions in the upper or lower hemisphere
## (z > 0 or z < 0). Not compatible with projectToUpper = T
##
## 
##
plotGradDirs <- function(bvals, bvecs, shells, projectToUpper = F, opacity = 0.5, whichShells = NA, whichHemi = "both") {

  
    if (is.na(whichShells[1])) {
        whichShells = c(2:length(shells))
    }
    
    if ( nrow(bvecs) != length(bvals) ) {
        stop("Number of b-vectors and b-values must match")
    }
    
    if ( !(whichHemi %in% c("both", "upper", "lower") ) ) {
        stop("valid whichHemi options are both, upper, lower")
    }
    
    rgl::rgl.bg(color = "##FFFFFF")
    rgl::par3d(windowRect = c(100, 100, 600, 600))
    
    shellIDs = bvalsToShells(bvals, shells)
    
    shellColors = c("orange", "blue", "green", "red", "yellow", "cyan", "magenta")
    
    ## points on the unit sphere
    pointRadius = 0.03

    spheres3d(c(0,0,0), radius = 1, color = "black", front = "line", back = "line", lit = F)
    spheres3d(c(0,0,0), radius = 0.999, color = "gray", alpha = opacity, back = "cull")
    
    rgl.pop("lights")
    
    light3d(ambient = "##FFFFFF", specular = "##999999", diffuse = "##DDDDDD")
    
    
    for (r in 1:nrow(bvecs)) {
        
        if (shellIDs[r] %in% whichShells ) {

            uVec = bvecs[r,] / sqrt(sum(bvecs[r,]^2))
            
            if (projectToUpper) {
                if (uVec[3] < 0) {
                    uVec = -uVec
                }
            }
            
            plotPoint = T      
            
            if (whichHemi == "upper") {
                if (uVec[3] < 0) {
                    plotPoint = F
                }
            }
            else if (whichHemi == "lower") {
                if (uVec[3] >= 0) {
                    plotPoint = F
                }
                
            }
            
            if (plotPoint) {
                spheres3d(uVec, radius = pointRadius, color = shellColors[shellIDs[r]], back = "cull")
            }
            
        }
        
    }


}

