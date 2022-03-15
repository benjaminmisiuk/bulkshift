#' Overlap
#' 
#' Extract the union of overlapping SpatRasters with different extents, origins, and projections.
#' 
#' @param x SpatRaster. Backscatter dataset undergoing correction.
#' @param y SpatRaster. Will be algined and/or projected to x if extent, origin, or crs do not match.
#' 
#' @return SpatRaster
#' 
#' @import terra
#' @export

overlap <- function(x, y){
  
  #check for crs matches
  if(crs(y) != crs(x)){
    warning('projecting data to match CRS.')
    target <- project(y, x)
  }
  
  #resample SpatRaster y if necessary
  if(ext(y) != ext(x) | res(y)[1] != res(x)[1] | res(y)[2] != res(x)[2]){
    y <- resample(y, x)
  }
  
  #mask x and y at the area of overlap 
  ovlp <- rast(
    list(
      mask(x, y),
      mask(y, x)
    )
  )
  
  return(ovlp)
}
