#' Extract SpatRaster overlap
#' 
#' Extract the union of overlapping [terra] SpatRasters with different extents, origins, and projections.
#' 
#' @details SpatRaster y will be transformed (e.g., projected, resampled) with respect to x if 
#' extents, origins, or crs do not match.
#' 
#' @param x,y SpatRaster.
#' 
#' @return SpatRaster
#' 
#' @examples 
#' bb2016 <- rast(system.file('extdata', 'bb2016.tif', package='bulkshift'))
#' bb2017 <- rast(system.file('extdata', 'bb2017.tif', package='bulkshift'))
#' 
#' ovlp <- overlap(bb2016, bb2017)
#' plot(ovlp)
#' 
#' #this function could be used to facilitate your own bespoke bulkshift, or to use a different model
#' err <- c(ovlp$bb2016 - ovlp$bb2017)
#' names(err) <- "err"
#' model_rast <- c(err, ovlp$bb2017)
#' df <- as.data.frame(model_rast)
#' 
#' model <- glm(err ~ bb2017, data = df)
#' 
#' bb2017_shifted <- bb2017 + terra::predict(bb2017, model)
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
