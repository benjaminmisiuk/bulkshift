#function for rescaling SpatRaster layers

raScale <- function(r, min = 0.1, max = 0.9){
  rmin <- as.numeric(global(r, fun='min', na.rm = TRUE))
  rmax <- as.numeric(global(r, fun='max', na.rm = TRUE))
  
  r <- min + ((r - rmin)*(max-min)) / (rmax - rmin)
}

#function to draw random non-NA cell numbers from a single SpatRaster layer

rSample <- function(x, n, prob = FALSE){
  
  if(prob){
    p <- as.vector(x)
    p <- p[!is.na(p)]
    s <- which(!is.na(values(x)))
    
    out <- sample(s, size = n, prob = p)
  } else {
    s <- which(!is.na(values(x)))
    out <- sample(s, size = n)
  }
  return(out)
}

#' Subsampling methods for a SpatRaster
#' 
#' Advanced subsampling of a [terra] SpatRaster for spatial modelling
#' 
#' @details If `samplemethods = "stratify"`, samples are distributed equally among quartiles of x. 
#' `samplemethods = "autocorrelation"` weights sampling by the inverse of local Moran's I, which is
#' calculated used [terra::autocor()].
#' 
#' @param x SpatRaster layer
#' @param size Sample size.
#' @param samplemethods Vector of subsampling methods. One or multiple of "uniform", "stratify", "autocorrelation".
#' 
#' @return Vector of cell numbers
#' 
#' @examples 
#' bb2017 <- rast(system.file('extdata', 'bb2017.tif', package='bulkshift'))
#' 
#' s <- bSample(x = bb2017, size = 1000, samplemethods = c("stratify"))
#' plot(bb2017)
#' plot(as.points(bb2017, na.rm = FALSE)[s], add = TRUE)
#' 
#' @import terra
#' @export

bSample <- function(x, size, samplemethods = c('uniform')){
  
  if(!all(samplemethods %in% c("uniform", "stratify", "autocorrelation"))) stop('all "samplemethods" must be one of "uniform", "stratify" or "autocorrelation"')
  
  xstrat <- x

  if('stratify' %in% samplemethods){
    q <- global(x, fun=quantile, na.rm = TRUE)
    
    mat <- matrix(
      c(
        q[[1]], q[[2]], 1,
        q[[2]], q[[3]], 2,
        q[[3]], q[[4]], 3,
        q[[4]], q[[5]], 4
      ),
      ncol = 3,
      byrow = TRUE
    )
    
    strat <- terra::classify(x, rcl = mat, include.lowest = TRUE); names(strat) <- 'strata'
    xstrat <- rast(
      list(
        q = mask(x, strat, inverse = TRUE, maskvalues = 1),
        q = mask(x, strat, inverse = TRUE, maskvalues = 2),
        q = mask(x, strat, inverse = TRUE, maskvalues = 3),
        q = mask(x, strat, inverse = TRUE, maskvalues = 4)
      )
    )
    
    size <- round(size / nlyr(xstrat))
  }
  
  if('autocorrelation' %in% samplemethods){
    xstrat <- lapply(xstrat, FUN = function(r) autocor(r, global = FALSE))
    xstrat <- lapply(xstrat, FUN = function(r) 1 - raScale(r))
    
    s <- lapply(xstrat, FUN = function(r) rSample(r, n = size, prob = TRUE))
    s <- do.call(c, s)
  } else {
    s <- lapply(xstrat, FUN = function(r) rSample(r, n = size))
    s <- do.call(c, s)
  }
  
  return(as.vector(s))
}
