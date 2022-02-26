#' @param r SpatRaster

raScale <- function(r, min = 0.1, max = 0.9){
  rmin <- as.numeric(global(r, fun='min', na.rm = TRUE))
  rmax <- as.numeric(global(r, fun='max', na.rm = TRUE))
  
  r <- min + ((r - rmin)*(max-min)) / (rmax - rmin)
}

#' Advanced subsampling of a SpatRaster for spatial modelling
#' 
#' @param x SpatRaster layer
#' @param size sample size
#' @param method vector of subsampling methods
#' 
#' #' @import terra
#' #' @export

bSample <- function(x, size, samplemethods = c('uniform'), ...){
  
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
    xstrat <- lapply(xstrat, FUN = function(r) autocor(r, global = FALSE, ...))
    xstrat <- lapply(xstrat, FUN = function(r) 1 - raScale(r))
    
    s <- lapply(xstrat, FUN = function(r) spatSample(r, size = size, method = 'weights', values = FALSE, cells = TRUE))
    s <- do.call(rbind, s)
    
    s <- s$cell
  } else {
    s <- lapply(xstrat, FUN = function(r) spatSample(r, size = size, na.rm = TRUE, values = FALSE, cells = TRUE))
    s <- do.call(rbind, s)
  }
  
  return(as.vector(s))
}
