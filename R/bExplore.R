#' Exploratory backscatter layer analysis
#' 
#' Explore backscatter layer statistics and relationships prior to harmonization to inform the bulk shift procedure
#' 
#' @details 
#' 
#' @param x SpatRaster layer
#' 
#' @return
#' 
#' @examples
#' 
#' @import terra
#' @export

bExplore <- function(x, y, boxplot = TRUE, loess = TRUE, ...){
  
  #extract the overlap between backscatter layers
  ovlp <- overlap(x, y)
  df <- as.data.frame(ovlp)
  
  #compare full data distributions to distributions at the overlap 
  if(boxplot){
    boxplot(
      list(as.vector(x), as.vector(y), df[ ,1], df[ ,2]),
      names = c(names(x), names(y), paste(names(x), 'overlap'), paste(names(y), 'overlap'))
    )
  }
  
  #observe the error as a function of x
  if(loess){
    s = sample(nrow(df), 10000)
    plot(df[ ,1], df[ ,2] - df[ ,1])
    lo <- loess((df[s,2] - df[s,1]) ~ df[s,1], ...)
    lines(
      sort(df[ ,1]), 
      predict(lo, sort(df[ ,1])),
      col = 'red'
    )
    
  } else {
    plot(x, y - x)
  }
}