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

bExplore <- function(x, y, preds = NULL, error_map = TRUE, boxplot = TRUE, error_plot = TRUE, loess = TRUE, loess_samp = 10000, ...){
  
  #extract the overlap between backscatter layers
  ovlp <- overlap(x, y)
  err <- ovlp[[1]] - ovlp[[2]]; names(err) <- 'error'
  
  if(error_map){
    readline(prompt="Press [enter] for next plot")
    plot(err, main = "Error map", reset = TRUE)
  }
  
  if(!is.null(preds)){
    ovlp_p <- rast(
      list(
        err,
        overlap(ovlp[[1]], preds)
      )
    )
  } else {
    ovlp_p <- rast(
      list(
        err,
        ovlp[[1]]
      )
    )
  }
  
  df <- as.data.frame(ovlp)
  df_p <- as.data.frame(ovlp_p)
  
  #compare full data distributions to distributions at the overlap 
  if(boxplot){
    readline(prompt="Press [enter] for next plot")
    boxplot(
      list(as.vector(x), as.vector(y), df[ ,1], df[ ,2]),
      names = c(names(x), names(y), paste(names(x), 'overlap'), paste(names(y), 'overlap'))
    )
  }
  
  #observe the error as a function of x
  if(error_plot){
    for(i in 2:ncol(df_p)){
      readline(prompt="Press [enter] for next plot")
      if(loess){
        s = sample(nrow(df_p), min(nrow(df_p), loess_samp))
        plot(df_p[ ,i], df_p$error, xlab = colnames(df_p)[i], ylab = "Error")
        lo <- loess(df_p[s,'error'] ~ df_p[s,i], ...)
        lines(
          sort(df_p[ ,i]), 
          predict(lo, sort(df_p[ ,i])),
          col = 'red'
        )
      } else {
        plot(df_p[ ,i], df_p$error, xlab = colnames(df_p)[i], ylab = "Error")
      }
    }
  }
  
}
