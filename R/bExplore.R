#' Exploratory backscatter layer analysis
#' 
#' Explore backscatter layer statistics and relationships prior to harmonization to inform the bulk shift procedure
#' 
#' @details 
#' `boxplot = TRUE` draws box plots comparing the full distributions of x and y, and also the distributions where they overlap. THis
#' is useful for assessing the representativeness of the area overlap (i.e., whether the error model will need to extrapolate).
#' 
#' Local polynomial fitting using loess can be slow with large sample sizes. Subsampling using `loess_samp` is used here only to fit 
#' the smoother; it is predicted on the entire dataset. Additional arguments can be passed to `...` (see [loess] for more information).
#' 
#' @param x,y SpatRasters. Backscatter layers to compare.
#' @param preds SpatRaster. Predictors to explore for explaining error between `x` and `y`.
#' @param error_map Logical. Whether to plot a map of error between `x` and `y`.
#' @param boxplot Logical. Whether to draw boxplots comparing `x` and `y`. See Details.
#' @param error_plot Logical. Whether to plot the error as a function of `x`, and also `preds` if specified.
#' @param loess Logical. Whether to use [loess] to visualized the relationship between `preds` (including `x`) and the error.
#' @param loess_samp Numeric. Maximum number of samples used to estimate the [loess] smoother. See Details.
#' @param ... Additional arguments to pass to [loess]. See Details.
#' 
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
