#' Backscatter error visualization
#' 
#' Visualize the modelled error between backscatter datasets, and optionally, additional predictors, in 2 and 3D.
#' 
#' @details 
#' This function generates 2D or 3D scatterplots from the output of [bulkshift] using the [plot3D] and [plot3Drgl] packages to assist in diagnosing the quality of model fit. 
#' Sensible plotting parameters are implemented by default but see those packages for additional plotting options.
#' 
#' The output shows the error between backscatter datasets on the y-axis as a function of the "shift" dataset on the x-axis. 
#' If no other predictors were used for [bulkshift], this is a 2D scatterplot. 
#' If additional predictors were supplied, they are plotted on the z-axis and a 3D scatterplot is generated.
#' 
#' In two dimensions, the backscatter error model is shown as a red line, which should pass through the center of the points.
#' In three dimensions, the error model is shown as a plane.
#' 
#' 
#' @param x An (unmodified) list output from [bulkshift].
#' @param interactive Logical. If additional predictors were supplited to [bulkshift], should an interactive 3D plot be returned?
#' @param cex Size of scatterplot points.
#' @param pch Scatterplot point symbols.
#' @param alpha Alpha (transparency) of points in scatterplot.
#' @param theta Azimuthal perspective on the 3D plot if `interactive = FALSE`.
#' @param phi Colatitude perspective on the 3D plot if `interactive = FALSE`.
#' 
#' @examples
#' bb2016 <- rast(system.file('extdata', 'bb2016.tif', package='bulkshift'))
#' bb2017 <- rast(system.file('extdata', 'bb2017.tif', package='bulkshift'))
#' bbdepth <- rast(system.file('extdata', 'bbdepth.tif', package='bulkshift'))
#' 
#' #2D example
#' b <- bulkshift(bb2017, bb2016)
#' bError(b)
#' 
#' #3D example
#' b <- bulkshift(bb2017, bb2016, bbdepth)
#' bError(b)
#' 
#' @export
#' 
bError <- function(x, interactive = TRUE, cex = 0.75, pch = 16, alpha = 1, theta = 45, phi = 0){
  
  if(!exists('data', where = b)) stop("No dataframe found. Re-run bulkshift() with savedata = TRUE.")
  
  
  #if there are amy additional predictors in the dataframe plot a 3D model
  if(any(!(names(b$data) %in% c('target', 'error', 'shift', 'shifted')))){
    pred <- names(b$data)[!(names(b$data) %in% c('target', 'error', 'shift', 'shifted'))]
    
    if(length(pred) > 1) stop("Too many predictors to plot in three dimensions")
    
    grid.lines <- 100
    x.pred <- seq(min(b$data$shift), max(b$data$shift), length.out = grid.lines)
    y.pred <- seq(min(b$data[ ,pred]), max(b$data[ ,pred]), length.out = grid.lines)
    xy <- expand.grid(x.pred, y.pred)
    names(xy) <- c('shift', pred)
    
    z.pred <- matrix(predict(b$errorModel, xy), 
                     nrow = grid.lines, ncol = grid.lines)
    
    plot3D::scatter3D(
      x = b$data$shift, 
      y = b$data[ ,pred], 
      z = b$data$error, 
      bty = "g", 
      surf = list(x = x.pred, y = y.pred, z = z.pred, col = "gray40", facets = NA, lwd = 0.25),
      xlab = "Backscatter (dB)",
      ylab = pred, 
      zlab = "Error (dB)", 
      theta = theta, 
      phi = 0, 
      clab = "Error", 
      cex = cex, 
      pch = pch,
      col = plot3D::ramp.col(c('red2', 'lightgoldenrod', 'green4'), alpha = alpha),
      ticktype = "detailed", 
      plot = !interactive
    )
    
    if(interactive) plot3Drgl::plotrgl()
  } else {
    plot3D::scatter2D(
      x = b$data$shift,
      y = b$data$error,
      xlab = "Shift backscatter (dB)",
      ylab = "Error (dB)",
      col = 'black',
      cex = cex,
      pch = pch,
      alpha = alpha
    )
    plot3D::lines2D(
       x = b$data$shift[order(b$data$shift)],
       y = predict(b$errorModel, b$data[order(b$data$shift), ]),
       add = TRUE,
       col = 'red'
     )
  }
}
