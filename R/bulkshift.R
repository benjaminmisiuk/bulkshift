#' Helper functions to calculate model evaluation statistics
#'
#' @param y vector of observed values
#' @param y_hat vector of predicted values

ve <- function(y, y_hat){
  1 - (sum((y - y_hat)^2) / sum((y - mean(y))^2))
}

mae <- function(y, y_hat){
  mean(abs(y - y_hat))
}

#' @import terra
#' @importFrom raster raster
#' @export

bulkshift <- function(shift, target, preds = NULL, model = "glm", mosaic = FALSE, mosaicmethod = "bilinear", savedata = TRUE, sample = NULL, samplemethods = "uniform", crossvalidate = NULL, ...){
  
  #separate the package into exploratory and modelling
  #strat still not working
  #two step model?
  #gam
  #penalized regression
  #move stratification to helper function?
  #weight stratification by autocorrelation?
  #spatial blocking?
  #pre and post stats? ks?
  #test using 1 RasterLayer, 1 SpatRaster, 1 RasterStack, all RasterLayer, 1 RasterStack other RasterLayer
  #spatial subsample
  
  #check for supported models
  if(!model %in% c('mean', 'glm', 'randomForest')) stop('model must be one of "mean", "glm", or "randomForest"')
  
  #note initial classes then convert any RasterLayers to SpatRasters
  classes <- c(class(shift)[1], class(target)[1], class(preds)[1]); classes <- classes[classes != 'NULL']
  isRasterLayer <- (sum(classes == 'RasterLayer' | classes == 'RasterStack') / length(classes)) == 1
  
  if(class(shift)[1] == 'RasterLayer') shift <- rast(shift)
  if(class(target)[1] == 'RasterLayer') target <- rast(target)
  if(class(preds)[1] == 'RasterLayer' | class(preds)[1] == "RasterStack") preds <- rast(preds)
  
  if(class(shift)[1] != 'SpatRaster' | class(target)[1] != 'SpatRaster') stop('all inputs must be rasters.')
  if(class(preds)[1] != 'SpatRaster' & !is.null(preds)) stop('predictors must be rasters.')

  #check for crs matches
  if(crs(shift) != crs(target)) stop('crs of shift and target layers do not match.')

  if(!is.null(preds)){
    if(crs(preds) != crs(shift)) stop('crs of preds does not match the shift layer.')
  }
  
  #resample the target if necessary
  if(ext(shift) != ext(target) | res(shift)[1] != res(target)[1] | res(shift)[2] != res(target)[2]){
    target_re <- resample(target, shift)
  } else {
    target_re <- target
  }
  
  #calculate the error between surveys
  err_ras <- (target_re - shift); names(err_ras) <- 'error'
  
  if(is.na(minmax(err_ras)[1]) & is.na(minmax(err_ras)[2])) stop('no overlap detected.')
    
  #create a SpatRaster of predictors
  if(!is.null(preds)){
    if(ext(preds) != ext(shift) | res(preds)[1] != res(shift)[1] | res(preds)[2] != res(shift)[2]){
      preds <- resample(preds, shift)
    }
    preds <- rast(
      list(
        shift,
        preds
      )
    )
  } else {
    preds <- shift
  }
  
  #mask SpatRasters at the area of overlap 
  ovlp <- rast(
    list(
      err_ras,
      mask(preds, err_ras)
    )
  )
  
  #extract values at the area of overlap, using subsampling if specified
  if(!is.null(sample)){
    sample <- floor(length(cells(ovlp$error)) * sample)
    s <- bSample(ovlp[[2]], size = sample, samplemethods = samplemethods)
    df <- ovlp[s]
  } else {
    df <- as.data.frame(ovlp)
  }
  
  df <- df[complete.cases(df), ]
  rm(ovlp, err_ras)
  
  #set aside validation data if "crossvalidate" is specified
  if(!is.null(crossvalidate)){
    cv_i <- sample(c(TRUE, FALSE), nrow(df), prob = c(crossvalidate, 1-crossvalidate), replace = TRUE)
    df_out <- df[cv_i, ]
    df <- df[!cv_i, ]
  }
  
  #model the error
  form <- as.formula(
    paste0(
      names(df)[1], 
      '~',
      paste0(names(df)[-1], collapse="+")
    )
  )
  
  if(model == 'mean') err_mod <- glm(df[ ,1] ~ 1)
  if(model == 'glm') err_mod <- glm(form, data = df, ...)
  if(model == 'randomForest') err_mod <- randomForest::randomForest(form, data = df, ...)
  
  #use the model to predict the error over the 'shift' dataset
  err_pred <- predict(preds, err_mod, type = 'response')
  shifted <- shift + err_pred; names(shifted) <- paste0(names(shift), '_shifted')
  
  #calculate fitted model statistics
  p <- predict(err_mod, df, type = 'response')
  fitStats <- data.frame(
    VE = ve(y = df$error + df$bb2017, y_hat = p + df$bb2017),
    MAE = mae(y = df$error + df$bb2017, y_hat = p + df$bb2017),
    r = cor(df$error + df$bb2017, p + df$bb2017)
  )
  
  #prepare output
  #if all of the original inputs were from the "raster" package, export class "RasterLayer"; otherwise, SpatRaster
  if(isRasterLayer){
    out <- list(
      shifted = raster(shifted),
      errorModel = err_mod,
      fitStats = fitStats
    )
  } else {
    out <- list(
      shifted = shifted,
      errorModel = err_mod,
      fitStats = fitStats
    )
  }
  
  #validate using witheld data if specified by "crossvalidate"
  if(!is.null(crossvalidate)){
    p <- predict(err_mod, df_out, type = 'response')
    
    testStats <- data.frame(
      VE = ve(y = df_out$error + df_out$bb2017, y_hat = p + df_out$bb2017),
      MAE = mae(y = df_out$error + df_out$bb2017, y_hat = p + df_out$bb2017),
      r = cor(df_out$error + df_out$bb2017, p + df_out$bb2017)
    )
    out <- c(out, list(testStats = testStats))
  }
  
  if(savedata) out <- c(out, list(data=df))

  if(mosaic){
    if(ext(shift) != ext(target) | res(shift)[1] != res(target)[1] | res(shift)[2] != res(target)[2]){ #if the original "target" and "shift" did not align, choose to resample the "shifted" with respect to the original "target"
      target <- extend(target, shifted)
      bs_mosaic <- terra::mosaic(
        resample(shifted, target, method = mosaicmethod),
        target
      )
      names(bs_mosaic) <- paste0(names(target), '_', names(shift), '_m')
    }
    #if all of the original inputs were from the "raster" package, export mosaic as class "RasterLayer"
    if(isRasterLayer){
      out <- c(out, mosaic =  raster(bs_mosaic))
    } else {
      out <- c(out, mosaic =  bs_mosaic)
    }
  }
  return(out)
}