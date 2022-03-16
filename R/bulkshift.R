#Functions to calculate model evaluation statistics
ve <- function(y, y_hat){
  1 - (sum((y - y_hat)^2) / sum((y - mean(y))^2))
}

mae <- function(y, y_hat){
  mean(abs(y - y_hat))
}

#' Sonar backscatter bulk shift
#' 
#' Relative calibration of two partially overlapping sonar backscatter datasets using the methods presented in Misiuk et al. (2020). This is largely a wrapper for the [terra] package to facilitate relative calibration of backscatter rasters. Supports several modelling methods and optionally, spatially explicit subsampling and validation.
#' 
#' @param shift SpatRaster. Backscatter dataset undergoing correction.
#' @param target SpatRaster. Backscatter dataset used as reference.
#' @param preds SpatRaster. One or more layers to use as additional predictor variables for backscatter calibration model.
#' @param model character. Method used to model error between backscatter datasets. Currently supported models are "mean" (i.e., intercept-only), "glm", and "randomForest". Model defaults are retained; use "..." to specify additional parameters.
#' @param mosaic logical. Mosaic and output the backscatter layers?
#' @param mosaicmethod character. Which method used to resample the corrected backscatter layer for mosaicking? See [terra::resample()] for details.
#' @param savedata logical. Whether to output the model data.frame.
#' @param sample numeric. Proportion of overlapping data to sample for modelling using [samplemethods()].
#' @param samplemethods character vector. The method used to subsample the data. One or several of "uniform" (the default), "stratify", or "autocorrelation". Specifying multiple methods results in a combined output. See Details.
#' @param crossvalidate numeric. Proportion of data to use for validation. Validation data are drawn from the dataset following saubsampling if [sample] is used.
#' @param ... Additional parameters to pass to models or subsampling methods.
#' 
#' @return List of bulkshift objects
#' 
#' @references 
#' Misiuk, B., Brown, C.J., Robert, K., Lacharite, M., 2020. Harmonizing Multi-Source Sonar Backscatter Datasets for Seabed Mapping Using Bulk Shift Approaches. Remote Sensing 12, 601. https://doi.org/10.3390/rs12040601
#' 
#' @import terra
#' @importFrom raster raster
#' @export

bulkshift <- function(shift, target, preds = NULL, model = "glm", mosaic = FALSE, mosaicmethod = "bilinear", savedata = TRUE, sample = NULL, samplemethods = "uniform", crossvalidate = NULL, ...){
  
  #separate the package into exploratory and modelling
  #the sample index?
  #gam
  #penalized regression
  #spatial blocking?
  #test using 1 RasterLayer, 1 SpatRaster, 1 RasterStack, all RasterLayer, 1 RasterStack other RasterLayer
  #'details'
  #bSample docs
  #examples

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
  
  #calculate error and get all SpatRasters masked at the overlap
  ovlp <- overlap(shift, target); names(ovlp) <- c("shift", "target")
  err_ras <- ovlp$target - ovlp$shift; names(err_ras) <- 'error'
  
  if(is.na(minmax(err_ras)[1]) & is.na(minmax(err_ras)[2])) stop('no overlap detected.')
  
  if(!is.null(preds)){
    ovlp_p <- overlap(shift, preds); names(ovlp_p[[1]]) <- "shift"
  } else {
    ovlp_p <- shift; names(ovlp_p) <- "shift"
  }
  
  ovlp <- rast(
    list(
      err_ras,
      mask(ovlp_p, err_ras)
    )
  )
    
  #extract values at the area of overlap, using subsampling if specified
  if(!is.null(sample)){
    sample <- floor(length(cells(ovlp$error)) * sample)
    s <- bSample(ovlp$shift, size = sample, samplemethods = samplemethods)
    df <- ovlp[s]
  } else {
    df <- as.data.frame(ovlp)
  }
  
  df <- df[complete.cases(df), ]
  rm(ovlp)
  
  #set aside validation data if "crossvalidate" is specified
  if(!is.null(crossvalidate)){
    cv_i <- rep(FALSE, nrow(df)); cv_i[1:round(crossvalidate*nrow(df))] <- TRUE
    cv_i <- sample(cv_i)
    df_out <- df[cv_i, ]
    df <- df[!cv_i, ]
  }
  
  #model the error
  form <- as.formula(
    paste0(
      'error', 
      '~',
      paste0(names(df)[-1], collapse="+")
    )
  )
  
  if(model == 'mean') err_mod <- glm(df[ ,1] ~ 1)
  if(model == 'glm') err_mod <- glm(form, data = df, ...)
  if(model == 'randomForest') err_mod <- randomForest::randomForest(form, data = df, ...)
  
  #use the model to predict the error over the 'shift' dataset
  err_pred <- predict(ovlp_p, err_mod, type = 'response')
  shifted <- shift + err_pred; names(shifted) <- paste0(names(shift), '_shifted')
  
  #calculate fitted model statistics before and after correction
  fitStats <- data.frame(
    layer = "Original",
    VE = ve(y = df$error + df$shift, y_hat = df$shift),
    MAE = mae(y = df$error + df$shift, y_hat = df$shift),
    r = cor(df$error + df$shift, df$shift)
  )
  p <- predict(err_mod, df, type = 'response')
  fitStats <- rbind(
    fitStats,
    data.frame(
      layer = "Corrected",
      VE = ve(y = df$error + df$shift, y_hat = p + df$shift),
      MAE = mae(y = df$error + df$shift, y_hat = p + df$shift),
      r = cor(df$error + df$shift, p + df$shift)
    )
  )
  
  #prepare output
  
  #if all original inputs were from the "raster" package, export class "RasterLayer"; otherwise, SpatRaster
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
    testStats <- data.frame(
      layer = "Original",
      VE = ve(y = df_out$error + df_out$shift, y_hat = df_out$shift),
      MAE = mae(y = df_out$error + df_out$shift, y_hat = df_out$shift),
      r = cor(df_out$error + df_out$shift, df_out$shift)
    )
    p_out <- predict(err_mod, df_out, type = 'response')
    testStats <- rbind(
      testStats,
      data.frame(
        layer = "Corrected",
        VE = ve(y = df_out$error + df_out$shift, y_hat = p_out + df_out$shift),
        MAE = mae(y = df_out$error + df_out$shift, y_hat = p_out + df_out$shift),
        r = cor(df_out$error + df_out$shift, p_out + df_out$shift)
      )
    )
    out <- c(out, list(testStats = testStats))
  }
  
  if(savedata){
    out <- c(
      out, 
      list(data = data.frame(target = df$shift + df$error, df, shifted = df$shift + p))
    )
    
    if(!is.null(crossvalidate)){
      out <- c(
        out, 
        list(
          dataVal = data.frame(target = df_out$shift + df_out$error, df_out, shifted = df_out$shift + p_out)
        )
      )
    }
  }

  if(mosaic){
    #resample the "shift" with respect to the "target" for mosaicking 
    bs_mosaic <- terra::mosaic(
      project(shifted, target, method = mosaicmethod, align = TRUE),
      target
    )
    names(bs_mosaic) <- paste0(names(target), '_', names(shift), '_m')
    
    #if all of the original inputs were from the "raster" package, export mosaic as class "RasterLayer"
    if(isRasterLayer){
      out <- c(out, mosaic =  raster(bs_mosaic))
    } else {
      out <- c(out, mosaic =  bs_mosaic)
    }
  }
  return(out)
}