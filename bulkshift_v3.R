bulkshift <-
  function(shift,                   #dataset being corrected
           target,                  #reference dataset
           preds=NULL,              #list of optional explanatory variables
           shift.method="lm",       #the model used for correction: "mean", "lm", "brt"
           mosaic=FALSE,            #mosaic the corrected "shift" layer and the target?
           mosaic.method="bilinear",#how should the shift layerr be resampled for mosaicking? "bilinear" (bilinear interpolation), or "ngb"
           save.data=FALSE,         #output the bulk shift data?
           err.plots=FALSE,         #output 2d (from bivariate models) or 3d (from multivariate models) error plots?
           dist.plots=FALSE,        #output distribution (cdf and pdf) plots for shifted data? 
           save.dist.plots=FALSE,   #if TRUE, save the plots in the list output rather than plotting to the graphics device
           save.err.plots=FALSE,    #if TRUE, save plots in the list output rather than plotting to the graphics device
           sample=FALSE,            #use a subsample of the dataset? may be necessary for BRTs or plotting with large datasets
           samp.size=0.25,          #subsample size as a proportion of all overlapping raster cells
           
           #arguments passed to glm() when shift.method="lm"
           family=gaussian,  #error distribution
           link="identity",  #link function
           
           #arguments passed to gbm.step() when shift.method="brt"
           tc=2,            #tree complexity
           lr=0.01,         #learning rate
           bf=0.5,          #bag fraction
           verb=FALSE       #verbose BRT?
  ){
    
    #checking for required packages; installing if not present
    required.packages <- c("raster", "rgdal", "dismo", "ggplot2", "plot3Drgl")
    for(i in 1:length(required.packages)){
      if(!(required.packages[i] %in% installed.packages())) install.packages(required.packages[i])
    }
    
    require(raster)
    require(rgdal)
    require(dismo)
    require(ggplot2)
    
    #aligning extents and origins of target and shift rasters; extracting the area of overlap
    cat("Aligning rasters... \n")
    names(shift) <- "shift"
    names(target) <- "target"
    
    if(xres(shift) != xres(target)) warning("Rasters have different resolutions. Resampling the shift layer...", call.=FALSE, immediate.=TRUE)
    
    shift.mask <- mask(resample(shift, target), target)
    target.mask <- mask(target, shift.mask)
    
    #for corrections using the shift layer and additional covariates (listed in the "preds" argument)
    if(!is.null(preds)) {
      
      back_data <- data.frame(as.data.frame(target.mask), as.data.frame(shift.mask))
      names(back_data) <- c("target", "shift")
      preds <- as.list(preds)
      
      #extract the covariate values for area of overlap
      for (i in 1:length(preds)){
        pred.mask <- mask(resample(preds[[i]], shift.mask), shift.mask)
        back_data <- data.frame(back_data, as.data.frame(pred.mask))
      }
      
      #remove NAs from the working data.frame and calculate the error between layers
      back_data <- back_data[complete.cases(back_data),]
      back_data <- data.frame(back_data, err=back_data$target-back_data$shift)
      
      #for subsampling large datasets
      if(sample==TRUE){
        if((samp.size*nrow(back_data))>nrow(back_data)){
          stop("Subsample size is greater than number of raster cells; use a proportion < 1.")
        } else back_data <- back_data[seq(1, nrow(back_data), length.out=samp.size*nrow(back_data)), ]
      }
      
      #for mean corrections
      if(shift.method == "mean") {
        
        #shift data in working data.frame
        cat("Shifting layer by the mean of the error... \n")
        error.model <- mean(back_data$err)
        back_data <- data.frame(back_data, shifted=mean(back_data$err)+back_data$shift)
        
        #apply the same correction to the original shift raster
        cat("Predicting new raster... \n")
        shifted <- shift + mean(back_data$err)
        
        #for error plotting
        if(err.plots == "TRUE"){
          
          cat("Plotting error... \n")
          
          if(save.err.plots){
            
            p1 <- ggplot(back_data, aes(x=shift, y=err)) + 
              theme_minimal(11) +
              geom_point(col="black", size=0.5, alpha=0.3, shape=16) +
              geom_hline(aes(yintercept=error.model, lty="Error model"), col="red", alpha=1, lwd=1) +
              scale_linetype(name=NULL) +
              xlab("Shift backscatter (dB)") + 
              ylab("Error (dB)") +
              theme(axis.text = element_text(color = "black", family = "sans"), legend.position=c(0.89,0.885), panel.border=element_rect(fill=NA, colour="gray80", size=1.1), legend.text = element_text(color = "black", family = "sans", size=12), legend.background = element_blank())
            
            plot(p1)
            
          } else{
            
            plot(ggplot(back_data, aes(x=shift, y=err)) + 
                   theme_minimal(11) +
                   geom_point(col="black", size=0.5, alpha=0.3, shape=16) +
                   geom_hline(aes(yintercept=error.model, lty="Error model"), col="red", alpha=1, lwd=1) +
                   scale_linetype(name=NULL) +
                   xlab("Shift backscatter (dB)") + 
                   ylab("Error (dB)") +
                   theme(axis.text = element_text(color = "black", family = "sans"), legend.position=c(0.89,0.885), panel.border=element_rect(fill=NA, colour="gray80", size=1.1), legend.text = element_text(color = "black", family = "sans", size=12), legend.background = element_blank())
            ) 
          }
        }
      }
      
      #for lm corrections
      if(shift.method == "lm"){
        
        #create data.frame of predictors and fit the model
        cat("Fitting multiple regression... \n")
        mod.vars <- back_data[2:length(back_data)]
        error.model <- glm(err ~ ., family=family(link=link), data=mod.vars)
        
        #align and stack the covariate rasters to the extent of the shift layer for predicting 
        cat("Predicting new raster... \n")
        pred.stack <- stack()
        for (i in 1:length(preds)){
          pred.mask <- mask(resample(preds[[i]], shift), shift)
          pred.stack <- stack(pred.stack, pred.mask)
        }
        
        pred.stack <- stack(mask(shift, pred.stack[[1]]), pred.stack)
        names(pred.stack) <- names(mod.vars[1:(length(mod.vars)-1)])
        
        #predict the shifted raster using the covariates in the stack; predict corrections to the working data.frame
        shifted <- (predict(pred.stack, error.model)) + shift
        back_data <- data.frame(back_data, shifted=(predict(error.model, back_data, type = "response")) + back_data$shift)
        
        #for 3d error plots when there are two predictors (e.g., the shift layer and one covariate)
        if(err.plots == "TRUE") {
          if(length(preds)==1){
            
            cat("Plotting error... \n")
            
            #checking for required packages; installing if not present
            required.packages <- c("plot3D", "plot3Drgl")
            for(i in 1:length(required.packages)){
              if(!(required.packages[i] %in% installed.packages())) install.packages(required.packages[i])
            }
            
            require(plot3D)
            require(plot3Drgl)
            
            #create grids of predictors for plotting the model in 3d
            grid.lines <- 100
            x.pred <- seq(min(mod.vars[1]), max(mod.vars[1]), length.out = grid.lines)
            y.pred <- seq(min(mod.vars[2]), max(mod.vars[2]), length.out = grid.lines)
            xy <- expand.grid(x.pred, y.pred)
            names(xy) <- names(c(mod.vars[1],mod.vars[2]))
            
            #predict z-values of the and y grids using the model
            if(shift.method == "brt"){
              z.pred <- matrix(predict(error.model, xy, n.trees = error.model$n.trees, type = "response"), 
                               nrow = grid.lines, ncol = grid.lines)
            } else {
              z.pred <- matrix(predict(error.model, xy), 
                               nrow = grid.lines, ncol = grid.lines)
            }
            
            if(save.err.plots){
              
              plot.new()
              #create a 3d plot of the data with predictors on x- and y-axes, and error on z; plot the model as a surface
              scatter3D(mod.vars[[1]], mod.vars[[2]], mod.vars[[3]], 
                        bty="g", surf=list(x = x.pred, y = y.pred, z = z.pred, col="gray40", facets=NA, lwd=0.25),
                        xlab="Backscatter (dB)", ylab="Var 1", zlab="Error (dB)", 
                        theta=45, phi=0, 
                        clab="Error", 
                        cex=0.5, pch=16,
                        col=ramp.col(c('red2', 'lightgoldenrod', 'green4'), alpha=0.25),
                        ticktype = "detailed", 
                        plot=TRUE)
              
              plotrgl()
              
              p1 <- recordPlot()
            } else {
              
              #create a 3d plot of the data with predictors on x- and y-axes, and error on z; plot the model as a surface
              scatter3D(mod.vars[[1]], mod.vars[[2]], mod.vars[[3]], 
                        bty="g", surf=list(x = x.pred, y = y.pred, z = z.pred, col="gray40", facets=NA, lwd=0.25),
                        xlab="Backscatter (dB)", ylab="Var 1", zlab="Error (dB)", 
                        theta=45, phi=0, 
                        clab="Error", 
                        cex=0.5, pch=16,
                        col=ramp.col(c('red2', 'lightgoldenrod', 'green4'), alpha=0.25),
                        ticktype = "detailed", 
                        plot=FALSE)
              
              #plot in an RGL window
              plotrgl()
            }
            
          } else warning("Too many dimensions; error not plotted.", call.=FALSE, immediate.=TRUE)
        }
      }
      
      #for BRT corrections
      if(shift.method == "brt"){
        
        #create data.frame of predictors and fit the model using default arguments or those passed from the initial call
        cat("Fitting BRT model... \n")
        mod.vars <- back_data[2:length(back_data)]
        error.model <- gbm.step(data=mod.vars,
                                gbm.x = 1:(length(mod.vars)-1),
                                gbm.y = length(mod.vars),
                                family = "gaussian",
                                tree.complexity = tc,
                                learning.rate = lr,
                                bag.fraction = bf,
                                n.folds = 10,
                                max.trees = 20000,
                                plot.main = FALSE,
                                verbose = verb,
                                silent = !verb,
                                keep.fold.vector = FALSE,
                                keep.fold.fit = FALSE)
        
        #predict corrections to the working data.frame
        back_data <- data.frame(back_data, shifted=(predict(error.model, back_data, n.trees = error.model$n.trees, type = "response")) + back_data$shift)
        
        #align and stack the covariate rasters to the extent of the shift layer for predicting 
        cat("Predicting new raster... \n")
        pred.stack <- stack()
        for (i in 1:length(preds)){
          pred.mask <- mask(resample(preds[[i]], shift), shift)
          pred.stack <- stack(pred.stack, pred.mask)
        }
        
        pred.stack <- stack(mask(shift, pred.stack[[1]]), pred.stack)
        names(pred.stack) <- names(mod.vars[1:(length(mod.vars)-1)])
        
        #predict the error using the covariates in the stack, and add that to the shift layer
        shifted <- (predict(pred.stack, error.model, n.trees=error.model$n.trees, type="response")) + shift
        
        #for 3d error plots when there are two predictors (e.g., the shift layer and one covariate)
        if(err.plots == "TRUE") {
          if(length(preds)==1){
            
            cat("Plotting error... \n")
            
            #checking for required packages; installing if not present
            required.packages <- c("plot3D", "plot3Drgl")
            for(i in 1:length(required.packages)){
              if(!(required.packages[i] %in% installed.packages())) install.packages(required.packages[i])
            }
            
            require(plot3D)
            require(plot3Drgl)
            
            #create grids of predictors for plotting the model in 3d
            grid.lines <- 100
            x.pred <- seq(min(mod.vars[1]), max(mod.vars[1]), length.out = grid.lines)
            y.pred <- seq(min(mod.vars[2]), max(mod.vars[2]), length.out = grid.lines)
            xy <- expand.grid(x.pred, y.pred)
            names(xy) <- names(c(mod.vars[1],mod.vars[2]))
            
            #predict z-values of the and y grids using the model
            if(shift.method == "brt"){
              z.pred <- matrix(predict(error.model, xy, n.trees = error.model$n.trees, type = "response"), 
                               nrow = grid.lines, ncol = grid.lines)
            } else {
              z.pred <- matrix(predict(error.model, xy), 
                               nrow = grid.lines, ncol = grid.lines)
            }
            
            if(save.err.plots){
              
              plot.new()
              #create a 3d plot of the data with predictors on x- and y-axes, and error on z; plot the model as a surface
              scatter3D(mod.vars[[1]], mod.vars[[2]], mod.vars[[3]], 
                        bty="g", surf=list(x = x.pred, y = y.pred, z = z.pred, col="gray40", facets=NA, lwd=0.25),
                        xlab="Backscatter (dB)", ylab="Var 1", zlab="Error (dB)", 
                        theta=45, phi=0, 
                        clab="Error", 
                        cex=0.5, pch=16,
                        col=ramp.col(c('red2', 'lightgoldenrod', 'green4'), alpha=0.25),
                        ticktype = "detailed", 
                        plot=TRUE)
              
              plotrgl()
              
              p1 <- recordPlot()
            } else {
              
              #create a 3d plot of the data with predictors on x- and y-axes, and error on z; plot the model as a surface
              scatter3D(mod.vars[[1]], mod.vars[[2]], mod.vars[[3]], 
                        bty="g", surf=list(x = x.pred, y = y.pred, z = z.pred, col="gray40", facets=NA, lwd=0.25),
                        xlab="Backscatter (dB)", ylab="Var 1", zlab="Error (dB)", 
                        theta=45, phi=0, 
                        clab="Error", 
                        cex=0.5, pch=16,
                        col=ramp.col(c('red2', 'lightgoldenrod', 'green4'), alpha=0.25),
                        ticktype = "detailed", 
                        plot=FALSE)
              
              #plot in an RGL window
              plotrgl()
            }
            
          } else warning("Too many dimensions; error not plotted.", call.=FALSE, immediate.=TRUE)
        }
      }
      
    } else {
      
      #for corrections with no additional covariates (i.e., using the shift layer only; preds=NULL)
      
      back_data <- data.frame(as.data.frame(target.mask), as.data.frame(shift.mask))
      names(back_data) <- c("target", "shift")
      
      #remove NAs from the working data.frame, create a dummy variable in case shift.method="brt", calculate the error between layers
      back_data <- back_data[complete.cases(back_data),]
      back_data <- data.frame(back_data, dummy=1, err=back_data$target-back_data$shift)
      
      #for subsampling large datasets
      if(sample==TRUE){
        if((samp.size*nrow(back_data))>nrow(back_data)){
          stop("Subsample size is greater than number of raster cells; use a proportion < 1.")
        } else back_data <- back_data[seq(1, nrow(back_data), length.out=samp.size*nrow(back_data)), ]
      }
      
      #For mean corrections
      if(shift.method == "mean") {
        
        #shift data in working data.frame
        cat("Shifting layer by the mean of the error... \n")
        error.model <- mean(back_data$err)
        back_data <- data.frame(back_data, shifted=mean(back_data$err) + back_data$shift)
        
        #apply the same correction to the original shift raster
        cat("Predicting new raster... \n")
        shifted <- shift + mean(back_data$err)
        
        #for error plotting
        if(err.plots == "TRUE"){
          
          cat("Plotting error... \n")
          
          if(save.err.plots){
            
            p1 <- ggplot(back_data, aes(x=shift, y=err)) + 
              theme_minimal(11) +
              geom_point(col="black", size=0.5, alpha=0.3, shape=16) +
              geom_hline(aes(yintercept=error.model, lty="Error model"), col="red", alpha=1, lwd=1) +
              scale_linetype(name=NULL) +
              xlab("Shift backscatter (dB)") + 
              ylab("Error (dB)") +
              theme(axis.text = element_text(color = "black", family = "sans"), legend.position=c(0.89,0.885), panel.border=element_rect(fill=NA, colour="gray80", size=1.1), legend.text = element_text(color = "black", family = "sans", size=12), legend.background = element_blank())
            
            plot(p1)
          } else {
            
            plot(ggplot(back_data, aes(x=shift, y=err)) + 
                   theme_minimal(11) +
                   geom_point(col="black", size=0.5, alpha=0.3, shape=16) +
                   geom_hline(aes(yintercept=error.model, lty="Error model"), col="red", alpha=1, lwd=1) +
                   scale_linetype(name=NULL) +
                   xlab("Shift backscatter (dB)") + 
                   ylab("Error (dB)") +
                   theme(axis.text = element_text(color = "black", family = "sans"), legend.position=c(0.89,0.885), panel.border=element_rect(fill=NA, colour="gray80", size=1.1), legend.text = element_text(color = "black", family = "sans", size=12), legend.background = element_blank())
            )
          }
        }
      }
      
      #For lm corrections
      if(shift.method == "lm"){
        
        #fit the model
        cat("Fitting simple regression... \n")
        error.model <- glm(err~shift, family=family(link=link), data=back_data)
        
        #predict the shifted raster; predict corrections to the working data.frame
        cat("Predicting new raster... \n")
        shifted <- (predict(shift, error.model)) + shift
        back_data <- data.frame(back_data, shifted=(predict(error.model, back_data, type = "response")) + back_data$shift)
        
        #for plotting error between target and shift layers, with the bulk shift model
        if(err.plots == "TRUE"){
          
          cat("Plotting error... \n")
          
          if(save.err.plots){
            
            p1 <- ggplot(back_data, aes(x=shift, y=err)) + 
              theme_minimal(11) +
              geom_point(col="black", size=0.5, alpha=0.3, shape=16) +
              geom_line(aes(y=shifted-shift, lty="Error model"), col="red", alpha=1, lwd=1) +
              scale_linetype(name=NULL) +
              xlab("Shift backscatter (dB)") + 
              ylab("Error (dB)") +
              theme(axis.text = element_text(color = "black", family = "sans"), legend.position=c(0.89,0.885), panel.border=element_rect(fill=NA, colour="gray80", size=1.1), legend.text = element_text(color = "black", family = "sans", size=12), legend.background = element_blank())
            
            plot(p1)
          } else {
            
            plot(ggplot(back_data, aes(x=shift, y=err)) + 
                   theme_minimal(11) +
                   geom_point(col="black", size=0.5, alpha=0.3, shape=16) +
                   geom_line(aes(y=shifted-shift, lty="Error model"), col="red", alpha=1, lwd=1) +
                   scale_linetype(name=NULL) +
                   xlab("Shift backscatter (dB)") + 
                   ylab("Error (dB)") +
                   theme(axis.text = element_text(color = "black", family = "sans"), legend.position=c(0.89,0.885), panel.border=element_rect(fill=NA, colour="gray80", size=1.1), legend.text = element_text(color = "black", family = "sans", size=12), legend.background = element_blank())
            )
          }
        }
      }
      
      #For BRT corrections
      if(shift.method == "brt"){
        
        #fit the model using default arguments or those passed from the initial call
        cat("Fitting BRT model... \n")
        error.model <- gbm.step(data=back_data,
                                gbm.x = c(2,3),
                                gbm.y = 4,
                                family = "gaussian",
                                tree.complexity = tc,
                                learning.rate = lr,
                                bag.fraction = bf,
                                n.folds = 10,
                                max.trees = 20000,
                                plot.main = FALSE,
                                verbose = verb,
                                silent = !verb,
                                keep.fold.vector = FALSE,
                                keep.fold.fit = FALSE)
        
        #predict corrections to the working data.frame
        back_data <- data.frame(back_data, shifted=(predict(error.model, back_data, n.trees = error.model$n.trees, type = "response")) + back_data$shift)
        
        #create a raster stack using a dummy layer with the shift layer for predicting 
        cat("Predicting new raster... \n")
        dummy <- mask(shift^0, shift)
        pred.stack <- stack(shift, dummy)
        names(pred.stack) <- c("shift", "dummy")
        
        #predict the error using the raster the stack, and add that to the shift layer
        shifted <- (predict(pred.stack, error.model, n.trees=error.model$n.trees, type="response")) + shift
        
        #for plotting error between target and shift layers, with the bulk shift model
        if(err.plots == "TRUE"){
          
          cat("Plotting error... \n")
          
          if(save.err.plots){
            
            p1 <- ggplot(back_data, aes(x=shift, y=err)) + 
              theme_minimal(11) +
              geom_point(col="black", size=0.5, alpha=0.3, shape=16) +
              geom_line(aes(y=shifted-shift, lty="Error model"), col="red", alpha=1, lwd=1) +
              scale_linetype(name=NULL) +
              xlab("Shift backscatter (dB)") + 
              ylab("Error (dB)") +
              theme(axis.text = element_text(color = "black", family = "sans"), legend.position=c(0.89,0.885), panel.border=element_rect(fill=NA, colour="gray80", size=1.1), legend.text = element_text(color = "black", family = "sans", size=12), legend.background = element_blank())
            
            plot(p1)
          } else {
            
            plot(ggplot(back_data, aes(x=shift, y=err)) + 
                   theme_minimal(11) +
                   geom_point(col="black", size=0.5, alpha=0.3, shape=16) +
                   geom_line(aes(y=shifted-shift, lty="Error model"), col="red", alpha=1, lwd=1) +
                   scale_linetype(name=NULL) +
                   xlab("Shift backscatter (dB)") + 
                   ylab("Error (dB)") +
                   theme(axis.text = element_text(color = "black", family = "sans"), legend.position=c(0.89,0.885), panel.border=element_rect(fill=NA, colour="gray80", size=1.1), legend.text = element_text(color = "black", family = "sans", size=12), legend.background = element_blank())
            )
          }
        }
      }
    }
    
    #calculate evaluation statistics for the bulk shift
    fit.stats <- c(init_MAE=mean(abs(back_data$target-back_data$shift)),
                   shifted_MAE=mean(abs(back_data$target-back_data$shifted)),
                   init_KS=suppressWarnings(ks.test(back_data$target, back_data$shift)$statistic),
                   shifted_KS=suppressWarnings(ks.test(back_data$target, back_data$shifted)$statistic)
    )
    
    #for plotting data distributions (pdf and cdf) resulting from the bulk shift
    if(dist.plots == TRUE){
      
      #checking for required packages; installing if not present
      required.packages <- c("reshape2", "cowplot")
      for(i in 1:length(required.packages)){
        if(!(required.packages[i] %in% installed.packages())) install.packages(required.packages[i])
      }
      require(reshape2)
      require(cowplot)
      
      cdf <- data.frame(Target=sort(back_data$target), Shift=sort(back_data$shift), 
                        Shifted=sort(back_data$shifted), fn.x=(1:length(back_data$shifted))/length(back_data$shifted))
      
      cdf.melt <- melt(cdf, id.vars="fn.x")
      
      cdf.plot <- ggplot(cdf.melt, aes(x=value, y=fn.x, col=variable)) +
        geom_line(lwd=0.75) +
        scale_color_manual(values=c("black", "green2", "dodgerblue")) +
        xlab(NULL) +
        ylab("CDF") +
        theme_minimal(11) +
        coord_cartesian(xlim=c(min(cdf.melt$value), max(cdf.melt$value)), expand=FALSE) +
        theme(axis.text = element_text(color = "black", family = "sans"), legend.position=c(0.9,0.5), panel.border=element_rect(fill=NA, colour="gray80", size=1.1), legend.text = element_text(color = "black", family = "sans", size=12), legend.background = element_blank(), legend.title=element_blank())
      
      pdf.plot <- ggplot() +
        geom_line(aes(x=density(back_data$target)$x, y=density(back_data$target)$y), lwd=0.75, col="black") +
        geom_line(aes(x=density(back_data$shift)$x, y=density(back_data$shift)$y), lwd=0.75, col="green2") +
        geom_line(aes(x=density(back_data$shifted)$x, y=density(back_data$shifted)$y), lwd=0.75, col="dodgerblue") +
        xlab("Backscatter (dB)") +
        ylab("PDF") +
        theme_minimal(11) +
        coord_cartesian(xlim=c(min(cdf.melt$value), max(cdf.melt$value)), expand=FALSE) +
        theme(axis.text = element_text(color = "black", family = "sans"), panel.border=element_rect(fill=NA, colour="gray80", size=1.1))
      
      if(save.dist.plots==TRUE){
        
        distributions <- plot_grid(cdf.plot, pdf.plot, ncol = 1, align="v", rel_heights = c(1,1))
        plot(distributions)
      } else {
        
        plot(plot_grid(cdf.plot, pdf.plot, ncol = 1, align="v", rel_heights = c(1,1)))
      }
    }
    
    #collate and return the results of the bulk shift
    if(save.data == TRUE){
      if(mosaic == TRUE){
        cat("Creating harmonized mosaic using", mosaic.method, "interpolation... \n")
        shifted.ext.re <- resample(shifted, extend(target, shifted), method=mosaic.method)
        mosaic <- mosaic(shifted.ext.re, target, fun=mean)
        
        result <- list(shifted=shifted, mosaic=mosaic, model=error.model, fit_stats=fit.stats, data=back_data)
      } else {
        
        result <- list(shifted=shifted, model=error.model, fit_stats=fit.stats, data=back_data)
      }
    } else {
      
      if(mosaic == TRUE){
        cat("Creating harmonized mosaic using", mosaic.method, "interpolation... \n")
        shifted.ext.re <- resample(shifted, extend(target, shifted), method=mosaic.method)
        mosaic <- mosaic(shifted.ext.re, target, fun=mean)
        
        result <- list(shifted=shifted, mosaic=mosaic, model=error.model, fit_stats=fit.stats)
      } else {
        
        result <- list(shifted=shifted, model=error.model, fit_stats=fit.stats)
      }
    }
    
    if(save.dist.plots & exists("distributions", inherits=FALSE)) result[["dist_plots"]] <- distributions
    if(save.err.plots & exists("p1", inherits=FALSE)) result[["err_plot"]] <- p1
    
    cat("Complete. Output list contains: ")
    cat(names(result), sep=", ")
    
    return(result)
  }