# bulkshift
### Backscatter raster harmonization for seabed mapping in R

This package facilitates relative calibration for processed sonar backscatter raster datasets using the methods described in Misiuk et al. (2020). 
It also extends these methods by providing functionality for advanced subsampling to handle dataset imbalance and spatial autocorrelation, and support for additional modelling methods. 

The functions in this package rely heavily on the [terra](https://github.com/rspatial/terra) package. RasterLayers are still supported, but the use of terra is highly 
recommended - it is meant to replace the [raster](https://github.com/rspatial/raster) package.

## Installation
You can use the remotes package to install bulkshift from github. Install remotes first if you do not have it. 
```
install.packages("remotes")
remotes::install_github("benjaminmisiuk/bulkshift")
```
## Basic functionality
Load the bulkshift and terra libraries. The terra package should have installed along with bulkshift. You can install
it explicitly using `install.packages("terra")`. While terra is not required to load bulkshift, you will need it
to read rasters into R.
```
library(terra)
library(bulkshift)
```
Example datasets from the Bedford Basin (Misiuk et al., 2020) are provided as part of the package. These are backscatter
datasets collected during two different surveys, and also a depth raster.
```
bb2016 <- rast(system.file('extdata', 'bb2016.tif', package='bulkshift'))
bb2017 <- rast(system.file('extdata', 'bb2017.tif', package='bulkshift'))
bbdepth <- rast(system.file('extdata', 'bbdepth.tif', package='bulkshift'))
```
These datasets are not calibrated, and mosaicking them as-is produces obvious artefacts where the values do not match between surveys:
```
plot(mosaic(bb2016, bb2017), col = gray.colors(100))
```
![](images/bshift_eg1.png)

At its most basic, the bulk shift performs relative calibration using mutual overlap between surveys:
```
b <- bulkshift(shift = bb2017, target = bb2016)
```
The result is a list with several objects, which may differ depending on the function arguments (e.g., `mosaic`, `savedata`, `crossvalidate`). We can inspect the model that was used to calibrate bb2017 just like any other model in R:
```
b$errorModel

#Call:  glm(formula = form, data = df)
#
#Coefficients:
#(Intercept)        shift  
#     3.7043      -0.2301  
#
#Degrees of Freedom: 16419 Total (i.e. Null);  16418 Residual
#Null Deviance:	    132200 
#Residual Deviance: 112300 	AIC: 78180
```
We can see that the model is a GLM. Recall that the response variable is the error between datasets, therefore, the coefficient for "shift" describes the change in error per unit of the "shift" backscatter dataset. The corrected backscatter layer can be plotted:
```
plot(b$shifted)
```
![](images/bshift_eg2.png)

The bathymetry layer can be added as a covariate in the model, and we can create a mosaic of the backscatter layers. All additional covariates must be included as layers in a single SpatRaster and passed to the `preds` argument:
```
b <- bulkshift(shift = bb2017, target = bb2016, preds = bbdepth, mosaic = TRUE)
plot(b$mosaic, col = gray.colors(100))
```
![](images/bshift_eg3.png)

We can then observe statistics such as the variance explained (VE), mean absolute error (MAE), and Pearson correlation between the datasets before and after correction. If a proportion of data were set aside for validation using the `crossvalidate` argument, test statistics will additionally be returned. The model data are also returned by default, and relationships between backscatter datasets before and after correction can be compared. These data can additionally be used to calculate your own bespoke validation statistics, or to apply modelling methods that are not implemented in the package.
```
b$fitStats
#     layer         VE      MAE         r
#1  Original -2.1262196 7.469381 0.8151604
#2 Corrected  0.8556397 1.268409 0.9250079

par(mfrow = c(1,2))
plot(b$data$target, b$data$shift)
plot(b$data$target, b$data$shifted)
```
![](images/bshift_eg4.png)

# References
Misiuk, B., Brown, C.J., Robert, K., Lacharite, M., 2020. Harmonizing Multi-Source Sonar Backscatter Datasets for Seabed Mapping Using Bulk Shift Approaches. Remote Sensing 12, 601. https://doi.org/10.3390/rs12040601
