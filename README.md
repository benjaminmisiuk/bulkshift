# bulkshift
### Bulk shift backscatter harmonization for seabed mapping in R

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

par(mfrow=c(2,1))
plot(bb2016, col = gray.colors(100))
plot(bb2017, col = gray.colors(100))
```
![](images/bshift_eg1.png)

At its most basic, the bulk shift performs relative calibration using mutual overlap between surveys:
```
b <- bulkshift(shift = bb2017, target = bb2016)
```
The result is a list with several important objects. We can inspect the model that was used to calibrate bb2017 just like
any other model in R:
```
b$errorModel
```
We can see that the model is a GLM, and the corrected backscatter layer can be plotted:
```
#reset the graphics device
dev.off()
plot(b$shifted)
```
![](images/bshift_eg2.png)

The bathymetry layer can be added as a covariate in the model, and we can create a mosaic of the backscatter layers. All additional covariates must be included as layers in a single SpatRaster and passed to the `preds` arggument:
```
bshift <- bulkshift(shift = bb2017, target = bb2016, preds = bbdepth, mosaic = TRUE)
plot(b$mosaic, col = gray.colors(100))
```
![](images/bshift_eg3.png)
# References
Misiuk, B., Brown, C.J., Robert, K., Lacharite, M., 2020. Harmonizing Multi-Source Sonar Backscatter Datasets for Seabed Mapping Using Bulk Shift Approaches. Remote Sensing 12, 601. https://doi.org/10.3390/rs12040601
