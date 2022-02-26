# bulkshift
### Bulk shift backscatter harmonization for seabed mapping in R

This package facilitates relative calibration for processed sonar backscatter raster datasets using the methods presented in Misiuk et al. (2020). 
It also extends these methods by providing functionality for advanced subsampling to handle dataset imbalance and spatial autocorrelation, and supports multiple modelling methods. 

The functions in this package rely heavily on the [terra](https://github.com/rspatial/terra) package. RasterLayers are still supported, but the use of terra is highly 
recommended. It is meant to replace the [raster](https://github.com/rspatial/raster) package.

# References
Misiuk, B., Brown, C.J., Robert, K., Lacharite, M., 2020. Harmonizing Multi-Source Sonar Backscatter Datasets for Seabed Mapping Using Bulk Shift Approaches. Remote Sensing 12, 601. https://doi.org/10.3390/rs12040601
