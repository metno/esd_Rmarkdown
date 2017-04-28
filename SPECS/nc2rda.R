## Generate rda-files for the decadal forecasts - organise the data in the same way as the GFDL data
##> summary(psl_DECA)
##                    Length   Class  Mode     
##Variable                   2 -none- list     
##Data                79704000 -none- numeric  
##xyCoords                   2 -none- list     
##Dates                      2 -none- list     
##InitializationDates        9 -none- list     
##Members                   10 -none- character
##
##> dim(psl_DECA$Data)
## [1]   9  10 480  41  45
##      ini mem tim lat lon
#> range(tas_DECA$xyCoords$x)
#[1] -88.75  21.25
#> range(tas_DECA$xyCoords$y)
#[1] -1.011236 79.887640
#> 
  
## Each netCDF file: float tas(time, ensemble, latitude, longitude) - global
#dimensions:
#  latitude = 73 ;
#  longitude = 96 ;
#  time = UNLIMITED ; // (120 currently)
#  ensemble = 10 ;
  
library(ncdf4)

getvar <- function(fname,varid,lons=c(-88.75, 21.25),lats=c(-1.0, 80.0)) {
  ## Extract the region
  ncid <- nc_open(fname)
  
  lon <- ncvar_get(ncid,'longitude')
  lat <- ncvar_get(ncid,'latitude')
  ens <- ncvar_get(ncid,'ensemble')
  tim <- ncvar_get(ncid,'time')
  ix <- (1:length(lon))[(lon >= lons[1]) & (lon <= lons[2])]
  iy <- (1:length(lat))[(lat >= lats[1]) & (lat <= lats[2])]
  lon <- lon[ix]; lat <- lat[iy]
  x <- ncvar_get(ncid,varid,start=c(1,1,min(iy),min(ix)),count(length(tim),ength(ens),length(iy),length(ix)))
  
}

varid <- 'tas'
files <- list.files(path='~/Downloads/SPECS-decadal-fcst_hadcm3_i2p1',pattern=varid,full.names=TRUE)
#files <- list.files(path='~/Downloads/SPECS-decadal-fcst_MPI',pattern=varid,full.names=TRUE)
