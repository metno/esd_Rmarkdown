### Present maps of downscaled scenarios for AACA - the Barents region
### temperature and precipitation.
### Based on PCA-representation for the predictand.
### Use 'pcafill' to fill in missing values.
### Rasmus E. Benestad, Oslo, November, 2015

library(esd)
library(ncdf4)

distill <- function(x,it=2050,probs=c(0.05,0.95)) {
  y <- subset(x,it=it)
  attributes(y) <- NULL
  y <- quantile(y,probs=probs,na.rm=TRUE)
  return(y)
}


grid2ncdf <- function(x,fname='grid2ncdf.nc') {
  print('grid2ncdf')
  dimlon <- ncdim_def( "longitude", "degree_east", lon(x) )
  dimlat <- ncdim_def( "latitude", "degree_north", lat(x) )
  x4nc <- ncvar_def(varid(x)[1], unit(x)[1], list(dimlon,dimlat), -1, 
                    longname=attr(x,'longname'), prec="integer")
     
     # Create a netCDF file with this variable
  ncnew <- nc_create( fname, x4nc )

  y <- coredata(x); attributes(y) <- NULL
  meany <- mean(y,na.rm=TRUE)
  y <- round(100*(y-meany))
  y[!is.finite(y)] <- -99
  # Write some values to this variable on disk.
  ncvar_put( ncnew, x4nc, round(y) )
  ncatt_put( ncnew, x4nc, "add_offset", meany, prec="float" )
  ncatt_put( ncnew, x4nc, "scale_factor", 0.01, prec="float" ) 
  ncatt_put( ncnew, x4nc, "_FillValue", -99, prec="float" ) 
  ncatt_put( ncnew, x4nc, "missing_value", -99, prec="float" ) 
  ncatt_put( ncnew, 0, "description", 
             "Gridded downscaled results using LaticeKrig")
  nc_close(ncnew)
}

propchange <- function(X,it0=c(1986,2015)) {
  coredata(X) <- 100*t(t(coredata(X))/apply(coredata(subset(X,it=it0)),2,'mean',na.rm=TRUE))
  X
}

gridmap <- function(fname='dse.aaca.t2m.rcp45.djf.rda',param=NULL,
                    fast=FALSE,it=2050,it0=2010,FUN='change',
                    probs=c(0.05,0.95),verbose=FALSE,prop=FALSE,
                    breaks=NULL,pal=NULL,col=NULL,rev=TRUE,scal=1) {

  if (verbose) print('gridmap')
  library(LatticeKrig)
  if (is.character(fname)) {
    load(fname)

  dse.results <- switch(substr(fname,1,18),
                       'dse.aaca.t2m.rcp26'=dse.aaca.t2m.rcp26,
                       'dse.aaca.pr.rcp26.'=dse.aaca.pr.rcp26,
                       'dse.aaca.t2m.rcp45'=dse.aaca.t2m.rcp45,
                       'dse.aaca.pr.rcp45.'=dse.aaca.pr.rcp45,
                       'dse.aaca.t2m.rcp85'=dse.aaca.t2m.rcp85,
                       'dse.aaca.pr.rcp85.'=dse.aaca.pr.rcp85)
  } else if (inherits(fname,'list')) dse.results <- fname
  
  if (is.null(param)) param <- varid(dse.results) 
  dse.station <- as.station(dse.results,verbose=verbose)
  
  if (prop) dse.station <- lapply(dse.station,propchange)

  Z0 <- lapply(dse.station,distill,it=it0,probs=probs)
  Z <- lapply(dse.station,distill,it=it,probs=probs)
  lon <- unlist(lapply(dse.station,lon))
  lat <- unlist(lapply(dse.station,lat))
  alt <- unlist(lapply(dse.station,function(x) alt(attr(x,'station'))))
  X <- attr(dse.station[[1]],'station')
  z0 <- matrix(unlist(Z0),length(probs),length(dse.station))*scal
  z <- matrix(unlist(Z),length(probs),length(dse.station))*scal
  data(etopo5)
  etopo5 <- subset(etopo5,is=list(lon=range(lon),lat=range(lat)))
  etopo5[etopo5<=0] <- NA

  ## Set the grid to be the same as that of etopo5:
  grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')

  #3 Flag dubplicated stations:
  ok <- !(duplicated(lon) & duplicated(lat))

  ## Spread in the  90-percente interval changing
  obj <- switch( FUN,
                'change'=LatticeKrig( x=cbind(lon[ok],lat[ok]),
                                      y=z[1,ok] - z0[1,ok],Z=alt[ok]),
                'spread'=LatticeKrig( x=cbind(lon[ok],lat[ok]),
                     y=100*(z[2,ok] - z[1,ok])/(z0[2,ok] - z0[1,ok]),Z=alt[ok]) )

  ##  obj <- LatticeKrig( x=cbind(lon[ok],lat[ok]), y=z[2,ok],Z=alt[ok])
  w <- predictSurface(obj, grid.list = grid,Z=etopo5)
  w$z[is.na(etopo5)] <- NA
  #surface( w )
  #points(cbind(lon[ok],lat[ok]))
  
  ## coldspells(x,dse=NULL,it='djf',threshold=0,
  ##                     verbose=FALSE,plot=TRUE,...)

  
  ## coldwinterdays(,dse=NULL,it='djf',threshold=0,
  ##                  verbose=FALSE,plot=TRUE,nmin=90,...)

  ## Get rid of packages that have functions of same name:
  detach("package:LatticeKrig")
  detach("package:fields")
  detach("package:spam")
  detach("package:grid")
  detach("package:maps")

  ## Convert the results from LatticeKrig to esd:
  W <- w$z
  W[W > max(breaks)] <- max(breaks)
  W[W < min(breaks)] <- min(breaks)
  attr(W,'variable') <- param
  attr(W,'unit') <- unit(dse.results)
  attr(W,'longitude') <- w$x
  attr(W,'latitude') <- w$y
  class(W) <- class(etopo5)

  ## Make a projection that zooms in on the Barents region
  dev.new()

  colbar <- list(breaks=breaks,rev=rev,pal=pal,col=col)
  if (verbose) print(colbar)
  if (is.null(rev)) rev <- switch(param,'t2m'=FALSE,'pr'=TRUE)
  Wx <- max(abs(W),na.rm=TRUE)
  if (is.null(breaks)) breaks <- round(seq(-Wx,Wx,length=31),2)
  if (fast)
    map(W,xlim=range(lon(W)),ylim=range(lat(W)),
        colbar=colbar) else {
    attr(W,'variable') <- NULL
    attr(W,'unit') <- NULL
    map(W,xlim=range(lon(W)),ylim=range(lat(W)),projection='sphere',
        colbar=colbar,verbose=verbose)
  }
  figlab(it)
  dev.copy2pdf(file=paste('dse_paper59_Fig_ESD_',
                           param,'_',FUN,'_',it,'map.pdf',sep=''))

  attr(W,'variable') <- param
  attr(W,'unit') <- unit(dse.results)
  attr(W,'longname') <- '2-meter temperature'
  invisible(W)
}

## Use LatticeKrig to grid the results - eg extract the 5 and 95 percentiles
## of projected temperature for 2050 from the ensembles.


breaks <- round(seq(-4,18,by=1),2)
col <- colscal(37)[16:37]

if (FALSE) {
  print('RCP4.5')
gridmap('dse.aaca.t2m.rcp45.djf.rda',it=2050,it0=2010,FUN='change',
         probs=c(0.95),breaks=breaks,col=col,verbose=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.djf.2050.pdf')

gridmap('dse.aaca.t2m.rcp45.djf.rda',it=2030,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.djf.2030.pdf')

gridmap('dse.aaca.t2m.rcp45.djf.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.djf.2100.pdf')

gridmap('dse.aaca.t2m.rcp45.jja.rda',it=2050,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.jja.2050.pdf')

gridmap('dse.aaca.t2m.rcp45.jja.rda',it=2030,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.jja.2030.pdf')

gridmap('dse.aaca.t2m.rcp45.jja.rda',it=2099,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.jja.2100.pdf')

gridmap('dse.aaca.t2m.rcp45.mam.rda',it=2050,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.mam.2050.pdf')

gridmap('dse.aaca.t2m.rcp45.mam.rda',it=2030,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.mam.2030.pdf')

gridmap('dse.aaca.t2m.rcp45.mam.rda',it=2099,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.mam.2100.pdf')

gridmap('dse.aaca.t2m.rcp45.son.rda',it=2050,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.son.2050.pdf')

gridmap('dse.aaca.t2m.rcp45.son.rda',it=2030,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.son.2030.pdf')

gridmap('dse.aaca.t2m.rcp45.son.rda',it=2099,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp45.son.2100.pdf')
}

if (FALSE) {
  print('RCP8.5')
gridmap('dse.aaca.t2m.rcp85.djf.rda',it=2050,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.djf.2050.pdf')

gridmap('dse.aaca.t2m.rcp85.djf.rda',it=2030,it0=2010,FUN='change',
        probs=c(0.95),breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.djf.2030.pdf')

gridmap('dse.aaca.t2m.rcp85.djf.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.djf.2100.pdf')

gridmap('dse.aaca.t2m.rcp85.jja.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.jja.2050.pdf')

gridmap('dse.aaca.t2m.rcp85.jja.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.jja.2030.pdf')

gridmap('dse.aaca.t2m.rcp85.jja.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.jja.2100.pdf')

gridmap('dse.aaca.t2m.rcp85.mam.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.mam.2050.pdf')

gridmap('dse.aaca.t2m.rcp85.mam.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.mam.2030.pdf')

gridmap('dse.aaca.t2m.rcp85.mam.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.mam.2100.pdf')

gridmap('dse.aaca.t2m.rcp85.son.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.son.2050.pdf')

gridmap('dse.aaca.t2m.rcp85.son.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.son.2030.pdf')

gridmap('dse.aaca.t2m.rcp85.son.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp85.son.2100.pdf')
}

if (FALSE) {
  print('RCP2.6')
gridmap('dse.aaca.t2m.rcp26.djf.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.djf.2050.pdf')

gridmap('dse.aaca.t2m.rcp26.djf.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.djf.2030.pdf')

gridmap('dse.aaca.t2m.rcp26.djf.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.djf.2100.pdf')

gridmap('dse.aaca.t2m.rcp26.jja.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.jja.2050.pdf')

gridmap('dse.aaca.t2m.rcp26.jja.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.jja.2030.pdf')

gridmap('dse.aaca.t2m.rcp26.jja.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.jja.2100.pdf')

gridmap('dse.aaca.t2m.rcp26.mam.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.mam.2050.pdf')

gridmap('dse.aaca.t2m.rcp26.mam.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.mam.2030.pdf')

gridmap('dse.aaca.t2m.rcp26.mam.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.mam.2100.pdf')

gridmap('dse.aaca.t2m.rcp26.son.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.son.2050.pdf')

gridmap('dse.aaca.t2m.rcp26.son.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.son.2030.pdf')

gridmap('dse.aaca.t2m.rcp26.son.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,rev=TRUE)
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.t2m.rcp26.son.2100.pdf')
}


breaks <- round(seq(-1.5,1.5,by=0.1),2)
col <- colscal(length(breaks)-1,rev=FALSE)

if (FALSE) {
gridmap('dse.aaca.pr.rcp26.wetmean.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,verbose=TRUE)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.mu.2050.pdf')

gridmap('dse.aaca.pr.rcp26.wetmean.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.mu.2030.pdf')

gridmap('dse.aaca.pr.rcp26.wetmean.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.mu.2100.pdf')

gridmap('dse.aaca.pr.rcp45.wetmean.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.mu.2050.pdf')

gridmap('dse.aaca.pr.rcp45.wetmean.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.mu.2030.pdf')

gridmap('dse.aaca.pr.rcp45.wetmean.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.mu.2100.pdf')

gridmap('dse.aaca.pr.rcp85.wetmean.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.mu.2050.pdf')

gridmap('dse.aaca.pr.rcp85.wetmean.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.mu.2030.pdf')

gridmap('dse.aaca.pr.rcp85.wetmean.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.mu.2100.pdf')
}

if (FALSE) {
gridmap('dse.aaca.pr.rcp26.wetfreq.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,verbose=TRUE)
figlab(expression(f[w]*(fraction)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.fw.2050.pdf')

gridmap('dse.aaca.pr.rcp26.wetfreq.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(f[w]*(fraction)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.fw.2030.pdf')

gridmap('dse.aaca.pr.rcp26.wetfreq.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(f[w]*(fraction)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.fw.2100.pdf')

gridmap('dse.aaca.pr.rcp45.wetfreq.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(f[w]*(fraction)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.fw.2050.pdf')

gridmap('dse.aaca.pr.rcp45.wetfreq.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(f[w]*(fraction)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.fw.2030.pdf')

gridmap('dse.aaca.pr.rcp45.wetfreq.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(f[w]*(fraction)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.fw.2100.pdf')

gridmap('dse.aaca.pr.rcp85.wetfreq.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(f[w]*(fraction)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.fw.2050.pdf')

gridmap('dse.aaca.pr.rcp85.wetfreq.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(f[w]*(fraction)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.fw.2030.pdf')

gridmap('dse.aaca.pr.rcp85.wetfreq.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(f[w]*(fraction)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.fw.2100.pdf')
}


breaks <- round(seq(-25,25,by=1),2)
col <- colscal(length(breaks)-1,rev=FALSE)

if (FALSE) {
gridmap('dse.aaca.pr.rcp26.amjjas.wetmean.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,verbose=TRUE)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.amjjas.mu.2050.pdf')

gridmap('dse.aaca.pr.rcp26.amjjas.wetmean.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.amjjas.mu.2030.pdf')

gridmap('dse.aaca.pr.rcp26.amjjas.wetmean.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.amjjas.mu.2100.pdf')

gridmap('dse.aaca.pr.rcp45.amjjas.wetmean.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.amjjas.mu.2050.pdf')

gridmap('dse.aaca.pr.rcp45.amjjas.wetmean.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.amjjas.mu.2030.pdf')

gridmap('dse.aaca.pr.rcp45.amjjas.wetmean.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.amjjas.mu.2100.pdf')

gridmap('dse.aaca.pr.rcp85.amjjas.wetmean.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.amjjas.mu.2050.pdf')

gridmap('dse.aaca.pr.rcp85.amjjas.wetmean.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.amjjas.mu.2030.pdf')

gridmap('dse.aaca.pr.rcp85.amjjas.wetmean.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.amjjas.mu.2100.pdf')
}


if (FALSE) {
gridmap('dse.aaca.pr.rcp26.jfmond.wetmean.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,verbose=TRUE)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.jfmond.mu.2050.pdf')

gridmap('dse.aaca.pr.rcp26.jfmond.wetmean.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.jfmond.mu.2030.pdf')

gridmap('dse.aaca.pr.rcp26.jfmond.wetmean.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.jfmond.mu.2100.pdf')

gridmap('dse.aaca.pr.rcp45.jfmond.wetmean.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.jfmond.mu.2050.pdf')

gridmap('dse.aaca.pr.rcp45.jfmond.wetmean.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.jfmond.mu.2030.pdf')

gridmap('dse.aaca.pr.rcp45.jfmond.wetmean.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.jfmond.mu.2100.pdf')

gridmap('dse.aaca.pr.rcp85.jfmond.wetmean.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.jfmond.mu.2050.pdf')

gridmap('dse.aaca.pr.rcp85.jfmond.wetmean.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.jfmond.mu.2030.pdf')

gridmap('dse.aaca.pr.rcp85.jfmond.wetmean.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col)
figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.jfmond.mu.2100.pdf')
}

#breaks <- round(seq(-30,30,by=1),2)
breaks <- round(seq(-10,10,by=1),2)
col <- colscal(length(breaks)-1,rev=FALSE)

if (FALSE) {
gridmap('dse.aaca.pr.rcp26.amjjas.wetfreq.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.amjjas.fw.2050.pdf')

gridmap('dse.aaca.pr.rcp26.amjjas.wetfreq.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.amjjas.fw.2030.pdf')

gridmap('dse.aaca.pr.rcp26.amjjas.wetfreq.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.amjjas.fw.2100.pdf')

gridmap('dse.aaca.pr.rcp45.amjjas.wetfreq.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.amjjas.fw.2050.pdf')

gridmap('dse.aaca.pr.rcp45.amjjas.wetfreq.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.amjjas.fw.2030.pdf')

gridmap('dse.aaca.pr.rcp45.amjjas.wetfreq.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.amjjas.fw.2100.pdf')

gridmap('dse.aaca.pr.rcp85.amjjas.wetfreq.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.amjjas.fw.2050.pdf')

gridmap('dse.aaca.pr.rcp85.amjjas.wetfreq.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.amjjas.fw.2030.pdf')

gridmap('dse.aaca.pr.rcp85.amjjas.wetfreq.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.amjjas.fw.2100.pdf')
}


if (FALSE) {
gridmap('dse.aaca.pr.rcp26.jfmond.wetfreq.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.jfmond.fw.2050.pdf')

gridmap('dse.aaca.pr.rcp26.jfmond.wetfreq.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.jfmond.fw.2030.pdf')

gridmap('dse.aaca.pr.rcp26.jfmond.wetfreq.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp26.jfmond.fw.2100.pdf')


gridmap('dse.aaca.pr.rcp45.jfmond.wetfreq.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.jfmond.fw.2050.pdf')

gridmap('dse.aaca.pr.rcp45.jfmond.wetfreq.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.jfmond.fw.2030.pdf')

gridmap('dse.aaca.pr.rcp45.jfmond.wetfreq.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp45.jfmond.fw.2100.pdf')

gridmap('dse.aaca.pr.rcp85.jfmond.wetfreq.rda',it=2050,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.jfmond.fw.2050.pdf')

gridmap('dse.aaca.pr.rcp85.jfmond.wetfreq.rda',it=2030,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.jfmond.fw.2030.pdf')

gridmap('dse.aaca.pr.rcp85.jfmond.wetfreq.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
        breaks=breaks,col=col,prop=TRUE)
figlab(expression(f[w]*('%')),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5',xpos = 0.8, ypos = 0.99)
dev.copy2pdf(file='paper59_dse.aaca.pr.rcp85.jfmond.fw.2100.pdf')
}

if (FALSE) {
  breaks <- round(seq(0,10,by=0.5),2)
  col <- colscal(20,col='warm')

  gridmap('dse.aaca.t2m.rcp45.djf.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
          breaks=breaks,col=col)
  figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
  figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
  figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
  dev.copy2pdf(file='fig2a.pdf')

  breaks <- round(seq(0,100,by=5),2)
  col <- colscal(length(breaks)-1,col='precip',rev=TRUE)

  gridmap('dse.aaca.t2m.rcp45.djf.rda',it=2099,it0=2010,FUN='change',probs=c(0.95),
          breaks=breaks,col=col,scal=8.75)
  figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
  figlab('95 percentile diff (2099 - 2010)',xpos = 0.001, ypos = 0.99)
  figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
  dev.copy2pdf(file='fig2b.pdf')
}


## AACA Chapter 4:
  breaks <- round(seq(0,10,by=0.5),2)
  col <- colscal(20,col='warm')

  W <- gridmap('dse.aaca.t2m.rcp45.djf.rda',it=2080,it0=2015,FUN='change',probs=c(0.95),
          breaks=breaks,col=col)
  figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
  figlab('95 percentile diff (2080 - 2015)',xpos = 0.001, ypos = 0.99)
  figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
  dev.copy2pdf(file='aacafig.esd1.pdf')
  dev2bitmap(file='aacafig.esd1.png',res=150, height = 14, width = 14)
  grid2ncdf(W,fname='t2m.rcp45.djf.2080-2015.nc')

if (FALSE) {
  breaks <- seq(0,100,by=5)
  col <- colscal(length(breaks)-1,col='precip',rev=TRUE)
  colbar <- list(col=col,breaks=breaks)
  source('mu2pt.R')

#  figlab(expression(paste(mu,' (%)')),xpos = 0.001, ypos = 0.1)
  figlab('95 percentile diff (2080 - 2015)',xpos = 0.001, ypos = 0.99)
  figlab('RCP4.5',xpos = 0.8, ypos = 0.99)
  dev.copy2pdf(file='aacafig.esd2.pdf')
  dev2bitmap(file='aacafig.esd2.png',res=150, height = 14, width = 14)
}
