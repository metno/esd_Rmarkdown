### Downscale scenarios for AACA - the Barents region
### temperature and precipitation.
### Based on PCA-representation for the predictand.
### Use 'pcafill' to fill in missing values.

library(esd)

distill <- function(x,it=2050,probs=c(0.05,0.95)) {
  y <- subset(x,it=it)
  attributes(y) <- NULL
  y <- quantile(y,probs=probs,na.rm=TRUE)
  return(y)
}

downscale <- function(Y,predictor,it='djf',param='t2m',FUN='mean',
                      period=c(1950,2015),plot=FALSE,rcp='rcp45',
                      lon=c(0,100),lat=c(65,90)) {

## Use a time and space window:  
  Y <- subset(Y,it=period)
  Y <- subset(Y,is=list(lat=lat,lon=lon))

## Estimate seasonal means & weed out stations with little data
  Y4 <- subset(as.4seasons(Y,FUN=FUN),it=it)
  ok <- apply(coredata(Y4),1,nv)
  Y4 <- subset(Y4,it=ok>0)
  nok <- apply(coredata(Y4),2,nv)
  Y4 <- subset(Y4,is=nok>15)

  if (plot) map(Y,FUN=FUN,cex=-2)

## Fill missing data using PCA-based regression
  Z <- pcafill(Y4)

  pca <- PCA(Z)
  if (plot) plot(pca)

## Downscale results
  dse.pca <- DSensemble(pca,predictor=predictor,
                        biascorrect=TRUE,rcp=rcp)
  invisible(dse.pca)
}

diagnoseESD <- function(fname='dse.aaca.t2m.rcp45.djf.rda',is=1) {
  load(fname)
  dse.station <- as.station(dse.aaca.t2m.rcp45,verbose=TRUE)

  vis(dse.station,verbose=TRUE)
  dev.copy(png,paste(substr(fname,1,nchar(fname)-7),".vis.png",sep=''))
  dev.off()

  diagnose(dse.station,verbose=TRUE)
  dev.copy(png,paste(substr(fname,1,nchar(fname)-7),".diagnose.png",sep=''))
  dev.off()

  diagnose(dse.station[[is]],verbose=TRUE)
  dev.copy(png,paste(substr(fname,1,nchar(fname)-7),".s",is,".diagnose.png",sep=''))
  dev.off()

  plot(dse.station[[is]],verbose=TRUE)
  dev.copy(png,paste(substr(fname,1,nchar(fname)-7),".s",is,".plot.png",sep=''))
  dev.off()
}

gridmap <- function(fname='dse.aaca.t2m.rcp45.djf.rda',it=2050,it0=2010,
                    what='change',probs=c(0.05,0.95),breaks=NULL,pal=NULL) {

  print('gridmap')
  library(LatticeKrig)
  load(fname)
  dse.station <- as.station(dse.aaca.t2m.rcp45,verbose=TRUE)

  Z0 <- lapply(dse.station,distill,it=it0,probs=probs)
  Z <- lapply(dse.station,distill,it=it,probs=probs)
  lon <- unlist(lapply(dse.station,lon))
  lat <- unlist(lapply(dse.station,lat))
  alt <- unlist(lapply(dse.station,function(x) alt(attr(x,'station'))))
  X <- attr(dse.station[[1]],'station')
  z0 <- matrix(unlist(Z0),length(probs),length(dse.station))
  z <- matrix(unlist(Z),length(probs),length(dse.station))
  data(etopo5)
  etopo5 <- subset(etopo5,is=list(lon=range(lon),lat=range(lat)))
  etopo5[etopo5<=0] <- NA

  ## Set the grid to be the same as that of etopo5:
  grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')

  #3 Flag dubplicated stations:
  ok <- !(duplicated(lon) & duplicated(lat))

  ## Spread in the  90-percente interval changing
  obj <- switch( what,
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
  attr(W,'variable') <- varid(X)[1]
  attr(W,'unit') <- unit(X)[1]
  attr(W,'longitude') <- w$x
  attr(W,'latitude') <- w$y
  class(W) <- class(etopo5)

  ## Make a projection that zooms in on the Barents region
  dev.new()

  rev <- switch(param,'t2m'=FALSE,'pr'=TRUE)
  Wx <- max(abs(W),na.rm=TRUE)
  if (is.null(breaks)) breaks <- round(seq(-Wx,Wx,length=31),2) 
  map(W,xlim=range(lon(W)),ylim=range(lat(W)),projection='sphere',
      colbar=list(breaks=breaks,rev=rev,pal=pal))


  figlab(it)
  dev.copy2eps(file=paste('dse_paper59_Fig_ESD_',
                           param,'_',FUN,'_',it,'map.eps',sep=''))
  
}


#-----------------------------------------------------------------------
# Define season and parameter

param <- 't2m'
FUN <- 'mean'
it <- 'jja'
reanalysis <- 'air.mon.mean.nc'

if (param=='precip') FUNX <- 'C.C.eq' else FUNX <- 'mean'

## Get the predictand -> Y
if (!exists('Y')) {
  load(paste(param,'.aaca.rda',sep=''))
  Y <- eval(parse(text=paste(param,'.aaca',sep='')))
}

## Get the large-scale predictor:
if (!exists('predictor')) {
  T2M <- retrieve(reanalysis,lon=c(0,100),lat=c(60,90))
  predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it)
} else if (length(month(subset(predictor,it=it)))==0)
  predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it) 

## Carry out the downscaling:
if (!file.exists(paste('dse.aaca.t2m.rcp45.',it,'.rda',sep=''))) {
  dse.aaca.t2m.rcp45 <- downscale(Y,predictor,it,param,FUN=FUN,FUNX=FUNX)
  save(file=paste('dse.aaca.t2m.rcp45.',it,'.rda',sep=''),dse.aaca.t2m.rcp45)
}

if (!file.exists(paste('dse.aaca.t2m.rcp26.',it,'.rda',sep=''))) {
  dse.aaca.t2m.rcp26 <- downscale(Y,predictor,it,param,rcp='rcp26',FUN=FUN,FUNX=FUNX)
  save(file=paste('dse.aaca.t2m.rcp26.',it,'.rda',sep=''),dse.aaca.t2m.rcp26)
}

if (!file.exists(paste('dse.aaca.t2m.rcp85.',it,'.rda',sep=''))) {
  dse.aaca.t2m.rcp85 <- downscale(Y,predictor,it,param,rcp='rcp85',FUN=FUN,FUNX=FUNX)
  save(file=paste('dse.aaca.t2m.rcp85.',it,'.rda',sep=''),dse.aaca.t2m.rcp85)
}

print('--- Completed downscaling ---')

## Use LatticeKrig to grid the results - eg extract the 5 and 95 percentiles
## of projected temperature for 2050 from the ensembles.

breaks <- round(seq(0,12,length=25),2)

if (TRUE) {
gridmap('dse.aaca.t2m.rcp45.djf.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.djf.2050.eps')


gridmap('dse.aaca.t2m.rcp45.djf.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.djf.2030.eps')

gridmap('dse.aaca.t2m.rcp45.djf.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.djf.2100.eps')

gridmap('dse.aaca.t2m.rcp45.jja.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.jja.2050.eps')

gridmap('dse.aaca.t2m.rcp45.jja.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.jja.2030.eps')

gridmap('dse.aaca.t2m.rcp45.jja.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.jja.2100.eps')

gridmap('dse.aaca.t2m.rcp45.mam.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.mam.2050.eps')

gridmap('dse.aaca.t2m.rcp45.mam.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.mam.2030.eps')

gridmap('dse.aaca.t2m.rcp45.mam.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.mam.2100.eps')

gridmap('dse.aaca.t2m.rcp45.son.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.son.2050.eps')

gridmap('dse.aaca.t2m.rcp45.son.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.son.2030.eps')

gridmap('dse.aaca.t2m.rcp45.son.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp45.son.2100.eps')
}

if (TRUE) {
gridmap('dse.aaca.t2m.rcp85.djf.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.djf.2050.eps')

gridmap('dse.aaca.t2m.rcp85.djf.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.djf.2030.eps')

gridmap('dse.aaca.t2m.rcp85.djf.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.djf.2100.eps')

gridmap('dse.aaca.t2m.rcp85.jja.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.jja.2050.eps')

gridmap('dse.aaca.t2m.rcp85.jja.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.jja.2030.eps')

gridmap('dse.aaca.t2m.rcp85.jja.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.jja.2100.eps')

gridmap('dse.aaca.t2m.rcp85.mam.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.mam.2050.eps')

gridmap('dse.aaca.t2m.rcp85.mam.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.mam.2030.eps')

gridmap('dse.aaca.t2m.rcp85.mam.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.mam.2100.eps')

gridmap('dse.aaca.t2m.rcp85.son.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.son.2050.eps')

gridmap('dse.aaca.t2m.rcp85.son.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.son.2030.eps')

gridmap('dse.aaca.t2m.rcp85.son.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp85.son.2100.eps')
}

if (TRUE) {
gridmap('dse.aaca.t2m.rcp26.djf.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.djf.2050.eps')

gridmap('dse.aaca.t2m.rcp26.djf.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.djf.2030.eps')

gridmap('dse.aaca.t2m.rcp26.djf.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.djf.2100.eps')

gridmap('dse.aaca.t2m.rcp26.jja.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.jja.2050.eps')

gridmap('dse.aaca.t2m.rcp26.jja.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.jja.2030.eps')

gridmap('dse.aaca.t2m.rcp26.jja.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (JJA)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.jja.2100.eps')

gridmap('dse.aaca.t2m.rcp26.mam.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.mam.2050.eps')

gridmap('dse.aaca.t2m.rcp26.mam.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.mam.2030.eps')

gridmap('dse.aaca.t2m.rcp26.mam.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (MAM)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.mam.2100.eps')

gridmap('dse.aaca.t2m.rcp26.son.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.son.2050.eps')

gridmap('dse.aaca.t2m.rcp26.son.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.son.2030.eps')

gridmap('dse.aaca.t2m.rcp26.son.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='warm')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (SON)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.t2m.rcp26.son.2100.eps')
}


breaks <- round(seq(0,5,length=21),2)
if (TRUE) {
gridmap('dse.aaca.pr.rcp26.wetmean.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='precip')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.pr.rcp26.mu.2050.eps')

gridmap('dse.aaca.pr.rcp26.wetmean.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='precip')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.pr.rcp26.mu.2030.eps')

gridmap('dse.aaca.pr.rcp26.wetmean.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='precip')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP2.6 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.pr.rcp26.mu.2100.eps')

gridmap('dse.aaca.pr.rcp45.wetmean.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='precip')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.pr.rcp45.mu.2050.eps')

gridmap('dse.aaca.pr.rcp45.wetmean.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='precip')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.pr.rcp45.mu.2030.eps')

gridmap('dse.aaca.pr.rcp45.wetmean.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='precip')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP4.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.pr.rcp45.mu.2100.eps')

gridmap('dse.aaca.pr.rcp85.wetmean.rda',it=2050,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='precip')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2050 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.pr.rcp85.mu.2050.eps')

gridmap('dse.aaca.pr.rcp85.wetmean.rda',it=2030,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='precip')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2030 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.pr.rcp85.mu.2030.eps')

gridmap('dse.aaca.pr.rcp85.wetmean.rda',it=2100,it0=2010,what='change',probs=c(0.95),
        breaks=breaks,pal='precip')
figlab(expression(T[2*m]*(degree*C)),xpos = 0.001, ypos = 0.1)
figlab('95 percentile diff (2100 - 2010)',xpos = 0.001, ypos = 0.99)
figlab('RCP8.5 (DJF)',xpos = 0.8, ypos = 0.99)
dev.copy2eps(file='paper59_dse.aaca.pr.rcp85.mu.2100.eps')
}
