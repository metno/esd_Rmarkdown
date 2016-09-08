### Downscale scenarios for AACA - the Barents region
### temperature and precipitation.
### Based on PCA-representation for the predictand.
### Use 'pcafill' to fill in missing values.
### Rasmus E. Benestad, Oslo, November, 2015

library(esd)


downscale <- function(Y,predictor,it='djf',param='t2m',FUN='mean',FUNX='mean',
                      period=c(1950,2015),plot=FALSE,rcp='rcp45',verbose=FALSE,
                      lon=c(0,100),lat=c(65,90),eofs=1:6,n=10,rel.cord=FALSE,select=NULL) {

  print('downscale')
  
## Use a time and space window:  
  Y <- subset(Y,it=period)
  Y <- subset(Y,is=list(lat=lat,lon=lon))

## Estimate seasonal means & weed out stations with little data
  Y4 <- subset(as.4seasons(Y,FUN=FUN),it=it)
  ok <- apply(coredata(Y4),1,nv)
  Y4 <- subset(Y4,it=ok>0)
  nok <- apply(coredata(Y4),2,nv)
  Y4 <- subset(Y4,is=nok>15)

  print(paste(round(100*sum(!is.finite(Y4))/length(!is.finite(Y4))),'% missing',sep=''))
  if (plot) map(Y,FUN=FUN,cex=-2)

  nmiss <- round(100*sum(!is.finite(Y4))/length(!is.finite(Y4)))
  #print(paste(nmiss,'% missing',sep=''))
  
## Fill missing data using PCA-based regression
  Z <- pcafill(Y4)

  pca <- PCA(Z,n=n)
  if (plot) plot(pca)

## Downscale results
  print('DSensemble')
  dse.pca <- DSensemble(pca,predictor=predictor,FUNX=FUNX,verbose=verbose,
                        biascorrect=TRUE,rcp=rcp,eofs=eofs,select=select,
                        lon=lon,lat=lat,rel.cord=rel.cord,it=it)

  attr(dse.pca,'N.missing') <- nmiss
  invisible(dse.pca)
}



#-----------------------------------------------------------------------
# Define season and parameter

param <- 't2m'
FUN <- 'mean'
reanalysis <- 'air.mon.mean.nc'
lon=c(0,100); lat=c(65,90)

if (param=='precip') FUNX <- 'C.C.eq' else FUNX <- 'mean'

## Get the predictand -> Y
if (!exists('Y')) {
  load(paste(param,'.aaca.rda',sep=''))
  Y <- eval(parse(text=paste(param,'.aaca',sep='')))
}

for (it in c('djf','mam','jja','son')) {

  ## Get the large-scale predictor:
  if (!exists('predictor')) {
    T2M <- retrieve(reanalysis,lon=lon,lat=lat)
    predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it)
  } else if (length(month(subset(predictor,it=it)))==0)
    predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it) 

  print(paste('Generating dse.aaca.t2m.rcp45.',it,'.rda',sep=''))
  
  ## Carry out the downscaling:
  if (!file.exists(paste('dse.aaca.t2m.rcp45.',it,'.rda',sep=''))) {
    dse.aaca.t2m.rcp45 <- downscale(Y,predictor,it,param,FUN=FUN,FUNX=FUNX,lon=lon,lat=lat)
    save(file=paste('dse.aaca.t2m.rcp45.',it,'.rda',sep=''),dse.aaca.t2m.rcp45)
  }

  print(paste('Generating dse.aaca.t2m.rcp26.',it,'.rda',sep=''))
  
  if (!file.exists(paste('dse.aaca.t2m.rcp26.',it,'.rda',sep=''))) {
    dse.aaca.t2m.rcp26 <- downscale(Y,predictor,it,param,rcp='rcp26',
                                    FUN=FUN,FUNX=FUNX,lon=lon,lat=lat)
    save(file=paste('dse.aaca.t2m.rcp26.',it,'.rda',sep=''),dse.aaca.t2m.rcp26)
  }

   print(paste('Generating dse.aaca.t2m.rcp85.',it,'.rda',sep=''))
  
  if (!file.exists(paste('dse.aaca.t2m.rcp85.',it,'.rda',sep=''))) {
    dse.aaca.t2m.rcp85 <- downscale(Y,predictor,it,param,rcp='rcp85',
                                    FUN=FUN,FUNX=FUNX,lon=lon,lat=lat)
    save(file=paste('dse.aaca.t2m.rcp85.',it,'.rda',sep=''),dse.aaca.t2m.rcp85)
  }
}
print('--- Completed downscaling ---')

