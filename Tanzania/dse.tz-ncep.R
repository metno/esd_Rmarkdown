
## Script to downscale temperature from Tanzania.
## Data from Habiba.

library(esd)
setwd('~/R/storeA')

predictands <- function(param='tmin',incl.ghcn=FALSE) {
  ## Prepare the predictand data
  #print("predictands")
  fname <- paste(param,'.tz.rda',sep='')
  path <- paste('Habiba/',param,sep='')
  if (!file.exists(fname)) {
  files <- list.files(path,pattern='rda',full.names=TRUE)
  locs <- list.files(path,pattern='rda')
  dse <- grep('dse.tm',files)
  files <- files[-dse]; locs <- locs[-dse]
  print(files)
  locs <- substr(locs,1,nchar(locs)-4)
  Y <- NULL
  for (i in 1:length(files)) {
    load(files[i])
    y1 <- eval(parse(text=locs[i]))
    if (is.null(Y)) Y <- y1 else
                    Y <- combine.stations(Y,y1)
  }
  save(file=fname,Y)
  plot(Y)
} else load(fname)
  #print('GHCN')
  if (incl.ghcn) {
    xstations <- switch(param,'tmin'='Habiba/tn.eafrica.rda',
                            'tmax'='Habiba/tx.eafrica.rda')
    load(xstations)
    if (param=='tmin') y <- tn.eafrica else y <- tx.eafrica
    y <- subset(y,it=c(1950,2015))
    Y <- combine(Y,y)
  }
  map(Y,FUN='mean',cex=-2)
  invisible(Y)
}

downscale <- function(Y,predictor,it='djf',param='t2m',FUN='mean',FUNX='mean',
                      period=c(1950,2015),plot=FALSE,rcp='rcp45',verbose=FALSE,
                      lon=c(0,100),lat=c(65,90),eofs=1:5,n=3,rel.cord=FALSE,select=1:5) {

  print('downscale')
  
## Use a time and space window:  
  Y <- subset(Y,it=period)
  Y <- subset(Y,is=list(lat=lat,lon=lon))

## Estimate seasonal means & weed out stations with little data
  Y4 <- subset(as.4seasons(Y,FUN=FUN,nmin=60),it=it)
  ok <- apply(coredata(Y4),1,nv)
  Y4 <- subset(Y4,it=ok>0)
  nok <- apply(coredata(Y4),2,nv)
  Y4 <- subset(Y4,is=nok>15)

  print(paste(round(100*sum(!is.finite(Y4))/length(!is.finite(Y4))),'% missing',sep=''))
  if (plot) map(Y,FUN=FUN,cex=-2)

  nmiss <- round(100*sum(!is.finite(Y4))/length(!is.finite(Y4)))
  print(paste(nmiss,'% missing',sep=''))
  
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

param <- 'tmax'
FUN <- 'mean'
reanalysis <- 'air.mon.mean.nc'
#reanalysis <- 'ERAINT_t2m_mon.nc'
FUNX <- 'mean'
data.path <- 'dse.Tz.NCEP/'

predictands(param) -> Y

Y <- subset(Y,it=c(1961,2011))
lon=range(lon(Y))+c(-5,5); lat=range(lat(Y))+c(-5,5)


## Get the predictand -> Y

for (it in c('djf','mam','jja','son')) {

  ## Get the large-scale predictor:
  if (!exists('predictor')) {
    T2M <- retrieve(reanalysis,lon=lon,lat=lat)
    predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it)
  } else if (length(month(subset(predictor,it=it)))==0)
    predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it) 

  #print(paste('Generating dse.',param,'.tz.rcp45.',it,'.rda',sep=''))
  
  ## Carry out the downscaling:
  if (!file.exists(paste(data.path,'dse.',param,'.tz.rcp45.',it,'.rda',sep=''))) {
    #dev.new()
    dse.t2m.tz.rcp45 <- downscale(Y,predictor,it,param,FUN=FUN,FUNX=FUNX,lon=lon,lat=lat)
    save(file=paste(data.path,'dse.',param,'.tz.rcp45.',it,'.rda',sep=''),dse.t2m.tz.rcp45)
  }

  #print(paste('Generating dse.',param,'.tz.rcp26.',it,'.rda',sep=''))
  
  if (!file.exists(paste(data.path,'dse.',param,'.tz.rcp26.',it,'.rda',sep=''))) {
    dse.t2m.tz.rcp26 <- downscale(Y,predictor,it,param,rcp='rcp26',
                                    FUN=FUN,FUNX=FUNX,lon=lon,lat=lat)
    save(file=paste(data.path,'dse.',param,'.tz.rcp26.',it,'.rda',sep=''),dse.t2m.tz.rcp26)
  }

   #print(paste('Generating dse.',param,'.tz.rcp85.',it,'.rda',sep=''))
  
  if (!file.exists(paste(data.path,'dse.',param,'.tz.rcp85.',it,'.rda',sep=''))) {
    dse.t2m.tz.rcp85 <- downscale(Y,predictor,it,param,rcp='rcp85',
                                    FUN=FUN,FUNX=FUNX,lon=lon,lat=lat)
    save(file=paste(data.path,'dse.',param,'.tz.rcp85.',it,'.rda',sep=''),dse.t2m.tz.rcp85)
  }
}
print('--- Completed downscaling ---')
setwd('~/R')
