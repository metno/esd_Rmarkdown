#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

## R-job for downscaling Nordic temperatures for KSS 
print('--- R-job for downscaling Nordic temperatures for KSS ---')
print(args)

library(esd)

## Function definitions:

downscale <- function(Y,predictor,it='djf',param='t2m',FUN='mean',FUNX='mean',
                      period=c(1950,2015),plot=FALSE,rcp='rcp45',verbose=FALSE,
                      lon=c(-20,40),lat=c(50,80),ip=1:6,n=6,
                      pattern='tas_Amon_ens',
                      rel.cord=FALSE,select=NULL) {
  
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
  
  print(paste(round(100*sum(!is.finite(Y4))/length(!is.finite(Y4))),'% missing',
              sep=''))
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
                        biascorrect=TRUE,rcp=rcp,ip=ip,select=select,
                        pattern=pattern,
                        lon=lon,lat=lat,rel.cord=rel.cord,it=it,plot=plot)
  
  attr(dse.pca,'N.missing') <- nmiss
  invisible(dse.pca)
}

##-----------------------------------------------------------------------
## Define season and parameter: passed from the parameters

print('Settings')
param <- args[1]  # parameter
rcp <- args[2]
it <- args[3]
if (param=='t2m') {
  reanalysis <- 'air.mon.mean.nc'
  FUN='mean'
  FUNX='mean'
  pattern='tas_Amon_ens'
  predictand <- '/lustre/storeB/users/rasmusb/data/t2m.nordic.rda'
  varid <- 'tg'
}
if (param=='mu') {
  reanalysis <- 'air.mon.mean.nc'
  FUN='wetmean'
  FUNX='meanx'
  pattern='tas_Amon_ens'
  predictand <- '/lustre/storeB/users/rasmusb/data/rr.nordic.rda'
  varid <- 'precip'
}
if (param=='fw') {
  reanalysis <- 'slp.mon.mean.nc'
  FUN='wetfreq'
  FUNX='mean'
  pattern='psl_Amon_ens'
  predictand <- '/lustre/storeB/users/rasmusb/data/rr.nordic.rda'
  varid <- 'precip'
}

outdir <- '/lustre/storeB/users/rasmusb/dse4KSS/'
if (!file.exists(outdir)) dir.create(outdir)

lon=c(-30,30); lat=c(50,75)

verbose=FALSE

## Retrieving the data
## Get the predictand data: daily temperature from ECAD:

print('Retrieve data')

if (!file.exists(predictand)) {
  ss <- select.station(param=varid,lon=c(5,30),lat=c(55,72),src='ecad',nmin=50)
  Y <- station(ss)
  save(Y,file=predictand)
} else load(predictand)

## Only use the most recent 66 years 
Y <- subset(Y,it=c(1960,2015))
## Remove records with many missing datapoints
nv <- apply(coredata(Y),2,'nv')
Y <- subset(Y,is=nv > 14000)

print('predictand')

  ## Get the large-scale predictor:
  if (!exists('predictor')) {
    T2M <- retrieve(reanalysis,lon=lon,lat=lat)
    predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it)
  } else if (length(month(subset(predictor,it=it)))==0)
    predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it) 
  
  print(paste('Generating dse.kss.t2m.rcp45.',it,'.rda',sep=''))
  
  ## Carry out the downscaling:
  if (!file.exists(paste('dse.kss.t2m.',rcp,',.',it,'.rda',sep=''))) {
    Z <- downscale(Y,predictor,it,param,rcp='rcp26',FUN=FUN,FUNX=FUNX,
                   lon=lon,lat=lat,verbose=verbose)
    save(file=paste(outdir,'dse.kss.t2m.',rcp,'.',it,'.rda',sep=''),Z)
  }
  
 
print('--- Completed downscaling ---')

