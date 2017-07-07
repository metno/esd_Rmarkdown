#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

## R-job for downscaling Nordic temperatures for KSS 
print('--- R-job for downscaling Nordic temperatures for KSS ---')
print(args)

#library(devtools)
#install.packages('esd',repos=NULL)
library(esd)

## Function definitions:

downscale <- function(Y,predictor,it='djf',param='t2m',FUN='mean',FUNX='mean',
                      period=c(1950,2015),plot=FALSE,rcp='rcp45',verbose=FALSE,
                      lon=c(-20,40),lat=c(50,80),ip=1:7,n=5,
                      pattern='tas_Amon_ens',
                      rel.cord=FALSE,select=NULL) {
  
  if (verbose) print(paste('downscale',param,it,FUN,FUNX))
  
  ## Use a time and space window:  
  if (verbose) print(paste('subset',paste(period,collapse='-')))
  Y <- subset(Y,it=period)
  
  ## Estimate seasonal means & weed out stations with little data
  if (verbose) print('season')
  Y4 <- subset(as.4seasons(Y,FUN=FUN),it=it)
  if (verbose) print('weed out missing data')
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
  if (verbose) print('pcafill')
  Z <- pcafill(Y4)
  ## Negative precipitation is impossible - clip to zero
  if (is.precip(Z)) coredata(Z)[Z<0]<- 0
  
  if (verbose) print('pca')
  pca <- PCA(Z,n=n)
  if (plot) plot(pca)
  
  ## Downscale results
  print(paste('DSensemble',varid(predictor)))
  dse.pca <- DSensemble(pca,predictor=predictor,FUNX=FUNX,verbose=verbose,
                        biascorrect=TRUE,rcp=rcp,ip=ip,select=select,
                        nmin=1,pattern=pattern,lon=lon,lat=lat,rel.cord=rel.cord,it=it,plot=plot)
  attr(dse.pca,'N.missing') <- nmiss
  invisible(dse.pca)
}

##-----------------------------------------------------------------------
## Define season and parameter: passed from the parameters

print('Settings')
param <- args[1]  # parameter
rcp <- args[2]
it <- args[3]
landmask <- FALSE
ip <- 1:5

sessionInfo()

if (param=='t2m') {
  reanalysis <- 'air.mon.mean.nc'
  #reanalysis <- 'ERA40_t2m_mon.nc'
  FUN='mean'
  FUNX='mean'
  pattern='tas_Amon_ens'
  predictand <- '/lustre/storeB/users/rasmusb/data/t2m.nordic.rda'
  varid <- 'tg'
  lon=c(-20,40); lat=c(55,80)
  n=3
  itp <- c(1960,2015)
}
if (param=='mu') {
  reanalysis <- 'air.mon.mean.nc'
  #reanalysis <- 'ERA40_t2m_mon.nc'
  FUN='wetmean'
  FUNX='C.C.eq'
  #FUNX <- 'mean'
  pattern='tas_Amon_ens'
  predictand <- '/lustre/storeB/users/rasmusb/data/rr.nordic.rda'
  varid <- 'precip'
  lon=c(-90,10); lat=c(10,50)
  landmask <- TRUE
  if (it=='mam') {lon=c(-50,20); lat=c(40,65)}
  if (it=='jja') {lon=c(-20,30); lat=c(50,70); landmask <- FALSE}
  if (it=='son') {lon=c(-50,20); lat=c(40,65)}
  n=3
  itp <-c(1965,2015)
}
if (param=='fw') {
  reanalysis <- 'slp.mon.mean.nc'
  #reanalysis <- 'ERA40_slp_mon.nc'
  FUN='wetfreq'
  FUNX='mean'
  pattern='psl_Amon_ens'
  predictand <- '/lustre/storeB/users/rasmusb/data/rr.nordic.rda'
  varid <- 'precip'
  #lon=c(-30,30); lat=c(50,75)
  lon=c(0,30); lat=c(50,70)
  n=3
  itp <- c(1965,2015)
}

outdir <- '/lustre/storeB/users/rasmusb/dse4KSS/'
if (!file.exists(outdir)) dir.create(outdir)

verbose=TRUE

## Retrieving the data
## Get the predictand data: daily temperature from ECAD:

print('Retrieve data')

if (!file.exists(predictand)) {
  ss <- select.station(param=varid,lon=c(5,30),lat=c(55,72),src='ecad',nmin=50)
  Y <- station(ss)
  save(Y,file=predictand)
} else load(predictand)

## Only use the most recent 66 years 
Y <- subset(Y,it=itp)
## Remove records with many missing datapoints
nv <- apply(coredata(Y),2,'nv')
Y <- subset(Y,is=nv > 16000)

print(paste('predictor',reanalysis,paste(lon,colapse='-'),'/',paste(lat,collapse='-'),landmask))

  ## Get the large-scale predictor:
if (!exists('predictor')) {
    T2M <- retrieve(reanalysis,lon=lon,lat=lat)
    predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it)
} else if (length(month(subset(predictor,it=it)))==0)
    predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it) 

if (landmask) predictor <- mask(predictor,land=TRUE)

print(paste('Generating dse.kss.',param,'.',rcp,'.',it,'.rda',sep=''))
  
  ## Carry out the downscaling:
if (!file.exists(paste('dse.kss.t2m.',rcp,',.',it,'.rda',sep=''))) {
    Z <- downscale(Y,predictor,it,param,rcp=rcp,FUN=FUN,FUNX=FUNX,
                   lon=lon,lat=lat,n=n,verbose=verbose)
    attr(Z,'predictor_file') <- 'predictor'
    attr(Z,'predictor_lon') <- range(lon)
    attr(Z,'predictor_lat') <- range(lat)
    attr(Z,'predictor_landmask') <- landmask
    save(file=paste(outdir,'dse.kss.',param,'.',rcp,'.',it,'.rda',sep=''),Z)
}
  
 
print('--- Completed downscaling ---')

