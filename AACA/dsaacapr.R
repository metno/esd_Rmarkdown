### Downscale scenarios for AACA - the Barenst region
### temperature and precipitation.
### Based on PCA-representation for the predictand.
### Use 'pcafill' to fill in missing values.

library(esd)
#library(LatticeKrig)

distill <- function(x,it=2050,probs=c(0.05,0.95)) {
  y <- subset(x,it=it)
  attributes(y) <- NULL
  y <- quantile(y,probs=probs,na.rm=TRUE)
  return(y)
}

downscale <- function(Y,predictor,param='precip',it=NULL,FUN='wetmean',FUNX='C.C.eq',
                      period=c(1950,2015),plot=FALSE,rcp='rcp45',
                      pattern="tas_Amon_ens_",
                      verbose=FALSE,eofs=1:6,n=10,prop=FALSE,
                      lon=c(-30,20),lat=c(-20,10),rel.cord=FALSE) {

  if (verbose) {print('downscale'); if (!is.null(it)) print(it)}
## Use a time and space window:  
  Y <- subset(Y,it=period)
  #Y <- subset(Y,is=list(lat=lat,lon=lon))
  ## Extract summer or winter
  if (!is.null(it)) {
    print(paste('Extract',paste(it,collapse='-')))
    Y <- subset(Y,it=it)
    nmin <- 28*length(it)
  } else nmin=360
  if (verbose) {print(paste('nmin=',nmin)); print(dim(Y))}

## Estimate seasonal means & weed out stations with little data
  Y4 <- annual(Y,FUN=FUN,nmin=nmin)
  counts <- zoo(annual(Y,FUN='count',nmin=nmin))
  ok <- apply(coredata(Y4),1,nv)
  Y4 <- subset(Y4,it=ok>0)
  nok <- apply(coredata(Y4),2,nv)
  Y4 <- subset(Y4,is=nok>15)
  if (verbose) {print(ok); print(nok); str(Y4)}

  if (plot) map(Y4,FUN=FUN,cex=-2)

  ## For wet-mean 
  if (prop) {
    X <- coredata(Y4)
    X <- t(100*t(X)/apply(X,2,'mean',na.rm=TRUE))
    X -> coredata(Y4)
    attr(Y4,'unit') <- '%'
  }

  nmiss <- round(100*sum(!is.finite(Y4))/length(!is.finite(Y4)))
  print(paste(nmiss,'% missing',sep=''))

## Fill missing data using PCA-based regression
  Z <- pcafill(Y4)

  pca <- PCA(Z,n=n)
  if (plot) plot(pca)

## Downscale results
  dse.pca <- DSensemble(pca,predictor=predictor,
                        biascorrect=TRUE,rcp=rcp,FUNX=FUNX,pattern=pattern,
                        verbose=verbose,eofs=eofs,lon=lon,lat=lat,rel.cord=rel.cord)
  attr(dse.pca,'N.missing') <- nmiss
  attr(dse.pca,'precipitation_counts') <- counts
  invisible(dse.pca)
}

# Define season and parameter
#detach('package:LatticeKrig')
#detach('package:fields')
#detach('package:spam')
#detach('package:grid')

param <- 'pr'
eofs <- 1:4

for (FUN in c('wetmean','wetfreq')) {
for (i in c(3,2)) {
if (i==1) it <- NULL
if (i==2) it <- month.abb[4:9]
if (i==3) it <- month.abb[c(1:3,10:12)]

#FUN <- 'wetfreq'
if (FUN=='wetmean') {
  FUNX <- 'C.C.eq'
  pattern <- "tas_Amon_ens_"
  reanalysis <- 'air.mon.mean.nc'
  prop=TRUE
  lon <- c(-50,10); lat <- c(40,75)
} else {
  FUNX <- 'mean'
  pattern <- "psl_Amon_ens_"
  reanalysis <- 'slp.mon.mean.nc'
  prop=FALSE
  lon <- c(-40,70); lat <- c(50,85)
}

## Get the predictand -> Y
if (!exists('Y')) {
  load(paste(param,'.aaca.rda',sep=''))
  Y <- eval(parse(text=paste(param,'.aaca',sep='')))
}

## Get the large-scale predictor:
if (!exists('predictor')) {
  T2M <- retrieve(reanalysis,lon=c(0,100),lat=c(60,90))
  if (reanalysis=='air.mon.mean.nc') attr(T2M,'units') <- 'degC'
  if (FUNX=='C.C.eq') T2M <- C.C.eq(T2M)
  predictor <- annual(T2M)
} 

## Carry out the downscaling:
if (is.null(is)) ses <- '' else
                 ses <- tolower(paste('.',paste(paste(substr(it,1,1),collapse=''),'.',sep=''),sep=''))

if (!file.exists(paste('dse.aaca.pr.rcp45',ses,FUN,'.rda',sep=''))) {
  dse.aaca.pr.rcp45 <- downscale(Y,predictor,param,it=it,FUN=FUN,FUNX=FUNX,prop=prop,
                                  pattern=pattern,lon=lon,lat=lat,eofs=eofs)
  save(file=paste('dse.aaca.pr.rcp45',ses,FUN,'.rda',sep=''),dse.aaca.pr.rcp45)
}

if (!file.exists(paste('dse.aaca.pr.rcp26',ses,FUN,'.rda',sep=''))) {
  dse.aaca.pr.rcp26 <- downscale(Y,predictor,param,it=it,rcp='rcp26',FUN=FUN,prop=prop,
                                  FUNX=FUNX,pattern=pattern,lon=lon,lat=lat,eofs=eofs)
  save(file=paste('dse.aaca.pr.rcp26',ses,FUN,'.rda',sep=''),dse.aaca.pr.rcp26)
}

if (!file.exists(paste('dse.aaca.pr.rcp85',ses,FUN,'.rda',sep=''))) {
  dse.aaca.pr.rcp85 <- downscale(Y,predictor,param,it=it,rcp='rcp85',FUN=FUN,prop=prop,
                                  FUNX=FUNX,pattern=pattern,lon=lon,lat=lat,eofs=eofs)
  save(file=paste('dse.aaca.pr.rcp85',ses,FUN,'.rda',sep=''),dse.aaca.pr.rcp85)
}

}
}

source('dsaaca.R')

print('--- Completed downscaling ---')

source('paper59_Atlas.R')
