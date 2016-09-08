## Downscaling exercise to examine the link between large-scale predictors and the local
## precipitation parameters.

library(esd)

param <- 'pr'
FUN <- 'wetfreq'
#it <- NULL
#it <- month.abb[4:9]
it <- month.abb[c(1:3,10:12)]
period <- c(1950,2015)
eofs <- 1:4
n <- 10
verbose <- FALSE
plot <- FALSE

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

## The predictor:

  predictor <- retrieve(reanalysis,lon=lon,lat=lat)
  predictor <- subset(predictor,it=it)
  if (FUNX=='C.C.eq') predictor <- C.C.eq(predictor)
  eof <- EOF(annual(predictor,nmin=6))

## Get the predictand -> Y

  load(paste(param,'.aaca.rda',sep=''))
  Y <- eval(parse(text=paste(param,'.aaca',sep='')))
  Y <- subset(Y,it=period)

  ## Extract summer or winter
  if (!is.null(it)) {
    print(paste('Extract',paste(it,collapse='-')))
    Y <- subset(Y,it=it)
    nmin <- 28*length(it)
  } else nmin=360
  if (verbose) {print(paste('nmin=',nmin)); print(dim(Y))}
  ## Remove suspect data:
  cY <- coredata(Y); cY[cY > 500] <- NA; cY -> coredata(Y)


## Estimate seasonal means & weed out stations with little data
  Y4 <- annual(Y,FUN=FUN,nmin=6*28)
  counts <- zoo(annual(Y,FUN='count',nmin=6*28))
  ok <- apply(coredata(Y4),1,nv)
  Y4 <- subset(Y4,it=ok>0)
  nok <- apply(coredata(Y4),2,nv)
  Y4 <- subset(Y4,is=nok>15)
  if (verbose) {print(ok); print(nok); str(Y4)}

  if (plot) map.station(Y4,FUN='mean',cex=-2)

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

  ds <- DS(pca,eof)
  plot(ds)
  if (is.null(it)) s <- '' else
           s <- tolower(paste('.',paste(substr(it,1,1),collapse=''),sep=''))
  dev.copy2pdf(file=paste('paper59_ds.',FUN,'',s,'.pdf',sep=''))
