## Script to carry out sensitivity analysis:
## The effect of the mean seasonal variations
## Geographical variations
## Long-term trends

library(esd)

gridmap <- function(X,projection="sphere",do.pca=TRUE,
                            xlim=NULL,ylim=NULL,new=TRUE,
                            FUN="mean",pattern=1,save.pdf=FALSE,
                            verbose=FALSE,...) {
  if (verbose) print("gridmap")
  require(LatticeKrig)
  data(etopo5)
  etopo5 <- subset(etopo5,is=list(lon=range(lon(X)),lat=range(lat(X))))
  etopo5[etopo5<=0] <- NA
  if (do.pca) {
    if (verbose) print("Principle component analysis")
    pca <- PCA(X)
    z <- attr(pca,'pattern')[,pattern]
  } else {
    if (verbose) print(paste("Calculate",FUN))
    z <- apply(coredata(X),2,FUN,na.rm=TRUE)
  }
  if (verbose) print(paste("Grid the data with LatticeKrig"))
  grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')
  ok <- !(duplicated(lon(X)) & duplicated(lat(X)))
  obj <- LatticeKrig( x=cbind(lon(X)[ok],lat(X)[ok]),
                      y=z[ok],Z=alt(X)[ok])
  w <- predictSurface(obj, grid.list = grid,Z=etopo5)
  w$z[is.na(etopo5)] <- NA

  if (verbose) print("Convert from LatticeKrig to esd")
  W <- w$z
  param <- varid(X)[1]
  units <- attr(X,"unit")[1]
  if (do.pca) {
    attr(W,'variable') <- paste("PC",pattern," ",param,sep="")
    attr(W,'unit') <- NULL
  } else if (FUN=="trend.coef") {
    attr(W,'variable') <- paste(param,"trend")
    attr(W,'unit') <- paste(units," decade^{-1}")
  } else {
    attr(W,'variable') <- paste(FUN,param)
    attr(W,'unit') <- units
  }
  attr(W,'longitude') <- w$x
  attr(W,'latitude') <- w$y
  class(W) <- class(etopo5)
  if(verbose) print("Get rid of packages that have functions of same name")
  detach("package:LatticeKrig")
  detach("package:fields")
  detach("package:spam")
  detach("package:grid")
  detach("package:maps")
  if (verbose) print("Project the results on a map")
  if (save.pdf) {
    if (do.pca) {
      fname <- paste("map.",varid(X)[1],".PC",pattern,".pdf",sep="")
    } else {
      fname <- paste("map",varid(X)[1],FUN,"pdf",sep=".")
    }
    fname <- gsub("[^\\.[:^punct:]]","",fname,perl=T)
    pdf(fname,width=8,height=8)
  } else if (new) {
    dev.new()
  }
  map(W,projection=projection,xlim=xlim,ylim=ylim,new=FALSE,
      verbose=verbose,...)
  if (save.pdf) dev.off()
  invisible(W)
}

#param <- 't2m'
#it <- 'annual'
#FUN <- 'mean'
param <- 'pr'
it <- 'annual'
FUN <- 'sum'
xlim <- c(0,120)
ylim <- c(62,85)
save.pdf <- TRUE
colbar.trend <- list(pal="t2m",breaks=seq(-1,1,0.1))

load(paste(param,'.aaca.rda',sep=''))
X <- eval(parse(text=paste(param,'.aaca',sep='')))
Xc <- aggregate(X,month,FUN)

X <- subset(X,it=c(1960,2015))
nok <- apply(coredata(X),2,nv)
X <- subset(X,is=(nok > 15000))
if (param=='pr') nmin=300 else nmin=365
if (sum(is.element(season.abb(),it))>0) {
  X <- subset(as.4seasons(X,FUN=FUN),it=it)
  index(X) <- year(X)
} else
if (it=='annual') X <- annual(X,FUN=FUN,nmin=nmin) else
if (it=='seasonal') X <- aggregate(X,month,FUN=FUN)
## Fill in missing data
print('missing values: pcafill')
if (param=='t2m') X <- pcafill(X)
if ((param=='pr') & (FUN=='sum') & (it=='annual'))
  attr(X,'unit')[] <- 'mm/year'

attr(Xc,"variable") <- gsub("[^\\.[:^punct:]]","",varid(X)[1],perl=T)
attr(X,"variable") <- gsub("[^\\.[:^punct:]]","",varid(X)[1],perl=T)

## Make gridded maps of a, b, c) PC:s, d) mean and e) trend
## !!!ADD point observations in gridmap!!!

## Mean climatology
print('Grid the mean clmatology')
gridmap(X,do.pca=FALSE,xlim=xlim,ylim=ylim,new=FALSE)
stop()

gridmap(X,do.pca=TRUE,xlim=xlim,ylim=ylim,save.pdf=save.pdf,
        pattern=1,verbose=TRUE)
gridmap(X,do.pca=TRUE,xlim=xlim,ylim=ylim,save.pdf=save.pdf,
        pattern=2,verbose=TRUE)
gridmap(X,do.pca=TRUE,xlim=xlim,ylim=ylim,save.pdf=save.pdf,
        pattern=3,verbose=TRUE)
gridmap(X,do.pca=FALSE,FUN="mean",xlim=xlim,ylim=ylim,
        save.pdf=save.pdf)
gridmap(X,do.pca=FALSE,FUN="trend.coef",xlim=xlim,ylim=ylim,
        colbar=colbar.trend,save.pdf=save.pdf)

## Principle component analysis
pca <- PCA(X)
plot(pca)
dev.copy2pdf(file=paste("pca",varid(X),"pdf",sep="."))
dev.off()

## Seasonal cycle
plot(Xc)
dev.copy2pdf(file=paste("seasonalcycle",varid(Xc),"pdf",sep="."))
dev.off()

## pca of seasonal cycle?
uwv <- svd(coredata(Xc))
