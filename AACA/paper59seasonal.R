## Script to carry out sensitivity analysis:
## The effect of the mean seasonal variations
## Geographical variations
## Long-term trends

rm(list=ls())
library(esd)

ar1 <- function(x) {
  n <- length(x)
  x0 <- x[1:(n-1)]
  x1 <- x[2:n]
  ok <- is.finite(x0) & is.finite(x1)
  ar1 <- cor(x0[ok],x1[ok])
  return(ar1)
}

analyseseason <- function(param='pr',FUN='wetmean',prop=FALSE,pattern=1,
                          map.show=FALSE,grid=FALSE) {
  
load(paste(param,'.aaca.rda',sep=''))
X <- eval(parse(text=paste(param,'.aaca',sep='')))

if (FUN=='sd') X <- anomaly(X)
if ( (FUN!='wet.spell') & (FUN!='dry.spell'))
  Xc <- aggregate(X,month,FUN) else {
    ispell <- switch(FUN,'wet.spell'=1,'dry.spell'=2)
    ## For spell-length statistics, more a special treatment is needed
    ## Remove records with many missing data (more than 2.5%).
    nok <- (apply(X,2,nv) > 0.975*dim(X)[1])
    X <- subset(X,is=nok)
    x <- list()
    for (i in 1:dim(X)[2]) x[[i]] <- subset(spell(subset(X,is=i),threshold=1),is=ispell)
    xx <- do.call('merge',x)
    class(xx) <- class(X)
    xx <- attrcp(X,xx)
    attr(xx,'variable') <- FUN
    attr(xx,'unit') <- 'days'
    Xc <- aggregate(xx,month,"mean")
}

if (prop) {
  coredata(Xc) <- t(t(coredata(Xc))/apply(coredata(subset(Xc,it=c(1961,1990))),2,'mean',na.rm=TRUE))
  FUN <- paste(FUN,'prop',sep='')
  attr(Xc,'unit') <- '%'
}

varnm <- switch(paste(param,FUN),
                't2m mean'=expression(T[2*m]),
                't2m sd'=expression(sigma[T]),
                't2m ar1'=expression(R(tau)),
                'pr wetmean'=expression(mu),
                'pr wetfreq'=expression(f[w]),
                'pr wet.spell'=expression(n[c*w*d]),
                'pr dry.spell'=expression(n[c*d*d]),
                'pr wetmeanprop'=expression(mu))
unitnm <- switch(paste(param,FUN),
                't2m mean'=expression(degree*C),
                't2m sd'=expression(degree*C),
                't2m ar1'='correlation',
                'pr wetmean'='mm/s',
                'pr wetfreq'='',
                'pr wet.spell'='days',
                'pr dry.spell'='days',
                'pr wetmeanprop'="'%'")
longnm <- switch(paste(param,FUN),
                't2m mean'='Mean temperature',
                't2m sd'='Temperature variability',
                't2m ar1'='Temperature persistence',
                'pr wetmean'='Precipitation intensity',
                'pr wetfreq'='Wet-day frequency',
                'pr wet.spell'='Mean wet-spell length',
                'pr dry.spell'='Mean dry-spell length',
                'pr wetmeanprop'='Proportional precipitation intensity')

nv <- apply(X,2,nv)
xlim <- range(lon(X)); ylim <- range(lat(X))
col <- rgb(0.75*(1-(lat(X)-ylim[1])/diff(ylim)),0.75*(lon(X)-xlim[1])/diff(xlim),
           0.75*(lat(X)-ylim[1])/diff(ylim),0.75*nv/dim(X)[1])
attr(Xc,'variable') <- varnm
attr(Xc,'longname') <- longnm
attr(Xc,'unit') <- unitnm
plot(Xc,col=col,map.show=map.show,map.insert=FALSE,cex.main=2,errorbar=FALSE)

if (map.show) {
  dev.copy2pdf(file=paste('paper59_Fig_seasonalvar_',
               param,'_',FUN,'.pdf',sep=''))
  dev.off()
  dev.copy2pdf(file=paste('paper59_Fig_seasonalvar_',
               param,'_map.pdf',sep=''))
  dev.off()
} else {
  dev.copy2pdf(file=paste('paper59_Fig_seasonalvar_',
               param,'_',FUN,'.pdf',sep=''))
  dev.off()
}

if (grid) {
pca <- PCA(pcafill(Xc),n=12)
attr(pca,'variable') <- varnm
attr(pca,'unit') <- unitnm
attr(pca,'longname') <- FUN

plot(pca,pattern=pattern)
dev.copy2pdf(file=paste('paper59_Fig_pcaseasonalvar_',pattern,'_',
               param,'_',FUN,'.pdf',sep=''))

z <- attr(pca,'pattern')[,pattern]*attr(pca,'eigenvalues')[pattern]

## Get elevation data
data(etopo5)
etopo5 <- subset(etopo5,is=list(lon=range(lon(X)),lat=range(lat(X))))
etopo5[etopo5<=0] <- NA

## Grid the PCA pattern
require(LatticeKrig)

## Set the grid to be the same as that of etopo5:
grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')

## Flag dubplicated stations:
ok <- !(duplicated(lon(X)) & duplicated(lat(X)))

## 
obj <- LatticeKrig( x=cbind(lon(X)[ok],lat(X)[ok]),
                    y=z[ok],Z=alt(X)[ok])

##  obj <- LatticeKrig( x=cbind(lon[ok],lat[ok]), y=z[2,ok],Z=alt[ok])
w <- predictSurface(obj, grid.list = grid,Z=etopo5)
w$z[is.na(etopo5)] <- NA
dev.new()
#surface( w )

## Get rid of packages that have functions of same name:
detach("package:LatticeKrig")
detach("package:fields")
detach("package:spam")
detach("package:grid")
detach("package:maps")

## Convert the results from LatticeKrig to esd:
W <- w$z


attr(W,'variable') <- varnm
attr(W,'unit') <- unitnm
attr(W,'longitude') <- w$x
attr(W,'latitude') <- w$y
class(W) <- class(etopo5)

## Make a projection that zooms in on the Barents region
#map(W,xlim=range(lon(W)),ylim=range(lat(W)),projection='sphere',colbar=list(n=21))

Wx <- max(abs(W),na.rm=TRUE)
rev <- switch(param,'t2m'=TRUE,'pr'=FALSE)
Wx <- max(abs(W),na.rm=TRUE)
if (min(W,na.rm=TRUE) < 0) {
  breaks <- round(seq(-Wx,Wx,length=31),2)
  pal <- switch(param,'t2m'='t2m','pr'='precip')
} else {
  breaks <- round(seq(0,Wx,length=31),2)
  pal <- switch(param,'t2m'='warm','pr'='precip')
}
mapcol <- colscal(n=length(breaks)-1,col=pal)

#attr(W,'variable') <- varid(X)[1]
#attr(W,'unit') <- unit(X)[1]

map(W,xlim=range(lon(W)),ylim=range(lat(W)),projection='sphere',
    colbar=list(col=mapcol,breaks=breaks,pal=pal,rev=rev))

dev.copy2pdf(file=paste('paper59_Fig_seasonalvar_',
               param,'_',FUN,'map.pdf',sep=''))
}
}

#analyseseason(param='t2m',FUN='mean',prop=FALSE,pattern=1,map.show=TRUE)
#analyseseason(param='t2m',FUN='sd',prop=FALSE,pattern=1)
analyseseason(param='t2m',FUN='ar1',prop=FALSE,pattern=1)
#analyseseason(param='pr',FUN='wetmean',prop=FALSE,pattern=1,map.show=TRUE)
#analyseseason(param='pr',FUN='wetmean',prop=FALSE,pattern=2)
#analyseseason(param='pr',FUN='wetfreq',prop=FALSE,pattern=1)
#analyseseason(param='pr',FUN='wetfreq',prop=FALSE,pattern=2)
#analyseseason(param='pr',FUN='wet.spell',prop=FALSE,pattern=1)
#analyseseason(param='pr',FUN='wet.spell',prop=FALSE,pattern=2)
#analyseseason(param='pr',FUN='dry.spell',prop=FALSE,pattern=1)
#analyseseason(param='pr',FUN='dry.spell',prop=FALSE,pattern=2)

#while (dev.cur()>1) dev.off()
#x11()
#source('paper59trends.R')
