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



analysetrends <- function(param='pr',pal=NULL,col=NULL,FUN='wetmean',prop=FALSE,pattern=1) {

load(paste(param,'.aaca.rda',sep=''))
X <- eval(parse(text=paste(param,'.aaca',sep='')))
X <- subset(X,it=c(1960,2015))

if (FUN=='sd') X <- anomaly(X)
if ( (FUN!='wet.spell') & (FUN!='dry.spell'))
  Xc <- aggregate(X,year,FUN) else {
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
    Xc <- aggregate(xx,year,"mean")
}

if (FUN=='wetmean') {Xcc <- coredata(Xc); Xcc[Xcc > 20] <- NA; Xcc -> coredata(Xc); rm('Xcc')}
if (FUN=='wetfreq') {Xcc <- coredata(Xc); Xcc[Xcc <= 0] <- NA; Xcc -> coredata(Xc); rm('Xcc')}
if (FUN=='wet.spell') {Xcc <- coredata(Xc); Xcc[Xcc > 10] <- NA; Xcc -> coredata(Xc); rm('Xcc')}
if (FUN=='dry.spell') {Xcc <- coredata(Xc); Xcc[Xcc > 20] <- NA; Xcc -> coredata(Xc); rm('Xcc')}

if (prop) {
  coredata(Xc) <- 100*t(t(coredata(Xc))/apply(coredata(subset(Xc,it=c(1961,1990))),2,'mean',na.rm=TRUE))
  FUN <- paste(FUN,'prop',sep='')
  attr(Xc,'unit') <- '%'
}

trends <- apply(Xc,2,'trend.coef')
nval <- apply(Xc,2,'nv')
trends[nval < 30] <- NA

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
                't2m ar1'='',
                'pr wetmean'='mm/day',
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
attr(Xc,'longname') <- longnm
plot(Xc,col=col)
grid()
dev.copy2pdf(file=paste('paper59_Fig_trend_',
               param,'_',FUN,'.pdf',sep=''))


pca <- PCA(pcafill(Xc),n=12)
attr(pca,'variable') <- varnm
attr(pca,'unit') <- unitnm
attr(pca,'longname') <- FUN

plot(pca,pattern=pattern)
dev.copy2pdf(file=paste('paper59_Fig_pcatrend_',pattern,'_',
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

## If temperature and all values are of same sign, use a one-signed color scheme
map(W,xlim=range(lon(W)),ylim=range(lat(W)),projection='sphere',
    colbar=list(col=mapcol,breaks=breaks,pal=pal,rev=rev))

dev.copy2pdf(file=paste('paper59_Fig_trend_',param,'_',FUN,'map-pca',pattern,'.pdf',sep=''))

tx <- 1.2*quantile(abs(trends),0.97,na.rm=TRUE)
breaks <- seq(-tx,tx,length=31)
trends[abs(trends) > tx] <- NA

dev.new()
h <- hist(trends,breaks=breaks,col=colscal(30),freq=FALSE,
     main=paste(varid(Xc)[1],'trends'),xlab=paste(unit(Xc)[1],'/decade'))
grid()
polygon(c(-tx,0,0,-tx,-tx),c(0,0,2*rep(max(h$density),2),0),col=rgb(0.5,0.5,0.5,0.2))
lines(h$mids,dnorm(h$mids,mean=mean(trends,na.rm=TRUE),sd=sd(trends,na.rm=TRUE)),
      lwd=5,col=rgb(0.5,0.3,0.3,0.25))
p.gt.0 <- round(pnorm(0,mean=mean(trends,na.rm=TRUE),sd=sd(trends,na.rm=TRUE)),3)
figlab(paste('Pr(X > 0)=',1-p.gt.0))

dev.copy2pdf(file=paste('paper59_Fig_trend_',param,'_',FUN,'pdf-pca',pattern,'.pdf',sep=''))
}

## Inter-annual variability

iav <- function(param='pr',FUN='wetmean',FUNX='sd',pal=NULL,rev=NULL,col=NULL,prop=FALSE) {

load(paste(param,'.aaca.rda',sep=''))
X <- eval(parse(text=paste(param,'.aaca',sep='')))
X <- subset(X,it=c(1960,2015))

if (FUN=='sd') X <- anomaly(X)
if ( (FUN!='wet.spell') & (FUN!='dry.spell'))
  Xc <- aggregate(X,year,FUN) else {
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
    Xc <- aggregate(xx,year,"mean")
}

if (FUN=='wetmean') {Xcc <- coredata(Xc); Xcc[Xcc > 20] <- NA; Xcc -> coredata(Xc); rm('Xcc')}
if (FUN=='wetfreq') {Xcc <- coredata(Xc); Xcc[Xcc <= 0] <- NA; Xcc -> coredata(Xc); rm('Xcc')}
if (FUN=='wet.spell') {Xcc <- coredata(Xc); Xcc[Xcc > 10] <- NA; Xcc -> coredata(Xc); rm('Xcc')}
if (FUN=='dry.spell') {Xcc <- coredata(Xc); Xcc[Xcc > 20] <- NA; Xcc -> coredata(Xc); rm('Xcc')}

if (prop) {
  Xc <- round(100*t(t(Xc)/apply(coredata(subset(Xc,it=c(1961,1990))),2,'mean',na.rm=TRUE)),1)
  FUN <- paste(FUN,'prop',sep='')
  attr(Xc,'unit') <- "'%'"
}

nval <- apply(Xc,2,'nv')
z <- apply(Xc,2,FUNX,na.rm=TRUE)
z[nval < 30] <- NA

varnm <- switch(paste(param,FUN),
                't2m mean'=expression(T[2*m]),
                't2m sd'=expression(sigma[T]),
                't2m ar1'=expression(R(tau)),
                'pr wetmean'=expression(mu),
                'pr wetfreq'=expression(f[w]),
                'pr wet.spell'=expression(n[c*w*d]),
                'pr dry.spell'=expression(n[c*d*d]),
                'pr wetmeanprop'=expression(mu),
                'pr wetfreqprop'=expression(f[w]))
unitnm <- switch(paste(param,FUN),
                't2m mean'=expression(degree*C),
                't2m sd'=expression(degree*C),
                't2m ar1'='',
                'pr wetmean'='mm/day',
                'pr wetfreq'='fraction',
                'pr wet.spell'='days',
                'pr dry.spell'='days',
                'pr wetmeanprop'="'%'",
                'pr wetfreqprop'="'%'")
longnm <- switch(paste(param,FUN),
                't2m mean'='Mean temperature',
                't2m sd'='Temperature variability',
                't2m ar1'='Temperature persistence',
                'pr wetmean'='Precipitation intensity',
                'pr wetfreq'='Wet-day frequency',
                'pr wet.spell'='Mean wet-spell length',
                'pr dry.spell'='Mean dry-spell length',
                'pr wetmeanprop'='Proportional precipitation intensity',
                'pr wetfreqprop'='Proportional precipitation frequency')

## Get elevation data
data(etopo5)
etopo5 <- subset(etopo5,is=list(lon=range(lon(X)),lat=range(lat(X))))
etopo5[etopo5<=0] <- NA

## Grid the PCA pattern
require(LatticeKrig)

## Set the grid to be the same as that of etopo5:
grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')

## Flag dubplicated stations:
ok <- !(duplicated(lon(X)) & duplicated(lat(X))) & is.finite(z)

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
attr(W,'unit') <- attr(Xc,'unit')
attr(W,'longitude') <- w$x
attr(W,'latitude') <- w$y
class(W) <- class(etopo5)

## Make a projection that zooms in on the Barents region
#map(W,xlim=range(lon(W)),ylim=range(lat(W)),projection='sphere',colbar=list(n=21))

if (is.null(rev)) rev <- switch(param,'t2m'=TRUE,'pr'=FALSE)
Wx <- max(abs(W),na.rm=TRUE)
if (min(W,na.rm=TRUE) < 0) {
  breaks <- round(seq(-Wx,Wx,length=31),2)
  if (is.null(pal)) pal <- switch(param,'t2m'='t2m','pr'='precip')
} else {
  breaks <- round(seq(0,Wx,length=31),2)
  if (is.null(pal)) pal <- switch(param,'t2m'='warm','pr'='precip')
}
mapcol <- colscal(n=length(breaks)-1,col=pal)
if (rev) mapcol <- rev(mapcol)

#attr(W,'variable') <- varid(X)[1]
#attr(W,'unit') <- unit(X)[1]

## If temperature and all values are of same sign, use a one-signed color scheme
map(W,xlim=range(lon(W)),ylim=range(lat(W)),projection='sphere',
    colbar=list(col=mapcol,breaks=breaks,pal=pal))
if (FUNX=='sd') figlab('Magnitude of interannual variations',ypos=0.97) else
if (FUNX=='trend.coef') figlab(paste('Historic trends',start(X),'-',end(X)),ypos=0.97)
#lab <- eval(parse(text=paste('expression(',varnm,'*(',unitnm,'))')))
#figlab(lab)
dev.copy2pdf(file=paste('paper59_Fig_iav_',param,'_',FUN,'_',FUNX,'map.pdf',sep=''))
}


if (FALSE) {
analysetrends(param='t2m',FUN='mean',pattern=1)
analysetrends(param='t2m',FUN='mean',pattern=2)
analysetrends(param='t2m',FUN='sd',pattern=1)
analysetrends(param='t2m',FUN='sd',pattern=2)
analysetrends(param='t2m',FUN='ar1',pattern=1)
analysetrends(param='t2m',FUN='ar1',pattern=2)
analysetrends(param='pr',FUN='wetmean',prop=TRUE,pattern=1)
analysetrends(param='pr',FUN='wetmean',prop=TRUE,pattern=2)
analysetrends(param='pr',FUN='wetfreq',prop=TRUE,pattern=1)
analysetrends(param='pr',FUN='wetfreq',prop=TRUE,pattern=2)
analysetrends(param='pr',FUN='wet.spell',pattern=1)
analysetrends(param='pr',FUN='wet.spell',pattern=2)
analysetrends(param='pr',FUN='dry.spell',pattern=1)
analysetrends(param='pr',FUN='dry.spell',pattern=2)
}


#iav(param='t2m',FUN='mean')
#iav(param='t2m',FUN='sd')
#iav(param='pr',FUN='sum',pal='t2m',rev=TRUE)
#iav(param='pr',FUN='wetmean',pal='t2m',rev=TRUE)
#iav(param='pr',FUN='wetfreq',pal='t2m',rev=TRUE)
#iav(param='pr',FUN='wet.spell',pal='t2m',rev=TRUE)
#iav(param='pr',FUN='dry.spell',pal='t2m',rev=TRUE)
#iav(param='pr',FUN='sum',FUNX='trend.coef',pal='t2m',rev=TRUE)
iav(param='pr',FUN='wetmean',FUNX='trend.coef',pal='t2m',rev=TRUE,prop=TRUE)
iav(param='pr',FUN='wetfreq',FUNX='trend.coef',pal='t2m',rev=TRUE,prop=TRUE)
#iav(param='pr',FUN='wet.spell',FUNX='trend.coef',pal='t2m',rev=TRUE)
#iav(param='pr',FUN='dry.spell',FUNX='trend.coef',pal='t2m',rev=TRUE)

load('t2m.aaca.rda')
X <- coredata(annual(t2m.aaca,'mean'))
Y <- coredata(annual(anomaly(t2m.aaca),'sd'))
ok <- (apply(X,2,nv)>30) 
x <- apply(X[,ok],2,'trend.coef')
y <- apply(Y[,ok],2,'trend.coef')

dev.new()
par(bty='n')
mx <- mean(x,na.rm=TRUE); sx <- sd(x,na.rm=TRUE)
my <- mean(y,na.rm=TRUE); sy <- sd(y,na.rm=TRUE)
s <- sin(seq(0,2*pi,length.out=360)); c <- cos(seq(0,2*pi,length.out=360))
plot(x,y,pch=19,xlim=c(-2,2),ylim=c(-0.5,0.5),cex=0.5,
     main='Temperature trends: mean and variability',
     xlab=expression(bar(x)),ylab=expression(sigma[T]))
rect(-2,-0.5,0,0.5,col=rgb(0.5,0.5,1,0.2))
rect(-2,-0.5,2,0,col=rgb(0.5,0.5,0,0.2))
rect(0,0,2,0.5,col=rgb(1,0.5,0.5,0.1))

for (p in seq(0.9,0.1,by=-0.1)) {
  rx <- qnorm(p,sd=sx); ry <- qnorm(p,sd=sy)
  polygon(mx+rx*s,my+ry*c,col=rgb(0.5,0.5,+.5,0.2),border='grey')
}
points(x,y,pch=19,cex=0.5)
lines(c(-2,mx),rep(my,2),lty=2)
text(-1.25,my,pos=3,paste(round(my/mx,2),'degree/degree'))

dev.copy2pdf(file='paper59_Fig_trend_t2m_meansigma.pdf')

