## Analysis of trends in the liquid and solid the precipitation phases
## based on stations where both temperature and precipitation are recorded.
### Rasmus E. Benestad, Oslo, November, 2015

library(esd)

gridmap <- function(X,projection="sphere",
                    xlim=NULL,ylim=NULL,new=TRUE,
                    verbose=TRUE,...) {
  if (verbose) print("gridmap")
  require(LatticeKrig)
  data(etopo5)
  etopo5 <- subset(etopo5,is=list(lon=range(lon(X)),lat=range(lat(X))))
  etopo5[etopo5<=0] <- NA

  if (verbose) print(paste("Grid the data with LatticeKrig"))
  grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')
  ok <- !(duplicated(lon(X)) & duplicated(lat(X)))
  obj <- LatticeKrig( x=cbind(lon(X)[ok],lat(X)[ok]),
                      y=X[ok],Z=alt(X)[ok])
  w <- predictSurface(obj, grid.list = grid,Z=etopo5)
  w$z[is.na(etopo5)] <- NA

  if (verbose) print("Convert from LatticeKrig to esd")
  W <- w$z
  param <- varid(X)[1]
  units <- attr(X,"unit")[1]
  attr(W,'variable') <- paste(param,"trend")
  attr(W,'unit') <- " percent/decade"
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
  dev.new()
  colbar <- list(col=rev(colscal(n=20,col='t2m')),breaks=round(seq(-10,10,length.out=21),1))
  map(W,projection=projection,xlim=xlim,ylim=ylim,new=FALSE,colbar=colbar,
      verbose=verbose,...)
  invisible(W)
}



sl <- -2
rl <- 2

FUN <- 'wetfreq'
load('pr.aaca.rda')
load('t2m.aaca.rda')

## Remove duplicated records:
dt <- duplicated(stid(t2m.aaca))
t2m.aaca <- subset(t2m.aaca,is=!dt)
dp <- duplicated(stid(pr.aaca))
pr.aaca <- subset(pr.aaca,is=!dp)

## Find stations with both precipitation and temperature
ipr <- is.element(stid(pr.aaca),stid(t2m.aaca))
it2 <- is.element(stid(t2m.aaca),stid(pr.aaca))

pr <- subset(pr.aaca,is=ipr)
t2m <- subset(t2m.aaca,is=it2)
pr <- matchdate(pr,it=t2m)
t2m <- matchdate(t2m,it=pr)

PR <- coredata(pr); T2M <- coredata(t2m)

snow <- PR; rain <- PR; mixp <- PR
snow[T2M < sl] <- NA
rain[T2M > rl] <- NA
mixp[(T2M >= sl) & (T2M <= rl)] <- NA

Snow <- pr; coredata(Snow) <- snow
mu.s <- aggregate(Snow,year,'wetmean')
coredata(mu.s) <- round(100*t(t(coredata(mu.s))/apply(coredata(mu.s),2,'mean',na.rm=TRUE)),1)
trend.s <- apply(mu.s,2,'trend.coef')
nval <- apply(mu.s,2,'nv')
trend.s[nval < 30] <- NA
trend.s <- attrcp(mu.s,trend.s)
gridmap(trend.s,xlim=range(lon(mu.s)),ylim=range(lat(mu.s)))
figlab('Snow',ypos = 0.99)
dev.copy2pdf(file='paper59_mu.s_trend.pdf')

Rain <- pr; coredata(Rain) <- rain
mu.l <- aggregate(Rain,year,'wetmean')
coredata(mu.l) <- round(100*t(t(coredata(mu.l))/apply(coredata(mu.l),2,'mean',na.rm=TRUE)),1)
trend.l <- apply(mu.l,2,'trend.coef')
nval <- apply(mu.l,2,'nv')
trend.l[nval < 30] <- NA
trend.l <- attrcp(mu.l,trend.l)
gridmap(trend.l,xlim=range(lon(mu.l)),ylim=range(lat(mu.l)))
figlab('Rain',ypos = 0.99)
dev.copy2pdf(file='paper59_mu.l_trend.pdf')


Mixp <- pr; coredata(Mixp) <- mixp
mu.m <- aggregate(Mixp,year,'wetmean')
coredata(mu.l) <- round(100*t(t(coredata(mu.l))/apply(coredata(mu.l),2,'mean',na.rm=TRUE)),1)
trend.m <- apply(mu.m,2,'trend.coef')
nval <- apply(mu.m,2,'nv')
trend.m[nval < 30] <- NA
trend.m <- attrcp(mu.m,trend.m)
gridmap(trend.m,xlim=range(lon(mu.m)),ylim=range(lat(mu.m)))
figlab('Mixed phase',ypos = 0.99)
dev.copy2pdf(file='paper59_mu.m_trend.pdf')

## histograms of trend coefficients:
dev.new()
tx <- ceiling(quantile(abs(trend.s),0.97,na.rm=TRUE))
breaks <- seq(-tx,tx,length=31)
trend.s[abs(trend.s) > tx] <- NA

h <- hist(trend.s,breaks=breaks,col=colscal(30),freq=FALSE,
     main='Trends in snow amount',xlab='%/decade')
grid()
polygon(c(-tx,0,0,-tx,-tx),c(0,0,2*rep(max(h$density),2),0),col=rgb(0.5,0.5,0.5,0.2))
lines(h$mids,dnorm(h$mids,mean=mean(trend.s,na.rm=TRUE),sd=sd(trend.s,na.rm=TRUE)),
      lwd=5,col=rgb(0.5,0.3,0.3,0.25))
p.gt.0 <- round(pnorm(0,mean=mean(trend.s,na.rm=TRUE),sd=sd(trend.s,na.rm=TRUE)),3)
figlab(paste('Pr(X > 0)=',1-p.gt.0))
dev.copy2pdf(file='paper59_Fig_trend_snow_pdf-pca')

## histograms of trend coefficients:
dev.new()
tx <- ceiling(quantile(abs(trend.l),0.97,na.rm=TRUE))
breaks <- seq(-tx,tx,length=31)
trend.l[abs(trend.l) > tx] <- NA

h <- hist(trend.l,breaks=breaks,col=colscal(30),freq=FALSE,
     main='Trends in rain amount',xlab='%/decade')
grid()
polygon(c(-tx,0,0,-tx,-tx),c(0,0,2*rep(max(h$density),2),0),col=rgb(0.5,0.5,0.5,0.2))
lines(h$mids,dnorm(h$mids,mean=mean(trend.l,na.rm=TRUE),sd=sd(trend.l,na.rm=TRUE)),
      lwd=5,col=rgb(0.5,0.3,0.3,0.25))
p.gt.0 <- round(pnorm(0,mean=mean(trend.l,na.rm=TRUE),sd=sd(trend.l,na.rm=TRUE)),3)
figlab(paste('Pr(X > 0)=',1-p.gt.0))
dev.copy2pdf(file='paper59_Fig_trend_rain_pdf-pca')

## histograms of trend coefficients:
dev.new()
tx <- ceiling(quantile(abs(trend.m),0.97,na.rm=TRUE))
breaks <- seq(-tx,tx,length=31)
trend.m[abs(trend.m) > tx] <- NA

h <- hist(trend.m,breaks=breaks,col=colscal(30),freq=FALSE,
     main='Trends in mixed rain-snow amount',xlab='%/decade')
grid()
polygon(c(-tx,0,0,-tx,-tx),c(0,0,2*rep(max(h$density),2),0),col=rgb(0.5,0.5,0.5,0.2))
lines(h$mids,dnorm(h$mids,mean=mean(trend.m,na.rm=TRUE),sd=sd(trend.m,na.rm=TRUE)),
      lwd=5,col=rgb(0.5,0.3,0.3,0.25))
p.gt.0 <- round(pnorm(0,mean=mean(trend.m,na.rm=TRUE),sd=sd(trend.m,na.rm=TRUE)),3)
figlab(paste('Pr(X > 0)=',1-p.gt.0))
dev.copy2pdf(file='paper59_Fig_trend_mixed_pdf-pca')

T2m <- aggregate(subset(t2m,it=month.abb[5:9]),year,'mean')
x <- apply(T2m,2,'trend.coef')
nval <- apply(T2m,2,'nv')
x[nval < 30] <- NA
mu <- aggregate(subset(pr,it=month.abb[5:9]),year,'wetmean')
coredata(mu) <- round(100*t(t(coredata(mu))/apply(coredata(mu),2,'mean',na.rm=TRUE)),1)
y <- apply(mu,2,'trend.coef')
nval <- apply(mu,2,'nv')
y[nval < 30] <- NA

dev.new()
par(bty='n')

mx <- mean(x,na.rm=TRUE); sx <- sd(x,na.rm=TRUE)
my <- mean(y,na.rm=TRUE); sy <- sd(y,na.rm=TRUE)
s <- sin(seq(0,2*pi,length.out=360)); c <- cos(seq(0,2*pi,length.out=360))
plot(x,y,xlim=c(-1,1),ylim=c(-10,10),pch=19,cex=0.5,
     main='Temperature and wet-day mean precipitation trends',
     sub='May-September',
     xlab=expression(bar(x)),ylab=expression(paste(mu,' (%)')))
rect(-1,-10,0,10,col=rgb(0.5,0.5,1,0.2))
rect(-1,-10,1,0,col=rgb(0.5,0.5,0,0.2))
rect(0,0,1,10,col=rgb(1,0.5,0.5,0.1))

for (p in seq(0.9,0.1,by=-0.1)) {
  rx <- qnorm(p,sd=sx); ry <- qnorm(p,sd=sy)
  polygon(mx+rx*s,my+ry*c,col=rgb(0.5,0.5,+.5,0.2),border='grey')
}
points(x,y,pch=19,cex=0.5)
lines(c(-1,mx),rep(my,2),lty=2)
text(-0.5,my,pos=3,paste(round(my/mx,2),'mm/day per degree'))
  
dev.copy2pdf(file='paper59_Fig_trend_t2mmu.pdf')


