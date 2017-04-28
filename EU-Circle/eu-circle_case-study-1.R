## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE---------------------------------------------------------
## library(knitr)
## purl('~/git/esd_Rmarkdown/EU-Circle/eu-circle_case-study-1.Rmd', output='~/git/esd_Rmarkdown/EU-Circle/eu-circle_case-study-1.R')

## ----message=FALSE-------------------------------------------------------
library(esd)

## ----warning=FALSE-------------------------------------------------------
if (!file.exists('eu-circle_case-study-1_precio.rda')) {
  ss <- select.station(param='precip',lon=c(5,10),lat=c(42,45),nmin=50,src=c('ecad','ghcnd'))
  x <- station(ss)
  x <- subset(x,is=list(alt=-400))
  save(file='eu-circle_case-study-1_precio.rda',x)
} else load('eu-circle_case-study-1_precio.rda')
## Only keep the ECAD data
x <- subset(x,is=is.element(src(x),'ECAD'))
map(x,FUN='q95',new=TRUE)
diagnose(x)

## ----warning=FALSE-------------------------------------------------------
plot(x,new=FALSE)
grid()

## ----warning=FALSE-------------------------------------------------------
mu <- annual(x,'wetmean',nmin=250)
wq95 <- function(x) {x <- x[is.finite(x)]; x <- x[x >= 1]; wq95 <- quantile(x,probs=0.95); wq95}
q95 <- annual(x,'wq95')
plot(mu,new=FALSE)
for (i in 1:dim(mu)[2]) lines(trend(subset(mu,is=i)),lty=2)
grid()


## ------------------------------------------------------------------------
x1 <- zoo(subset(annual(x,'count',threshold=1),is=1))
y1 <- zoo(subset(mu,is=1))
plot(x1,y1,xlab=expression(n[X>1*mm]),ylab=expression(mu),xlim=c(30,150),ylim=c(5,30),
     main='Test: wet-day mean precipitation dependency on sample size')
for (i in 2:dim(mu)[2]) {
  x1 <- zoo(subset(annual(x,'count',threshold=1),is=i))
  y1 <- zoo(subset(mu,is=i))
  points(x1,y1)
}
grid()

## ----warning=FALSE-------------------------------------------------------
## check the relationship between wet-day 95-percentile and wet-day mean
plot(-log(0.05)*zoo(mu[,2]),zoo(q95[,1]),pch=19,xlim=c(0,120),ylim=c(0,120),
     main='Wet-day mean v.s. 95th wet percentile',
     xlab=expression(paste(-ln(1-0.95)*mu,' (mm/day)')), ylab=expression(paste(q[95],' (mm/day)')))
points(-log(0.05)*zoo(mu[,1]),zoo(q95[,2]),pch=19,col=rgb(1,0,0,0.2))
points(-log(0.05)*zoo(mu[,3]),zoo(q95[,3]),pch=19,col=rgb(0,0,1,0.2))
points(-log(0.05)*zoo(mu[,4]),zoo(q95[,4]),pch=19,col=rgb(0,1,0,0.2))
points(-log(0.05)*zoo(mu[,5]),zoo(q95[,5]),pch=19,col=rgb(0,0.7,0.7,0.2))
points(-log(0.05)*zoo(mu[,6]),zoo(q95[,6]),pch=19,col=rgb(0.7,0,0.7,0.2))
points(-log(0.05)*zoo(mu[,7]),zoo(q95[,7]),pch=19,col=rgb(0.5,0.5,0.5,0.2))
grid()
legend(0,120,loc(mu),pch=19,col=c('black',rgb(1,0,0,0.2),rgb(0,0,1,0.2),rgb(0,1,0,0.2),rgb(0,0.7,0.7,0.2),rgb(0.7,0,0.7,0.2),rgb(0.5,0.5,0.5,0.2)),bty='n',cex=0.5)

## ----warning=FALSE-------------------------------------------------------
fw <- annual(x,'wetfreq',nmin=250)
plot(fw,new=FALSE)
for (i in 1:dim(fw)[2]) lines(trend(subset(fw,is=i)),lty=2)
grid()

## ----warning=FALSE-------------------------------------------------------
x0 <- 100 #mm
Pr <- zoo(1-pbinom(0,size=365,prob=coredata(fw)*exp(-x0/coredata(mu))),order.by=year(fw))
class(Pr) <- class(mu)
Pr <- attrcp(mu,Pr)
attr(Pr,'variable') <- 'Pr'
attr(Pr,'unit') <- 'probability'
attr(Pr,'longname') <- paste('probability of at least one day with more than',x0,'mm') 
plot(Pr,ylab=expression(Pr(n>0)),xlab='',new=FALSE,errorbar=FALSE,
          main=paste('Probability of at least one day with more than',x0,'mm'))
for (i in 1:dim(Pr)[2]) lines(trend(subset(Pr,is=i)),lty=2)
grid()

## ----warning=FALSE-------------------------------------------------------
z <- subset(x,it=c(1900,2016),is=1)
ncd <- spell(z,threshold=1)
plot(ncd)
grid()

## ----warning=FALSE-------------------------------------------------------
hist(ncd)

## ----warning=FALSE-------------------------------------------------------
amncd <- annual(ncd)
plot(amncd,col=c(rgb(0,0,1),rgb(1,0.3,0)),new=FALSE,map.show=FALSE)
grid()
lines(trend(subset(amncd,is=1)),lty=2)
mcdd <- trend(subset(amncd,is=2))
lines(mcdd,lty=1)
print(attr(mcdd,'coefficients'))

## ----warning=FALSE-------------------------------------------------------
plot(1-pgeom(1:100,1/mean(mcdd)),type='l',xlab='consecutive dry days',ylab=expression(Pr(X > x)))
plot(index(mcdd),1-pgeom(30, 1/mcdd),type='l',xlab='',ylab=expression(Pr(X > 30*days)),
     main='Estimated probability for 30-day dry spell',ylim=c(0,0.05))
for (i in 1:dim(Pr)[2]) lines(trend(subset(Pr,is=i)),lty=2)
grid()

## ----warning=FALSE-------------------------------------------------------
X <- retrieve('data.ECAD/rr_0.25deg_reg.nc',lon=c(4,8),lat=c(43,46))
map(X,new=FALSE)

## ------------------------------------------------------------------------
mean.dry <- function(x,na.rm=TRUE) mean(subset(spell(x,threshold=1),is=2),na.rm=na.rm)
map(X,FUN='mean.dry',new=FALSE,colbar=list(pal='warm'))

## ------------------------------------------------------------------------
trend.dry <- function(x,na.rm=TRUE) trend.coef(subset(spell(x,threshold=1),is=2))
map(X,FUN='trend.dry',new=FALSE,colbar=list(pal='warm'))

## ---- fig.height=8-------------------------------------------------------
dryseason <- function(x,na.rm=TRUE) mean(subset(spell(x,threshold=1),is=2),na.rm=na.rm)
eof.dry <- EOF(annual(X,FUN='dryseason'))
plot(eof.dry)

