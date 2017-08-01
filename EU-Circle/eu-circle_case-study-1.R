## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE---------------------------------------------------------
## ## Extract just the R-code
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

cols <- c(rgb(1,0,0,0.2),rgb(0,0,1,0.2),rgb(0,1,0,0.2),rgb(0,0.7,0.7,0.2),rgb(0.7,0,0.7,0.2),rgb(0.5,0.5,0.5,0.2))
## check the relationship between wet-day 95-percentile and wet-day mean
plot(-log(0.05)*zoo(mu[,1]),zoo(q95[,1]),pch=19,xlim=c(0,120),ylim=c(0,120),
     main='Wet-day mean v.s. 95th wet percentile',
     xlab=expression(paste(-ln(1-0.95)*mu,' (mm/day)')), ylab=expression(paste(q[95],' (mm/day)')))
for (i in 1:dim(mu[2])) points(-log(0.05)*zoo(mu[,i]),zoo(q95[,i]),pch=19,col=cols[i])
grid()
legend(0,120,loc(mu),pch=19,col=c('black',cols[1:dim(mu)[2]]),bty='n',cex=0.5)

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
X <- retrieve('~/data/data.ECAD/rr_0.25deg_reg.nc',lon=c(4,8),lat=c(43,46))
map(X,new=FALSE)

## ------------------------------------------------------------------------
mean.dry <- function(x,na.rm=TRUE) mean(subset(spell(x,threshold=1),is=2),na.rm=na.rm)
map(X,FUN='mean.dry',new=FALSE,colbar=list(pal='warm'))

## ------------------------------------------------------------------------
trend.dry <- function(x,na.rm=TRUE) trend.coef(subset(spell(x,threshold=1),is=2))
map(X,FUN='trend.dry',new=FALSE,colbar=list(pal='warm'))

## ---- fig.height=8-------------------------------------------------------
dryseason <- function(x,na.rm=TRUE) if (sum(is.finite(x))>1) mean(subset(spell(x,threshold=1),is=2),na.rm=na.rm) else NA
if (!file.exists('case1-dryseason.rda')) {
  z <- aggregate(subset(X,it='mjjas'),year,FUN='dryseason')
  save(z,file='case1-dryseason.rda')
} else load('case1-dryseason.rda')
attr(z,'variable') <- 'duration'
attr(z,'longname') <- 'mean dry spell duration'
eof.dry <- EOF(z)
plot(eof.dry,new=FALSE)

## ------------------------------------------------------------------------
eof.dry <- subset(eof.dry,ip=1:4)
index(eof.dry) <- year(eof.dry)

## ------------------------------------------------------------------------
T2M <- aggregate(subset(retrieve('~/Downloads/air.mon.mean.nc',lon=c(-12,35),lat=c(27,50)),it='mjjas'),year,'mean')
SLP <- aggregate(subset(retrieve('~/Downloads/slp.mon.mean.nc',lon=c(-12,35),lat=c(27,50)),it='mjjas'),year,'mean')
eof.t2m <- EOF(T2M)
eof.slp <- EOF(SLP)

## ---- fig.height=8, warning=FALSE----------------------------------------
ds.t2m <- DS(eof.dry,eof.t2m)
plot(ds.t2m,new=FALSE)

## ---- fig.height=8, warning=FALSE----------------------------------------
ds.slp <- DS(eof.dry,eof.slp)
plot(ds.slp,new=FALSE)

## ---- fig.height=8, warning=FALSE----------------------------------------
cca.slp.t2m <- CCA(eof.slp,eof.t2m)
plot(cca.slp.t2m,new=FALSE)

## ---- fig.width=8, fig.height=8------------------------------------------
## Estimate the residuals
eof.res <- subset(eof.dry,it=ds.slp)
coredata(eof.res) <- coredata(eof.res) - coredata(ds.slp)

## ---- fig.height=8, warning=FALSE----------------------------------------
ds.res <- DS(eof.res,eof.t2m)
plot(ds.res,new=FALSE)

## ------------------------------------------------------------------------
T2M <- retrieve('~/Downloads/air.mon.mean.nc',lon=c(-12,35),lat=c(27,50))
SLP <- retrieve('~/Downloads/slp.mon.mean.nc',lon=c(-12,35),lat=c(27,50))
## Also - fix the index of eof - set to year.
index(eof.dry) <- year(eof.dry)

## ---- warning=FALSE------------------------------------------------------
if (!file.exists('dse.mld.slp.rda')) {
  dse.mld.slp <- DSensemble.eof(eof.dry,predictor=SLP,it='mjjas',ip=1:10,path="~/data/CMIP5.monthly",pattern = "psl_Amon_ens_")
  save(dse.mld.slp,file = 'dse.mld.slp.rda')
} else load('dse.mld.slp.rda')
plot(dse.mld.slp,new=FALSE)
grid()

## ------------------------------------------------------------------------
 zmap1 <- map(dse.mld.slp,FUN="mean",FUNX='mean',it=c(2000,2009),plot=FALSE)
 zmap2 <- map(dse.mld.slp,FUN="mean",FUNX='mean',it=c(2090,2099),plot=FALSE)
 zmap <- zoo(coredata(zmap2) - coredata(zmap1),order.by=index(zmap2))
 zmap <- attrcp(zmap1,zmap)
 class(zmap) <- class(zmap1)
 attr(zmap,'variable') <- 'mean dry duration'
 attr(zmap,'unit') <- 'days'
 rm('zmap1','zmap2')
 map(zmap,colbar=list(pal='warm'),main='Change in mean duration of dry spells from 2000-2009',new=FALSE)

## ------------------------------------------------------------------------
## Read wild fire data
colnames <- c('year','number','type','department','code-INSEE','municipality','something','DFCI.code','Alert','origin','area')
fires <- read.table('~/Downloads/liste_incendies_ du_09_05_2017.csv',sep=';',col.names=colnames,header=FALSE)
h <- hist(fires$year,col='grey',breaks=1973:2016)

## ------------------------------------------------------------------------
breaks <- seq(0,3000,by=250)
lambda <- mean(h$counts)
hist(h$counts,breaks=breaks,freq = FALSE,ylim=c(0,0.02))
lines(seq(0,3000,by=10),dpois(seq(0,3000,by=10),lambda=lambda),lwd=3,col='red')

## ------------------------------------------------------------------------
summary(fires$area)
plot(fires$area,log='y')

## ------------------------------------------------------------------------
## Number of annual fires
naf <- zoo(h$counts,order.by=trunc(h$mids))
naf[naf==0] <- NA
index(naf) <- year(naf)
xy <- merge(naf,zoo(eof.dry),all=FALSE)
cal.naf <- data.frame(y = coredata(xy$naf), x1=coredata(xy$x.1), x2=coredata(xy$x.2),
                      x3= coredata(xy$x.3), x4=coredata(xy$x.4))
naf.model <- lm(y ~ x1 + x2 + x3 + x4, data=cal.naf)
z.naf <- zoo(predict(naf.model),order.by=trunc(h$mids))
plot(naf,z.naf)

## ------------------------------------------------------------------------
firearea <- data.frame(i=rep(1,length(fires$year)),area=as.numeric(as.character(fires$area)),year=fires$year)
## Only include fires where more than 1 km^2 has burned down
Af <- aggregate(firearea,by=list(firearea$year,firearea$area > 100000),FUN=sum)

## Number of events with Af > 1km^2 per year
plot(Af$Group.1[Af$Group.2],Af$i[Af$Group.2])

## Number of events with Af > 1km^2 per year
plot(Af$Group.1[Af$Group.2],Af$area[Af$Group.2])

## Relationship between the number of events and total area.
col <- rgb((2015 - Af$Group.1[Af$Group.2])/42,0,1- (2015 - Af$Group.1[Af$Group.2])/42,0.4)
plot(Af$i[Af$Group.2],Af$area[Af$Group.2],col=col,pch=19)
## annual total area of fires
ataf <- as.station(zoo(Af$area[Af$Group.2],order.by=Af$Group.1[Af$Group.2]),loc='Southern France',param='area',unit='m^2',
                  longname='burned area')
plot(ataf,new=FALSE)

## ------------------------------------------------------------------------
index(ataf) <- year(ataf)
index(eof.dry) <- year(eof.dry)
eof.dry <- subset(eof.dry,ip=1:4)
ds.ataf.dry <- DS(log(ataf),eof.dry)
ds.ataf.t2m <- DS(log(ataf),eof.t2m)
ds.ataf.slp <- DS(ataf,eof.slp)

## ---- fig.height=8-------------------------------------------------------
plot(ds.ataf.dry,new=FALSE)

## ---- fig.height=8-------------------------------------------------------
plot(ds.ataf.t2m,new=FALSE)

## ---- fig.height=8-------------------------------------------------------
plot(ds.ataf.slp,new=FALSE)

## ---- fig.height=8-------------------------------------------------------
## Convert the projected data into an EOF for the multi-model ensemble
eof.dry.z <- as.eof(dse.mld.slp)
## Use this EOF as input with the downscaling model calibrated with mean dry interval duration to make projections of ln(fire area)
ds.Af <- exp(predict(ds.ataf.dry,newdata=eof.dry.z))
plot(ds.Af,new=FALSE)

## ------------------------------------------------------------------------
fevents <- data.frame(i=rep(1,length(fires$year)),area=as.numeric(as.character(fires$area)),year=fires$year,date=as.Date(fires$Alert))
## Only include fires where more than 1 km^2 has burned down
se <- aggregate(fevents$i,by=list(fevents$date,fevents$area > 100000),FUN=sum)

## Examine the statistical distribution of simultaneous events
hist(se$x, col='grey', main='Number of simultaneous wild fires')


plot(year(se$Group.1),se$x,pch=19,col=rgb(1,0,0,0.2),main='Number of simultaneous wild fires by year',ylab='Simultaneous events',xlab='')
points(year(se$Group.1),se$x,pch=21,col=rgb(0.5,0.5,0.5))
grid()

