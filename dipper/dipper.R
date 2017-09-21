## R-script to analyse dipper population

library(esd)
library(ncdf4)

examine <- FALSE

print('read dipper data')
dipper <- read.table('~/Dropbox/data/dipper.csv',header=TRUE)

dipper <- zoo(x=dipper[[2]],order.by=dipper[[1]])
y <- as.station(dipper,loc='Lyngdalselva',
                lon=7,lat=58.5,alt=520,
                param='population',unit='count',
                reference='Marlène Gamelon')

predictor <- retrieve("air.mon.mean.nc",lon=c(-1,17),lat=c(57,62))
## Only use the winter temperature:
predictor <- aggregate(subset(predictor,it='djf'),year,FUN='mean')

if (examine) {
  print('dipper statistical distribution')
  ## Check the distribution
  n <- seq(0,150,by=10); mu <- mean(dipper)
  hist(coredata(dipper),breaks=n,freq=FALSE,col='grey')
  lines(x,dpois(x,mu),lwd=2,col='red')

  ## Single ESD for inspection
  X <- EOF(predictor)
  z <- DS(y,X)
  plot(z)
}

## Downscale the dipper population directly based on the large-scale
## annual mean temperature

print('Downscaled results')
if (!file.exists('dipper.Z.rcp45.rda')) {
  Z.rcp45 <- DSensemble.annual(y,biascorrect=TRUE,
                               predictor=predictor,
                               lon=c(-1,17),lat=c(57,62),
                               abscoords=TRUE)
  save(file='dipper.Z.rcp45.rda',Z.rcp45)
} else load('dipper.Z.rcp45.rda')

if (!file.exists('dipper.Z.rcp85.rda')) {
  Z.rcp85 <- DSensemble.annual(y,biascorrect=TRUE,
                               rcp="rcp85",predictor=predictor,
                               lon=c(-1,17),lat=c(57,62),
                               abscoords=TRUE)
  save(file='dipper.Z.rcp85.rda',Z.rcp85)
} else load('dipper.Z.rcp85.rda')

if (!file.exists('dipper.Z.rcp26.rda')) {
  Z.rcp26 <- DSensemble.annual(y,biascorrect=TRUE,
                               rcp="rcp26",predictor=predictor,
                               lon=c(-1,17),lat=c(57,62),
                               abscoords=TRUE)
  save(file='dipper.Z.rcp26.rda',Z.rcp26)
} else load('dipper.Z.rcp26.rda')

print('ensemble statistics')
year <- year(Z.rcp45)
ci90.rcp45 <- apply(coredata(Z.rcp45),1,quantile,
                    probs=c(0.05,0.95),na.rm=TRUE)
ci90.rcp26 <- apply(coredata(Z.rcp26),1,quantile,
                    probs=c(0.05,0.95),na.rm=TRUE)
ci90.rcp85 <- apply(coredata(Z.rcp85),1,quantile,
                    probs=c(0.05,0.95),na.rm=TRUE)

print('plotting')
dev.new()
par(bty='n')
plot(range(year),range(ci90.rcp45,ci90.rcp26,ci90.rcp85,y),
     type='n',xlab='',ylab='Population',main='Dipper',
     sub=loc(y))
grid()
polygon(c(year,rev(year)),c(ci90.rcp85[1,],rev(ci90.rcp85[2,])),
                            col=rgb(1,0.5,0,0.3),border=rgb(1,0.5,0))
polygon(c(year,rev(year)),c(ci90.rcp45[1,],rev(ci90.rcp45[2,])),
                            col=rgb(0.5,1,0,0.3),border=rgb(0.5,1,0))
polygon(c(year,rev(year)),c(ci90.rcp26[1,],rev(ci90.rcp26[2,])),
                            col=rgb(0,0.5,1,0.3),border=rgb(0,0.5,1))
lines(year(y),coredata(y),lwd=4,pch=19,type='b')

legend(1900,200,c('RCP2.6','RCP4.5','RCP8.5'),
       col=c(rgb(0,0.5,1),rgb(0.5,1,0),rgb(1,0.5,0)),
       lwd=7,lty=1,bty='n')
legend(1905,170,'observed',lwd=4,pch=19,lty=1,bty='n')

dev.copy2pdf(file='dipper.pdf')
