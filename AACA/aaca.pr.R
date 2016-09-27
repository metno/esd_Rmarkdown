## R-script to prepare data for ESD for AACA

## Science questions:
## What are the recent trends in temperature, precipitation, and snow depth?
## Is there a connection between wet-day mean and temperature?
## What are the main characteristics of downscaled PCA for T(2m)?
## Can we infer implications for number of cold days, warm days?
## Can we infer trends in wet-day  mean from temperature?
## What is the relationship between length of snow season and temperature?
## What is the relationship between max snow depth and temperature?
## Are there differences between Tmax, Tmin and daily mean T(2m) in terms of trends?

## Maps based on LatticeKrig & etopo5

library(esd)

if (!file.exists('pr.aaca.rda')) {
  ## Get the Arctic rain gauge data
  ss <- select.station(src=c('ecad','ghcnd'),
                       cntr=c('Norway','Sweden','Finland','Russia'),
                       lat=c(65,90),lon=c(0,90),nmin=40,param='rr')
  pr.aaca <- station(ss)
  ## Weed out stations with large missing data gaps & emphasis 1961--2015:
  pr.aaca <- subset(pr.aaca,it=c(1961,2015))
  n <- apply(coredata(pr.aaca),2,nv)
  pr.aaca <- subset(pr.aaca,is=n > 365.25*35)
  diagnose(pr.aaca)
  save(file='pr.aaca.rda',pr.aaca)
} else load('pr.aaca.rda')

mu <- annual(pr.aaca,FUN='wetmean')
fw <- annual(pr.aaca,FUN='wetfreq')

pca.mu <- PCA(allgood(mu))
pca.fw <- PCA(allgood(fw))

if (!file.exists('t2m.aaca.rda')) {
  ## Get the Arctic rain gauge data
  ss <- select.station(src=c('ecad','ghcnd'),
                       cntr=c('Norway','Sweden','Finland','Russia'),
                       lat=c(65,90),lon=c(0,90),nmin=40,param='t2m')
  t2m.aaca <- station(ss)
  ## Weed out stations with large missing data gaps & emphasis 1961--2015:
  t2m.aaca <- subset(t2m.aaca,it=c(1961,2015))
  n <- apply(coredata(t2m.aaca),2,nv)
  t2m.aaca <- subset(t2m.aaca,is=n > 365.25*35)
  diagnose(t2m.aaca)
  save(file='t2m.aaca.rda',t2m.aaca)
} else load('t2m.aaca.rda')

T2M <- annual(t2m.aaca,FUN='mean')
T2S <- annual(anomaly(t2m.aaca),FUN='sd')
pca.t2m <- PCA(allgood(T2M))
pca.t2s <- PCA(allgood(T2S))

cca.mu.t2m <- CCA(pca.mu,pca.t2m)
plot(cca.mu.t2m)

if (!file.exists('t2m.aaca.rda')) {
  ## Get the Arctic rain gauge data
  ss <- select.station(src='ecad',
                       cntr=c('Norway','Sweden','Finland','Russia'),
                       lat=c(65,90),lon=c(0,90),nmin=40,param='sd')
  sd.aaca <- station(ss)
  ## Weed out stations with large missing data gaps & emphasis 1961--2015:
  sd.aaca <- subset(sd.aaca,it=c(1961,2015))
  n <- apply(coredata(sd.aaca),2,nv)
  sd.aaca <- subset(sd.aaca,is=n > 365.25*35)
  diagnose(sd.aaca)
  save(file='sd.aaca.rda',sd.aaca)
} else load('sd.aaca.rda')

SDM <- annual(sd.aaca,FUN='max')
SDN <- annual(anomaly(sd.aaca),FUN='count',threshold=1)
pca.sdm <- PCA(allgood(SDM))
pca.sdn <- PCA(allgood(SDN))
