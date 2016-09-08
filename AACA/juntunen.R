## R script to downscale temperature etc for kittila for the AACA report.

library(esd)

#https://no.wikipedia.org/wiki/Kittil%C3%A4
#Koordinater: 67°38′48,678″N 24°54′21,622″Ø

if (!file.exists('kittila.rda')) {
  kittila <- read.table('~/Dropbox/Public/AACA/kittila_timeseries.csv',sep=';',
                        as.is=TRUE,col.names=c('date','precip','t2m','SWE'))
  t2m <- zoo(kittila$t2m,order.by=as.Date(kittila$date))
  precip <- zoo(kittila$precip,order.by=as.Date(kittila$date))
  SWE <- zoo(kittila$SWE,order.by=as.Date(kittila$date))

  t2m <- as.station(t2m,loc='Kittila',lon=24.90,lat=67.63,alt=196,
                    unit='deg C',param='temperature')
  precip <- as.station(precip,loc='Kittila',lon=24.90,lat=67.63,alt=196,
                       unit='mm',param='precipitation')
  SWE <- as.station(SWE,loc='Kittila',lon=24.90,lat=67.63,alt=196,
                    unit='mm',param='Snow depth water equivalent')
  save(file='kittila.rda',t2m,precip,SWE)
} else load('kittila.rda')

## Projections for kittila temperature
if (!file.exists('dse.kittila.rda')) {
  dse.kittila <- DSensemble.t2m(t2m,biascorrect=TRUE,verbose=FALSE)
  save(file='dse.kittila.rda',dse.kittila)
} else load('dse.kittila.rda')


## Historic trends in temperature, precip (mu and fw), and SWE.

mu <- annual(precip,FUN='wetmean')
fw <- annual(precip,FUN='wetfreq')

plot(annual(t2m))
lines(trend(annual(t2m)))

plot(mu)
lines(trend(mu))

plot(fw)
lines(trend(fw))

## Annual maximum snow depth
plot(annual(SWE,FUN='max'))
lines(trend(annual(SWE,FUN='max')))

## Number of days with snow
plot(annual(SWE,FUN='count',threshold=1))
lines(trend(annual(SWE,FUN='count',threshold=1)))

## Use seasonal variations in mu or variations in geographic vaiations in
## mean climatology to project future mu.

Mu <- aggregate(precip,month,FUN='wetmean')
plot(Mu)

plot(zoo(aggregate(t2m,month)),zoo(Mu))
