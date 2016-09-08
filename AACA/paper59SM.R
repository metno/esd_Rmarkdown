## Rasmus E. Benestad
## Test of pcafill

library(esd)
param <- 't2m'
FUN <- 'mean'
lon=c(0,100); lat=c(65,90)
reanalysis <- 'air.mon.mean.nc'
period=c(1950,2015)

if (param=='precip') FUNX <- 'C.C.eq' else FUNX <- 'mean'

## Get the predictand -> Y
if (!exists('Y')) {
  load(paste(param,'.aaca.rda',sep=''))
  Y <- eval(parse(text=paste(param,'.aaca',sep='')))
}
Y <- subset(Y,is=list(lat=lat,lon=lon))

## Use a time and space window:  
for (it in c('djf','mam','jja','son')) {
  y <- subset(Y,it=period)

## Estimate seasonal means & weed out stations with little data
  Y4 <- subset(as.4seasons(y,FUN=FUN),it=it)
  ok <- apply(coredata(Y4),1,nv)
  Y4 <- subset(Y4,it=ok>0)
  nok <- apply(coredata(Y4),2,nv)
  Y4 <- subset(Y4,is=nok>15)

  map(y,FUN=FUN,cex=-1.5)
  dev.copy2pdf(file=paste('paper59-stations.mean',param,'.',it,'pdf',sep=''))

## Fill missing data using PCA-based regression
  dev.new()
  Z <- pcafill.test(Y4,max.miss=50)
  dev.copy2pdf(file=paste('paper59-pcafill.test.',param,'.',it,'.pdf',sep=''))
}
