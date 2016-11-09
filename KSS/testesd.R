param <- 't2m'
it <- 'djf'
period <- c(1950,2015)
n <- 3
ip <- 1:2
landmask <- FALSE

if (param=='t2m') {
  #reanalysis <- 'air.mon.mean.nc'
  reanalysis <- 'ERA40_t2m_mon.nc'
  FUN='mean'
  FUNX='mean'
  pattern='tas_Amon_ens'
  predictand <- '/lustre/storeB/users/rasmusb/data/t2m.nordic.rda'
  varid <- 'tg'
  #lon=c(-10,40); lat=c(50,75)
  lon=c(-5,40); lat=c(55,72)
}
if (param=='mu') {
  reanalysis <- 'air.mon.mean.nc'
  FUN='wetmean'
  FUNX='C.C.eq'
  pattern='tas_Amon_ens'
  predictand <- '/lustre/storeB/users/rasmusb/data/rr.nordic.rda'
  varid <- 'precip'
  lon=c(-90,10); lat=c(25,75)
  landmask <- TRUE
}
if (param=='fw') {
  reanalysis <- 'slp.mon.mean.nc'
  FUN='wetfreq'
  FUNX='mean'
  pattern='psl_Amon_ens'
  predictand <- '/lustre/storeB/users/rasmusb/data/rr.nordic.rda'
  varid <- 'precip'
  lon=c(-30,30); lat=c(50,75)
}

outdir <- '/lustre/storeB/users/rasmusb/dse4KSS/'
if (!file.exists(outdir)) dir.create(outdir)

verbose=FALSE

## Retrieving the data
## Get the predictand data: daily temperature from ECAD:

print('Retrieve data')

if (!file.exists(predictand)) {
  ss <- select.station(param=varid,lon=c(5,30),lat=c(55,72),src='ecad',nmin=50)
  Y <- station(ss)
  save(Y,file=predictand)
} else load(predictand)

## Only use the most recent 66 years 
Y <- subset(Y,it=c(1960,2015))
## Remove records with many missing datapoints
nv <- apply(coredata(Y),2,'nv')
Y <- subset(Y,is=nv > 14000)

print('predictor')

## Get the large-scale predictor:
T2M <- retrieve(reanalysis,lon=lon,lat=lat)
predictor <- subset(as.4seasons(T2M,FUNX=FUNX),it=it)

if (landmask) predictor <- mask(predictor,land=TRUE)
eof <- EOF(predictor)

Y <- subset(Y,it=period)
  
## Estimate seasonal means & weed out stations with little data
Y4 <- subset(as.4seasons(Y,FUN=FUN),it=it)
ok <- apply(coredata(Y4),1,nv)
Y4 <- subset(Y4,it=ok>0)
nok <- apply(coredata(Y4),2,nv)
Y4 <- subset(Y4,is=nok>15)
  
print(paste(round(100*sum(!is.finite(Y4))/length(!is.finite(Y4))),'% missing',
              sep=''))
  
nmiss <- round(100*sum(!is.finite(Y4))/length(!is.finite(Y4)))
print(paste(nmiss,'% missing',sep=''))
  
## Fill missing data using PCA-based regression
Z <- pcafill(Y4)
  
pca <- PCA(Z,n=n)

ds.test <- DS(pca,eof,ip=ip)
diag <- diagnose(ds.test)
xcorr <- round(100*cor(diag$xval$X.PCA.1,diag$xval$Z.PCA.1))
plot(ds.test)

save(ds.test,file=paste(outdir,'ds.test.',param,'.',it,
                 '.',lon[1],'-',lon[2],'.',lat[1],'-',lat[2],
                 '.xcorr',xcorr,'.rda',sep=''))
