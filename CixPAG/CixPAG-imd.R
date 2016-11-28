library(esd)
 
## Reading data predictand: maximum temperature

tx <- retrieve('data/IMD/maxT_IMD.1951-2013.nc')
attr(tx,'unit') <- 'degC'
map(tx)
 
## Some initial exploration of the data - the mean summer maximum daily maximum temperature

Tx.maxjja <- subset(as.4seasons(tx,FUN='max'),it='jja')
index(Tx.maxjja) <- year(Tx.maxjja)
eof.Y <- EOF(Tx.maxjja,n=5)
plot(eof.Y,new=FALSE)
 
## Trend in the men summer maximum daily maximum temperature

map(Tx.maxjja,FUN='trend',new=FALSE)
 
## The mean number of days with daily maximum temperature exceeding 40C. 

Tx.nhd40 <- annual(tx,FUN='count',threshold=40)
map(Tx.nhd40,new=FALSE)
 
## Downscaling

predictor <- retrieve('air.mon.mean.nc',lon=c(60,100),lat=c(0,25))
predictor <- annual(subset(predictor,it='jja'),nmin=3)
eof.X <- EOF(predictor,n=10)
plot(eof.X,new=FALSE)
 
## Check the link between large and small scales

class(eof.Y) <- c("eof", "field", "season", "zoo") # Fix some incorrect settings
class(eof.X) <- c("eof", "field", "season", "zoo") # Fix some incorrect settings
index(eof.Y) <- as.Date(paste(year(eof.Y),'-07-01',sep=''))
index(eof.X) <- as.Date(paste(year(eof.X),'-07-01',sep=''))
ds.test <- DS(eof.Y,eof.X,ip=1:6)
plot(ds.test,new=FALSE)
 

### Apply the downscaling to the CMIP5 ensemble RCP4.5 

## Carry out the downscaling once and save the results. 

if (!file.exists('dse.Tx.jja.cixpag.rda')) {
  dse.Tx <- DSensemble.eof(eof.Y,eof.X,ip=1:6,it='jja',
                           predictor='air.mon.mean.nc',verbose=FALSE)
  save(dse.Tx,file='dse.Tx.jja.cixpag.rda')
} else load('dse.Tx.jja.cixpag.rda')

