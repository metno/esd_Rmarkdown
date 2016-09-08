# Read the data from AARI/EALAT and store as 'esd' station data.

library(esd)

## The Russian data are stored in a annoying and cumbersome format
## without a consistent standard set-up
readEALAT <- function(fname,param='t2m',verbose=FALSE) {
  skip=0; eoh <- FALSE
  ## Check the length of the header
  if (verbose) {print(fname); print(param)}
  while (!eoh) {
    skip=skip + 1
    z <- readLines(fname,n=skip)
    z <- z[length(z)]
    eoh <- length(grep('/',z)>0)
  }
  ## Clean the data:
  h <- readLines(fname,n=skip)
  ZZ <- readLines(fname)
  zz <- ZZ
  for (i in 1:length(zz)) zz[i] <- gsub('*',' ',ZZ[i],fixed=TRUE,perl=TRUE)
  for (i in 1:length(zz)) zz[i] <- gsub('?',' ',ZZ[i],fixed=TRUE,perl=TRUE)
  for (i in 1:length(zz)) zz[i] <- gsub('+',' ',ZZ[i],fixed=TRUE,perl=TRUE)
  for (i in 1:length(zz))
    zz[i] <- gsub('9999.9','    NA',zz[i],fixed=TRUE,perl=TRUE)
  for (i in 1:length(zz))
    zz[i] <- gsub('999.9','   NA',zz[i],fixed=TRUE,perl=TRUE)
  for (i in 1:length(zz))
    zz[i] <- gsub(' 999','  NA',zz[i],fixed=TRUE,perl=TRUE)

  width <- switch(tolower(param),
                  't2m'=c(5,6,3,3,6,6,6),
                  'precip'=c(5,6,3,3,6),
                  'sd'=c(4,3,3,4),
                  'tx'=c(5,6,3,3,6,6,6),
                  'tn'=c(5,6,3,3,6,6,6))
  iyr <- switch(tolower(param),
                  't2m'=2,'precip'=2,'sd'=1,
                  'tx'=2,'tn'=2)
  nc <- nchar(zz[(skip+1):length(zz)])
  ok <- is.element(nc,max(nc)+c(-1,0))
  if (file.exists('readEALAT.txt')) file.remove('readEALAT.txt')
  writeLines(zz[(skip+1):length(zz)][ok],con='readEALAT.txt')
  if (verbose) print(zz[1:skip])
  x <- read.fwf('readEALAT.txt',width=width,as.is=TRUE)
  t <- as.Date(paste(x[[iyr]],x[[iyr+1]],x[[iyr+2]],sep='-'))

  ix <- switch(tolower(param),'t2m'=5,'precip'=5,'sd'=4,
               'tx'=6,'tn'=7)

  if (verbose) {
    print(ix)
    print(summary(x[[ix]]))
    print(as.character(x[[ix]][1:10]))
  }
  yz <- zoo(as.numeric(x[[ix]]),order.by=t)
  lat <- round(as.numeric(substr(h[3],37,39)) +
         as.numeric(substr(h[3],40,42))/60,4)
  minutes <- as.numeric(substr(h[4],51,54))/60
  if (is.na(minutes)) minutes <- 0
  lon <- round(as.numeric(substr(h[4],48,50)) + minutes,4)
  if (is.na(lon)) {print(h[4]); browser()}
         
  alt <- as.numeric(substr(h[5],32,35))

  unit <- switch(param, 't2m'='degC','precip'='mm/day','sd'='cm',
               'tx'='degC','tn'='degC')

  if (verbose) print(c(substr(h[1],7,nchar(h[1])),
                       param,unit,lon,lat,alt,
                       substr(h[2],32,nchar(h[2]))))
  if (verbose) {print(range(index(yz))); plot(yz)}
  y <- as.station(yz,param=param,unit=unit,
                  loc=substr(h[1],7,nchar(h[1])),
                  stid=as.numeric(substr(h[2],32,nchar(h[2]))),
                  lat=lat,lon=lon,alt=alt,src='AARI/EALAT')

  invisible(y)
}

## Files with oth3r format:
readothert2m <- function(param='t2m') {
  files <- list.files(path='data/EALAT_Russia',
                            full.names=TRUE)
  files <- files[c(48,50,58,63)]
  print(files)
  iyr=2
  ix <- switch(tolower(param),'t2m'=6,'tx'=7,'tn'=5)
  for (i in 1:length(files)) {
    fname <- files[i]
    print(fname)
    h <- readLines(fname,n=5)
    z <- read.fwf(fname,width=c(5,5,3,3,7,7,7),skip=5,as.is=TRUE)
    t <- as.Date(paste(z[[iyr]],z[[iyr+1]],z[[iyr+2]],sep='-'))
    yz <- zoo(as.numeric(z[[ix]]),order.by=t)
    lat <- round(as.numeric(substr(h[2],6,9)) +
      as.numeric(substr(h[2],10,12))/60,4)
    lon <- round(as.numeric(substr(h[3],6,9)) +
      as.numeric(substr(h[3],10,12))/60,4)
    alt <- NA
    
    unit <- switch(param, 't2m'='degC','precip'='mm/day','sd'='cm',
                   'tx'='degC','tn'='degC')
    y <- as.station(yz,param=param,unit=unit,
                    loc=substr(h[1],2,nchar(h[1])),
                    stid=as.numeric(z[[1]])[1],
                    lat=lat,lon=lon,alt=alt,src='AARI/EALAT')

    if (i==1) Y <- y else
              Y <- combine(Y,y)
  }
  invisible(Y)
}

## Extract the EALAT data and save as esd station object
EALAT2esd <- function(param='t2m',verbose=FALSE) {

  path=switch(tolower(param),
    't2m'='data/EALAT_Russia/Temp_dogn',
    'tx'='data/EALAT_Russia/Temp_dogn',
    'tn'='data/EALAT_Russia/Temp_dogn',
    'precip'='data/EALAT_Russia/Nedbor_Dogn',
    'sd'='data/EALAT_Russia/Sno_Dogn')

  print(path)
  files <- list.files(path=path,pattern='txt',full.names=TRUE)

  print(files)
  N <- length(files)

  data(geoborders)

  par(bty="n",xaxt="n",yaxt="n")
  plot(geoborders,col="grey",cex=0.5,xlab="",ylab="",
       xlim=c(-10,180),ylim=c(50,75),type='l')    

  for (i in 1:N) {
    y <- readEALAT(files[i],param=param,verbose=verbose)
    print(c(i,loc(y),lon(y),lat(y),alt(y),files[i]))
    points(attr(y,'longitude'),attr(y,'latitude'),col="red",pch=19)

    ## Remove entries with duplicate dates: 
    if (sum(duplicated(index(y)))>0) {
      print('Remove duplicates!')
      y <- subset(y,it=!duplicated(y))
      #browser()
    }
    
    if (i==1) Y <- y else
    if (length(y) > 3600) Y <- combine(Y,y)
  }
    
  if (sum(is.element(c('t2m','tn','tx'),tolower(param)))>0) {
    print('Pick up the other temperature data in seperate/different files')
    Y <- combine(Y,readothert2m(param))
  }
  
  save(file=paste('ealat.',param,'.rda',sep=''),Y)
  invisible(Y)
}

t2m <- EALAT2esd()
map(t2m,FUN='mean',cex=-2)

pr <- EALAT2esd(param='precip')
map(pr,FUN='mean',cex=-2)

#No good data: Tx
Tx <- EALAT2esd(param='Tx')
map(Tx,FUN='mean',cex=-2)

Tn <- EALAT2esd(param='Tn')
map(Tn,FUN='mean',cex=-2)

sd <- EALAT2esd(param='sd')
map(sd,FUN='mean',cex=-2)

ss <- select.station(param='t2m',lat=c(65,90),lon=c(0,100),src='ecad')
t2m.aaca <- station(ss)
save(file='t2m.aaca.rda',t2m.aaca)

### GHCND - too much missing data
#ss <- select.station(cntr='Russia',param='tmin',lat=c(65,90),lon=c(0,100))
#tn <- station(ss)
#ss <- select.station(cntr='Russia',param='tmax',lat=c(65,90),lon=c(0,100))
#tx <- station(ss)
#t2m <- 0.5*(tn + tx)
#t2m <- attrcp(tx,t2m)
#attr(t2m,'variable')[] <- 't2m'
#attr(t2m,'longname')[] <- 'daily temperature'
#attr(t2m,'info') <- 't2m <- 0.5*(tn + tx)'
#class(t2m)
#save(file='aaca.GHCND.t2m.rda',tn,tx,t2m)
