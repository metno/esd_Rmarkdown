## author: K Parding
## cyclone identification & tracking

library(esd)
force <- FALSE
path.era <- "/vol/fou/klima/rasmusb"

lon <- NULL
lat <- c(0,90)
for (yr in seq(1979,2012,1)) {
  if (!file.exists(paste("cyclones.NH.",yr,".rda",sep="")) | force) {
    t1 <- Sys.time()
    print(yr)
    fname <- file.path(path.era,"ERAINT-full-slp00_1979-2012.nc")
    Z1 <- retrieve(fname,type="ncdf",it=c(yr,yr),lon=lon,lat=lat)
    index(Z1) <- as.POSIXct(paste(index(Z1)," 00:00:00",sep=""))
    fname <- file.path(path.era,"ERAINT-full-slp06_1979-2012.nc")
    Z2 <- retrieve(fname,type="ncdf",it=c(yr,yr),lon=lon,lat=lat)
    index(Z2) <- as.POSIXct(paste(index(Z2)," 06:00:00",sep=""))
    fname <- file.path(path.era,"ERAINT-full-slp12_1979-2012.nc")
    Z3 <- retrieve(fname,type="ncdf",it=c(yr,yr),lon=lon,lat=lat)
    index(Z3) <- as.POSIXct(paste(index(Z3)," 12:00:00",sep=""))
    fname <- file.path(path.era,"ERAINT-full-slp18_1979-2012.nc")
    Z4 <- retrieve(fname,type="ncdf",it=c(yr,yr),lon=lon,lat=lat)
    index(Z4) <- as.POSIXct(paste(index(Z4)," 18:00:00",sep=""))
    Z <- rbind(Z1,Z2,Z3,Z4)
    Z <- as.field(Z,lon(Z1),lat(Z1),attr(Z1,"variable"),attr(Z1,"unit"))
    Z <- attrcp(Z1,Z)
    rm("Z1","Z2","Z3","Z4"); gc(reset=TRUE)
    t2 <- Sys.time()
    print(paste('retrieve took',round(as.numeric(t2-t1,units="secs")),'s'))
    t3 <- Sys.time()
    cci.yr <- CCI(Z,m=14,it=NULL,is=NULL,cyclones=TRUE,accuracy=NULL,
           label=NULL,lplot=FALSE,verbose=TRUE,
           fname=paste("cyclones.NH.",yr,".rda",sep=""))
    t4 <- Sys.time()
    print(paste('CCI of the field took',
            round(as.numeric(t4-t3,units="secs")),'s'))
    rm("Z","cci.yr","t1","t2","t3","t4"); gc(reset=TRUE)
  }
}

X0 <- NULL
files <- list.files(pattern="cyclones.NH.[0-9]{4}.rda")
for (f in files) {
  load(f)
  yr <- paste(unlist(regmatches(f,gregexpr("[0-9]",f))),collapse="")
  fname <- paste("cyclones.NH.",yr,".tracked.rda",sep="")
  if (!file.exists(fname) | force) {
    Xf <- X[,!colnames(X) %in% c("trajectory","tracklength",
             "trackcount","distance","timestep")]
    Xf <- attrcp(X,Xf)
    class(Xf) <- class(X)
    X <- Track.events(Xf,x0=X0,verbose=TRUE)
    save(file=fname,X)
  }
  X0 <- X
}

if (!file.exists("cyclone.NH.rda")) {
  files <- list.files(pattern="cyclones.NH.[0-9]{4}.tracked")
  Y <- NULL
  for (f in files) {
    print(f)
    load(f)
    if (is.null(Y)) {
      Y <- X
    } else {
      Y <- merge(Y,X,all=TRUE)
    }
  }
  Y <- attrcp(X,Y)
  class(Y) <- class(X)
  save(file="cyclones.NH.rda",Y)
} else {
  load("cyclones.NH.rda")
}
  

