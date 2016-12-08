# R.E. Benestad (rasmus.benestad@physics.org)
# R-script for reading CDCN daily precipitation data
# NCAR, Mesa Lab, Boulder March 31, 2011

# Look for similarity between exponential distribution and 24-hr precip world wide
# Look for systematic relationship between mu, meanP, and meanT, hence a dependence of shape of p.d.f. on
# local mean T(2m) and precip.



readGDCN <- function(name,param="precip", metadata="gdcn.meta.rda",
                     lon.rng=c(-180,360),lat.rng=c(-90,90),plot=FALSE) {

  if (is.character(metadata)) load(metadata)
  type <- switch(tolower(param),
                 "precip"="PRCP","t2m"="TMAX",
                 "prcp2"="PRCP",
                 "tmax"="TMAX","tmin"="TMIN")
  unit <- switch(tolower(param),
                 "precip"="mm/day","t2m"="deg C","tmax"="deg C","tmin"="deg C")
  # Reads the ASCII files with daily GDCN data:
  cnames <- c(paste("day.",1:31,sep=""),paste("flags.day.",1:31,sep=""))
  reshufle <- rep(0,62); ii <- 0
  for (i in seq(1,62,by=2)) {
    ii <- ii + 1
    reshufle[i:(i+1)] <- seq(ii,62,by=31)
  }
  cnames <- cnames[reshufle]

  ccn <- as.numeric(substr(name,nchar(name)-14,nchar(name)-12)) # Country code from file name
  x <- read.fwf(name,widths=c(3,8,4,2,4,rep(c(5,2),31)),
                col.names=c("country","stnr","year","month","type",cnames))
  imatch <- is.element(gdcn.meta$stnr,x$stnr[1]) &
            is.element(gdcn.meta$country.code,x$country[1]) &
            (gdcn.meta$lon >= min(lon.rng)) & (gdcn.meta$lon <= max(lon.rng)) &
            (gdcn.meta$lat >= min(lat.rng)) & (gdcn.meta$lat <= max(lat.rng))
  #str(x); print(cnames)
  ipick <- is.element(x$type,type)
  year <- sort(rep(x$year[ipick],31))
  month <- rep(x$month[ipick],each=31)
  day <- rep(1:31,sum(ipick))
  
  if ( (sum(ipick)>0) & (sum(imatch)==1) ) {

    location <- gdcn.meta$location[imatch]
    if (is.na(location)) location <- paste('stid',gdcn.meta$stnr[imatch],sep='_')
    while (substr(location,nchar(location),nchar(location))==" ")
      location <- substr(location,1,nchar(location)-1)
    dat <- c(t(as.matrix(x[ipick,seq(6,67,by=2)])))*0.1
    qua <- as.matrix(x[ipick,seq(7,68,by=2)])
    dat[dat < -99] <- NA
    OK <- is.finite(dat)
    dat <- dat[OK]
    qua <- qua[OK]
    t <- as.Date(paste(year[OK],month[OK],day[OK],sep='-'))
    ok <- is.finite(dat) & is.finite(t)
    #print(t)
    url <- "http://www.ncdc.noaa.gov/oa/climate/research/gdcn/gdcn.html"
    y <- as.station(zoo(dat[ok],order.by=t[ok]),loc=location,
                    param=param,unit=unit,lon=gdcn.meta$lon[imatch],
                    lat=gdcn.meta$lat[imatch],alt= gdcn.meta$alt[imatch],
                    cntr=gdcn.meta$country[imatch],longname=NA,
                    stid=x$stnr[1],quality=c(qua)[ok],
                    src="GDCN",url=url,method="observation")
    #str(y)
  } else {
    print(paste("station_id",x$stnr[1]," and element",type," search file",
                sum(ipick),"lines with matching type; search inventory",
                sum(imatch),"matche(s)"))
    #print(table(gdcn.meta$country[imatch]))
    #print(table(x$type))
    #print(table(round(gdcn.meta$lat[imatch]),round(gdcn.meta$lon[imatch])))
    #print(rbind(gdcn.meta$location[imatch],gdcn.meta$alt[imatch]))
    #print(summary(c(as.matrix(x[ipick,seq(6,67,by=2)])*0.1)))
    if (plot) {
      require(esd)
      par(bty="n")
      plot(gdcn.meta$lon[imatch],gdcn.meta$lat[imatch],pch=19,col="red",
           xlab="",ylab="",xlim=c(-180,180),ylim=c(-90,90))
      data(geoborders)
      lines(geoborders)
      text(gdcn.meta$lon[imatch],gdcn.meta$lat[imatch],
           gdcn.meta$location[imatch],cex=0.7)
    }
    
    y <- NULL
  }
  invisible(y)
}





