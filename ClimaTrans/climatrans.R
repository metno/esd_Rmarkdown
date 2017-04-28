## ----eval=FALSE----------------------------------------------------------
## library(rmarkdown); render('climatrans.Rmd',pdf_document())

## ---- eval=FALSE---------------------------------------------------------
## library(knitr)
## purl('~/git/esd_Rmarkdown/ClimaTrans/climatrans.Rmd', output='~/git/esd_Rmarkdown/ClimaTrans/climatrans.R')

## ----echo=TRUE-----------------------------------------------------------
library(esd)
library(ncdf4)
setwd('~/R')
#rm(list=ls())
## Information about the system and session

## ------------------------------------------------------------------------
tas2es <- function(x,mask=TRUE,land=TRUE,season=5:10,FUN='min') {
  if (mask) x <- mask(x,land=land)
  es <- subset(C.C.eq(x),it=month.abb[season])
  es <- aggregate(es,year,FUN=FUN)
  index(es) <- year(es)
  invisible(es)
}

## ------------------------------------------------------------------------
PrexpPr <- function(mu,fw,x0=10) {
  Pr <- zoo(coredata(fw)*exp(-x0/coredata(mu)),order.by=year(fw))

  Pr.mu <- zoo(mean(coredata(fw))*exp(-x0/coredata(mu)),order.by=year(fw))
  Pr.fw <- zoo(coredata(fw)*exp(-x0/mean(coredata(mu))),order.by=year(fw))
  attr(Pr,'prob.f(mu)') <- Pr.mu
  attr(Pr,'prob.f(fw)') <- Pr.fw
  Pr <- attrcp(mu,Pr)
  attr(Pr,'variable') <- 'Pr'
  attr(Pr,'unit') <- 'probability'
  attr(Pr,'longname') <- paste('Probability of exceeding',x0)
  class(Pr) <- class(mu)
  invisible(Pr)
}

## ------------------------------------------------------------------------
RMSE <- function(x,y) sqrt(sum((x-y)^2))/length(x)

precipcmip <- function(predictor='data/ERAINT/era-int-precip-mon.nc',
                    path="CMIP5.monthly/",rcp='rcp45',pattern="pr_Amon_ens_",
                    type='ncdf',it=c(1950,2015),plot=TRUE,verbose=FALSE,
                    lon=c(70,90),lat=c(7,30),select=NULL,xlim=NULL,fname=NULL) {
  
  if (is.null(fname)) fname <- paste('precipcmip.',paste(lon,lat,sep='.',collapse='-'),'.rda',sep='')
  
  print(paste('precipcmip: results (will be) stored in',fname))
  
  if (file.exists(fname)) {
    load(fname)
    return(results)
  }

  ## Diagnose the common EOFs and effect of bias-adjustment

  if (verbose) print(paste('getcmip, it=',it))

  if (verbose) print('predictor')
  if (is.character(predictor))
    ## Unit: "mm"
    rea <- retrieve(ncfile=predictor,lon=lon,lat=lat,
                    type=type,verbose=verbose) else
  if (inherits(predictor,'field'))
    rea <- subset(predictor,is=list(lon=lon,lat=lat))
  if (!is.null(it)) {
    rea <- subset(rea,it=it)
    it <- c(min(year(rea)),max(year(rea)))
  }
  
  ## The reanalysis
  if (verbose) print('aggregate the reanalysis')
  Rea <- aggregate.area(rea,FUN='mean')
  X <- aggregate(Rea,month,FUN='mean')
  Xt <- annual(Rea,FUN='sum')

  if (plot) plot(X,lwd=3,col='black',main='Climatology',ylim=c(0,2*round(max(X))),
                 ylab=ylab(X),xlab='Calendar month',map.type='rectangle')
  
  # Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
  N <- length(ncfiles)

  if (is.null(select)) select <- 1:N else
                       N <- length(select)
  if (verbose) {print('GCMs:'); print(path); print(ncfiles[select])}

  Z <- matrix(rep(NA,N*12),N,12)
  Zt <- matrix(rep(NA,N*201),N,201)
  gcmnm <- rep("",N); rmse <- rep(NA,N); r <- rep(NA,N)

  ## Set up a list environment to keep all the results
  if (verbose) print("loop...") 
  for (i in 1:N) {
    if (verbose) print(ncfiles[select[i]])

    ## unit: "kg m-2 s-1"
    gcm <- retrieve(ncfile = ncfiles[select[i]],type=type,
                          lon=range(lon(rea)),
                          lat=range(lat(rea)),verbose=verbose)
    GCM <- aggregate.area(gcm,FUN='mean')
    ## Annual total precipitation
    zt <- annual(GCM,FUN='sum')
    i1 <- is.element(year(zt),1900:2100)
    i2 <- is.element(1900:2100,year(zt))
    Zt[i,i2] <- coredata(zt)[i1]

    if (!is.null(it)) {
      if (verbose) print('Extract some months ot a time period')
      if (verbose) print(it)
      gcm <- subset(gcm,it=it)
    }
    ## The mean annual cycle
    ## The units of the GCMs are mm/day - multiply by 30
    z <- 30*aggregate(GCM,month,FUN='mean')
    Z[i,] <- coredata(z)
    #gcmnm[i] <- attr(gcm,'model_id')
    gcmnm.i <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-r")
    gcmnm.i <- gsub(' ','',gcmnm.i)
    gcmnm.i <- gsub('-','.',gcmnm.i)
    rmse[i] <- round(RMSE(coredata(X),coredata(z)),2)
    r[i] <- round(cor(coredata(X),coredata(z)),2)
    gcmnm[i] <- gcmnm.i
    print(paste(gcmnm.i,unit(GCM),unit(Rea),rmse[i],r[i]))
    if (plot) lines(z,col=rgb(i/N,0,0.5+0.5*r[i],0.1),lwd=3)
  }
  attr(Z,'GCMs') <- gcmnm
   if (plot) {
     lines(X,col='black',lwd=3)
     grid()
     figlab(paste(start(Rea),end(Rea)))
   }
  results <- list(reanalysis=X,GCMs=z,
                  obs=Xt,gcm=zoo(t(Zt),order.by=1900:2100),rmse=rmse,cor=r)
  dev.copy2pdf(file=paste('precipcmip.',
                 paste(lon,lat,sep='.',collapse='-'),'.pdf',sep=''))
  save(file=fname,results)
  return(results)
}

## ---- eval=FALSE---------------------------------------------------------
## hotsummerdays <- function (x, y = NULL, dse = NULL, it = "jja", threshold = 30,
##     verbose = FALSE, plot = TRUE, nmin = 90, new = TRUE, ...)
## {
##     if (verbose)
##         print("mildwinterdays")
##     stopifnot(inherits(x, "station"))
##     if (is.null(y))
##         y <- x
##     djf <- subset(x, it = it)
##     djfy <- subset(y, it = it)
##     nwd1 <- annual(djfy, FUN = "count", threshold = threshold,
##         nmin = nmin)
##     mwd1 <- annual(djf, FUN = "mean", nmin = nmin)
##     cal <- data.frame(x = c(coredata(mwd1)), y = c(coredata(nwd1)))
##     dfit <- glm(y ~ x, family = "poisson", data = cal)
##     if (plot) {
##         par(bty = "n")
##         plot(cal, pch = 19, ylim = c(0, 90), xlab = expression(paste("mean temperature ",
##             (degree * C))), ylab = "number of hot days", main = loc(x))
##         pre <- data.frame(x = seq(min(cal$x, na.rm = TRUE) -
##             1, max(cal$x, na.rm = TRUE) + 5, by = 0.1))
##         lines(pre$x, exp(predict(dfit, newdata = pre)), col = rgb(1,
##             0, 0, 0.3), lwd = 3)
##         djf.sd <- sd(coredata(djf), na.rm = TRUE)
##         qqnorm(coredata(djf))
##         qqline(coredata(coredata(djf)), col = "red")
##         grid()
##     }
##     if (is.null(dse))
##         dse <- DSensemble.t2m(x, biascorrect = TRUE, verbose = verbose,
##             plot = plot)
##     djf.dse <- subset(dse, it = "djf")
##     index(djf.dse) <- year(djf.dse)
##     ovl <- window(djf.dse, start = year(start(x)), end = year(end(x)))
##     djf.dse <- djf.dse - mean(coredata(ovl), na.rm = TRUE) +
##         mean(coredata(mwd1), na.rm = TRUE)
##     q1 <- data.frame(x = apply(coredata(djf.dse), 1, quantile,
##         probs = 0.05, na.rm = TRUE))
##     q2 <- data.frame(x = apply(coredata(djf.dse), 1, quantile,
##         probs = 0.95, na.rm = TRUE))
##     qm <- data.frame(x = apply(coredata(djf.dse), 1, mean, na.rm = TRUE))
##     obs <- data.frame(x = coredata(mwd1))
##     t <- year(index(djf.dse))
##     preq1 <- exp(predict(dfit, newdata = q1))
##     preq1[preq1 > 90] <- NA
##     tr1 <- predict(lm(preq1 ~ t + I(t^2) + I(t^3) + I(t^4) +
##         I(t^5)))
##     tr1[!is.finite(preq1)] <- NA
##     preq2 <- exp(predict(dfit, newdata = q2))
##     preq2[preq2 > 90] <- NA
##     tr2 <- predict(lm(preq2 ~ t + I(t^2) + I(t^3) + I(t^4) +
##         I(t^5)))
##     tr2[!is.finite(preq2)] <- NA
##     prem <- exp(predict(dfit, newdata = qm))
##     prem[prem > 90] <- NA
##     tr3 <- predict(lm(prem ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5)))
##     tr3[!is.finite(prem)] <- NA
##     Nwd <- zoo(cbind(preq1, preq2, prem, tr1, tr2, tr3), order.by = t)
##     nwd.pre <- zoo(exp(predict(dfit, newdata = obs)), order.by = year(mwd1))
##     if (plot) {
##         par(bty = "n")
##         plot(zoo(djf.dse, order.by = year(djf.dse)), plot.type = "single",
##             col = rgb(0.5, 0.5, 0.5, 0.2), ylab = expression(paste("mean temperature",
##                 (degree * C))), xlab = "", main = loc(x))
##         points(mwd1, pch = 19)
##         grid()
##         par(bty = "n")
##         plot(Nwd, plot.type = "single", lwd = 5, main = loc(x),
##             ylim = c(0, 90), xlab = "", ylab = paste("number of hot days: T(2m) > ",
##                 threshold, unit(x)), col = c(rgb(0.5, 0.5, 0.7,
##                 0.5), rgb(0.8, 0.5, 0.5, 0.5), rgb(0.8, 0.5,
##                 0.8, 0.5), rgb(0.3, 0.3, 0.6, 0.5), rgb(0.6,
##                 0.3, 0.3, 0.5), rgb(0.6, 0.3, 0.6, 0.5)), ...)
##         grid()
##         points(nwd1, pch = 19)
##         lines(nwd.pre, col = rgb(0.5, 0.5, 0.5, 0.5))
##     }
##     Nwd <- attrcp(x, Nwd)
##     attr(Nwd, "unit") <- "days"
##     attr(Nwd, "info") <- paste("number of hot days: t2m > ",
##         threshold)
##     attr(Nwd, "observation") <- nwd1
##     attr(Nwd, "nwd.pre") <- nwd.pre
##     index(Nwd) <- t
##     class(Nwd) <- c("nevents", "zoo")
##     invisible(Nwd)
## }

## ------------------------------------------------------------------------
showceofvar <-  function(x,N,i,gcmnm.i,add=FALSE,xlim=NULL,ylim=NULL,verbose=FALSE,...) {
  n <- length(x$mean.diff)
  j <- 1:n
  col <- rgb(j/n,abs(sin(pi*j/n)),(1-j/n),0.3)
  if (is.null(xlim)) xlim <- c(-1,1)
  if (is.null(ylim)) ylim <- c(1,N)
  
  if (!add) {
    par(bty="n")
    par0 <- par()
    plot(xlim,ylim,type="n",
         ylab="GCM",xlab=expression(1- sigma[p*r*e]/sigma[r*e*f]),
         main=paste("Diagnostics: common EOFs",attr(x,'variable')),
         xlim=xlim,ylim=ylim,col="grey")
    grid()
    #legend(xlim[1],ylim[2],c("same sign","different sign"),
    #       pch=c(19,21),bty="n",col="grey")
    par(xpd=TRUE)
    text(xlim[2],ylim[2],'EOF #',col='grey40',cex=0.8,pos=3)

    par(new=TRUE,fig=c(0.85,0.95,0.70,0.85),mar=c(0,3,0,0),
        cex.axis=0.7,yaxt="s",xaxt="n",las=1)
     colbar <- rbind(1:n,1:n)
    image(1:2,1:n,colbar,col=col)
    par(fig=par0$fig,mar=par0$mar,cex.axis=par0$cex.axis,
        yaxt=par0$yaxt,xaxt=par0$xaxt,las=par0$las,new=TRUE)
    plot(xlim,ylim,type="n",xlim=xlim,ylim=ylim,
         xlab='',ylab='',main='',sub='')
    par(par0$new)
  }
  cex <- abs(0.1 - x$mean.diff)/0.1
  pch <- rep(19,n); pch[cex < 0] <- 21
  cex <- abs(cex); cex[cex > 2] <- 2
  if (verbose) {
     print('Mean difference:');print(x$mean.diff)
     print('Ration of standard deviation');print(x$sd.ratio)
     print('Size');print(cex)
     print('col');print(col)
     #points(x$mean.diff,1-x$sd.ratio,pch=pch,col='grey75',cex=1)
  }
  
  points(1-x$sd.ratio,rep(i,n),pch=pch,col=col,cex=cex)
  #print(paste(' > ',i,i,gcmnm.i[i]))
  text(xlim[1],i,gcmnm.i[i],pos=4,cex=0.6)
}


biasdiag <- function(predictor,param='slp',it=NULL,
                    path="CMIP5.monthly/",pattern="psl_Amon_ens_",type='ncdf4',
                    period=c(1950,2015),rcp='rcp45',verbose=TRUE,
                    lon=c(0,100),lat=c(65,90),select=NULL,xlim=NULL) {

  ## Diagnose the common EOFs and effect of bias-adjustment

  if (verbose) print(paste('biasdiag, it=',it))

  # Ensemble GCMs
  path <- file.path(path,rcp,fsep = .Platform$file.sep)
  ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
  N <- length(ncfiles)

  if (is.null(select)) select <- 1:N else
                       N <- length(select)
  if (verbose) {print('GCMs:'); print(path); print(ncfiles[select])}
  
  nceofs <- 6 # Number of leading EOFs to evalue
  years <- 1900:2100
  m <- length(years)
  months <- rep(month(nceofs),m)
  X <- matrix(rep(NA,N*m*nceofs),N,m*nceofs)
  dim(X) <- c(nceofs,N,m)
  gcmnm.i <- rep("",N)
  scorestats <- matrix(rep(NA,N*9),N,9)
  colnames(scorestats) <- c("1-r.xval","mean.diff","sd.ratio","autocorr.ratio",
                            "res.trend","res.K-S","res.ar1",'amplitude.ration',
                            "1-R2")

  t <- as.Date(paste(years,months,'01',sep='-'))

  cols <- rgb(seq(1,0,length=100),rep(0,100),seq(0,1,length=100),0.15)

  flog <- file("DSensemble.ceof-log.txt","at")
  
  ## Set up a list environment to keep all the results
  results <- list() 
  if (!file.exists('diag.ceof.slp.rda')) {
    if (verbose) print("loop...") 
    for (i in 1:N) {
    
      if (verbose) print(ncfiles[select[i]])
      gcm <- retrieve(ncfile = ncfiles[select[i]],type=type,
                            lon=range(lon(predictor)),
                           lat=range(lat(predictor)),verbose=verbose)
      if (!is.null(it)) {
        if (verbose) print('Extract some months ot a time period')
        if (verbose) print(it)
        gcm <- subset(gcm,it=it)
        gcm <- subset(gcm,it=range(year(predictor)))
      }
      gcmnm <- paste(attr(gcm,'model_id'),attr(gcm,'realization'),sep="-r")
      gcmnm <- gsub(' ','',gcmnm)
      gcmnm <- gsub('-','.',gcmnm)
      gcmnm.i[i] <- gcmnm
      print(gcmnm.i[i])
      GCM <- aggregate(gcm,year,FUN='mean') 
      SLPGCM <- combine(predictor,GCM)
      Z <- try(EOF(SLPGCM))
      diag <- diagnose(Z)
      showceofvar(diag,N,i,gcmnm.i,add=(i!=1),xlim=xlim,verbose=verbose)
      eval(parse(text=paste('results$member.',i,'.',gcmnm.i[i],' <- diag',sep='')))
    } 
    save(file='diag.ceof.slp.rda',results,gcmnm.i)
  } else {
    load('diag.ceof.slp.rda')
    N <- length(gcmnm.i)
    for (i in 1:N) {
      #print(paste(' < ', i,gcmnm.i[i]))
      showceofvar(results[[i]],N,i,gcmnm.i,add=(i!=1),xlim=xlim,verbose=verbose)
    }
  }
  dev.copy2pdf(file=paste('diag.ceof.slp.india.pdf',sep=''))
  return(results)
}


## ----eval=FALSE----------------------------------------------------------
## ## Precipitation
## ## time origin: minutes since 1901-01-01 00:00
## library(ncdf4)
## ncid <- nc_open('~/Dropbox/Public/ClimaTrans/delhi.nc')
## delhi <- ncvar_get(ncid,varid='p')
## time <- ncvar_get(ncid,varid='time')
## nc_close(ncid)
## delhi <- zoo(x=delhi,order.by=as.Date(time/(24*60),origin='1901-01-01'))
## delhi <- as.station(delhi,loc='Delhi',param='precip',unit='mm/day',lon=77.2,lat=28.6,
##                     cntr='India',longname='precipitation',
##                     info='extracted rainfall data from IMD gridded daily data',
##                     ref='Madhusoodanan M.S. ("Madhusoodanan M.S." <madhusoodanan@gmail.com>)')
## 
## ncid <- nc_open('~/Dropbox/Public/ClimaTrans/bombay.nc')
## bombay <- ncvar_get(ncid,varid='p')
## time <- ncvar_get(ncid,varid='time')
## nc_close(ncid)
## bombay <- zoo(x=bombay,order.by=as.Date(time/(24*60),origin='1901-01-01'))
## bombay <- as.station(bombay,loc='Bombay',param='precip',unit='mm/day',lon=72.9,lat=19.0,
##                 cntr='India',longname='precipitation',
##                 info='extracted rainfall data from IMD gridded daily data',
##                 ref='Madhusoodanan M.S. ("Madhusoodanan M.S." <madhusoodanan@gmail.com>)')
## 
## ncid <- nc_open('~/Dropbox/Public/ClimaTrans/bangalore.nc')
## bangalore <- ncvar_get(ncid,varid='p')
## time <- ncvar_get(ncid,varid='time')
## nc_close(ncid)
## 
## bangalore <- zoo(x=bangalore,order.by=as.Date(time/(24*60),origin='1901-01-01'))
## bangalore <- as.station(bangalore,loc='Bangalore',param='precip',unit='mm/day',
##                 lon=77.6,lat=13.0,cntr='India',longname='precipitation',
##                 info='extracted rainfall data from IMD gridded daily data',
##                 ref='Madhusoodanan M.S. ("Madhusoodanan M.S." <madhusoodanan@gmail.com>)')
## 
## ## Combine the individual stations into a group of stations
## climatrans.pr <- combine(delhi,bombay,bangalore)
## ## Save the results for future use
## save(file='climatrans.pr.rda',climatrans.pr)
## 
## ## Maximum temperature
## ## time origin: minutes since 1969-01-01 00:00
## ncid <- nc_open('~/Dropbox/Public/ClimaTrans/delhi_tmax.nc')
## delhi <- ncvar_get(ncid,varid='temp')
## time <- ncvar_get(ncid,varid='time')
## nc_close(ncid)
## delhi <- zoo(x=delhi,order.by=as.Date(time/(24*60),origin='1969-01-01'))
## tx.delhi <- as.station(delhi,loc='Delhi',param='tmax',unit='degC',lon=77.2,lat=28.6,
##                cntr='India',longname='daily maximum temperature',
##                info='extracted from IMD gridded daily data',
##                ref='Madhusoodanan M.S. ("Madhusoodanan M.S." <madhusoodanan@gmail.com>)')
## 
## ncid <- nc_open('~/Dropbox/Public/ClimaTrans/bombay_tmax.nc')
## bombay <- ncvar_get(ncid,varid='temp')
## time <- ncvar_get(ncid,varid='time')
## nc_close(ncid)
## bombay <- zoo(x=bombay,order.by=as.Date(time/(24*60),origin='1969-01-01'))
## tx.bombay <- as.station(bombay,loc='Bombay',param='tmax',unit='degC',lon=72.9,lat=19.0,
##                 cntr='India',longname='daily maximum temperature',
##                 info='extracted from IMD gridded daily data',
##                 ref='Madhusoodanan M.S. ("Madhusoodanan M.S." <madhusoodanan@gmail.com>)')
## 
## ncid <- nc_open('~/Dropbox/Public/ClimaTrans/bangalore_tmax.nc')
## bangalore <- ncvar_get(ncid,varid='temp')
## time <- ncvar_get(ncid,varid='time')
## nc_close(ncid)
## bangalore <- zoo(x=bangalore,order.by=as.Date(time/(24*60),origin='1969-01-01'))
## tx.bangalore <- as.station(bangalore,loc='Bangalore',param='tmax',unit='degC',
##                 lon=77.6,lat=13.0,cntr='India',longname='daily maximum temperature',
##                 info='extracted from IMD gridded daily data',
##                 ref='Madhusoodanan M.S. ("Madhusoodanan M.S." <madhusoodanan@gmail.com>)')
## 
## ## Combine the individual stations into a group of stations
## climatrans.tx <- combine(tx.delhi,tx.bombay,tx.bangalore)
## ## Save the results for future use
## save(file='climatrans.tx.rda',climatrans.tx)

## ----echo=TRUE-----------------------------------------------------------
## Load precipitation data extracted from the gridded IMD data set:
load('climatrans.pr.rda')
## Display the structure of the data as
str(climatrans.pr)
plot(climatrans.pr,new=FALSE)

## ---- fig.width=6, fig.height=6------------------------------------------
## Wet-day frequency
plot(aggregate(climatrans.pr,month,'wetfreq'),new=FALSE)
grid()

## ---- fig.width=6, fig.height=6------------------------------------------
## Wet-day mean
plot(aggregate(climatrans.pr,month,'wetmean'),new=FALSE)
grid()

## ---- fig.width=6, fig.height=6------------------------------------------
## Monthly precipitation totals
plot(aggregate(climatrans.pr,month,'sum'),new=FALSE)
grid()

## ------------------------------------------------------------------------
wq95 <- function(x) {x <- x[is.finite(x)]; x <- x[x >= 1]; wq95 <- quantile(x,probs=0.95); wq95}
season<-5:10
fw <- aggregate(subset(climatrans.pr,it=month.abb[season]),year,'wetfreq')
mu <- aggregate(subset(climatrans.pr,it=month.abb[season]),year,'wetmean')
nd <- aggregate(subset(climatrans.pr,it=month.abb[season]),year,'nv')
ndhr100 <- aggregate(subset(climatrans.pr,it=month.abb[season]),year,'count',threshold=100)
ndhr50 <- aggregate(subset(climatrans.pr,it=month.abb[season]),year,'count',threshold=50)
q95 <- aggregate(subset(climatrans.pr,it=month.abb[season]),year,'wq95')

## ------------------------------------------------------------------------
## check the relationship between wet-day 95-percentile and wet-day mean
plot(-log(0.05)*zoo(mu[,2]),zoo(q95[,2]),pch=19,xlim=c(0,120),ylim=c(0,120),
     main='Wet-day mean v.s. 95th wet percentile',
     xlab=expression(paste(-ln(1-0.95)*mu,' (mm/day)')), ylab=expression(paste(q[95],' (mm/day)')))
points(-log(0.05)*zoo(mu[,1]),zoo(q95[,1]),pch=19,col=rgb(1,0,0,0.2))
points(-log(0.05)*zoo(mu[,3]),zoo(q95[,3]),pch=19,col=rgb(0,0,1,0.2))
grid()
legend(0,120,loc(climatrans.pr),pch=19,col=c(rgb(1,0,0,0.2),'black',rgb(0,0,1,0.2)),bty='n')

## ------------------------------------------------------------------------
z <- climatrans.pr[,3]; z <- z[z > 1]
hist(z,breaks=seq(0,500,by=2),xlim=c(0,100),freq=FALSE)
lines(seq(0,250,by=1),dexp(seq(0,250,by=1),rate=1/mean(z)),lwd=3,col='red')
grid()

## ---- fig.width=6, fig.height=6------------------------------------------
## Wet-day frequency
plot(fw,new=FALSE)
for(i in 1:3) lines(trend(subset(fw,is=i)))
grid()
print(as.numeric(trend(fw,result='coef'))) 
print(as.numeric(trend(fw,result='pval'))) 

## ---- fig.width=6, fig.height=6------------------------------------------
## Wet-day mean
plot(mu,new=FALSE)
for(i in 1:3) lines(trend(subset(mu,is=i)))
grid()
print(as.numeric(trend(mu,result='coef'))) 
print(as.numeric(trend(mu,result='pval')))

## ---- fig.width=6, fig.height=6------------------------------------------
## Probability of exceedin 50 mm/day
Pr <- PrexpPr(mu,fw,x0=50) 
plot(Pr,new=FALSE,errorbar=FALSE)
for(i in 1:3) lines(trend(subset(Pr,is=i)))
grid()
print(as.numeric(trend(Pr,result='coef')))
print(as.numeric(trend(Pr,result='pval')))

## ---- fig.width=6, fig.height=6------------------------------------------
## Aggregate to annual totals
pt <- annual(climatrans.pr,'sum')
attr(pt,'unit') <- 'mm/season'
plot(pt,new=FALSE)
for(i in 1:3) lines(trend(subset(pt,is=i)))
grid()
print(as.numeric(trend(pt,result='pval')))
print(as.numeric(trend(pt,result='pval')))

## ---- fig.width=6, fig.height=6------------------------------------------
plot(ndhr100,new=FALSE)
grid()

## ---- fig.width=6, fig.height=6------------------------------------------
plot(ndhr50,new=FALSE)
grid()

## ------------------------------------------------------------------------
print(summary(coredata(pt)))
print(summary(coredata(fw*nd)))
print(summary(coredata(mu)))

## ------------------------------------------------------------------------
## Retrieve the reanalysis from FTP
if (!file.exists("air.mon.mean.nc"))
  download.file('ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc',
              destfile = 'air.mon.mean.nc')
if (!file.exists("slp.mon.mean.nc"))
  download.file('ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/slp.mon.mean.nc',
              destfile = 'slp.mon.mean.nc')

## ------------------------------------------------------------------------
## Aggregate sub-season values to annual mean values
if (!exists("es")) {
  t2m <- retrieve('air.mon.mean.nc',lon=c(60,100),lat=c(0,23))
  es <- aggregate(tas2es(t2m,season=season),year,'mean')
}

## ------------------------------------------------------------------------
if (!exists("slp")) {
  slp <- retrieve('slp.mon.mean.nc',lon=c(60,100),lat=c(0,23))
  slp <- aggregate(subset(slp,it=month.abb[season]),year,'mean')
}


## ------------------------------------------------------------------------
## Correlation between the wet-day mean and the temperature-based predictor
y <- subset(mu,is=2)
corfield(y,es,colbar=list(pal='t2m',breaks=seq(-1,1,by=0.1)),new=FALSE)

## ------------------------------------------------------------------------
## Correlation between the wet-day frequency and the SLP
z <- subset(fw,is=2)
attr(z,'unit') <- 'correlation'
corfield(z,slp,colbar=list(pal='t2m',breaks=seq(-1,1,by=0.1)),new=FALSE)

## ------------------------------------------------------------------------
pca.mu <- PCA(mu)
pca.fw <-PCA(fw)
plot(pca.mu,new=FALSE)

## ------------------------------------------------------------------------
plot(pca.fw,new=FALSE)

## ------------------------------------------------------------------------
eof.es <- EOF(es)
eof.slp <- EOF(slp)
eof.t2m <- EOF(aggregate(subset(t2m,is=list(lon=c(65,80),lat=c(10,20))),year,'mean'))
index(eof.es) <- year(eof.es)
index(eof.slp) <- year(eof.slp)
index(eof.t2m) <- year(eof.t2m)

## ---- fig.width=6, fig.height=6------------------------------------------
plot(eof.es,new=FALSE)

## ---- fig.width=6, fig.height=6------------------------------------------
plot(eof.slp,new=FALSE)

## ---- fig.width=6, fig.height=6------------------------------------------
plot(eof.t2m,new=FALSE)

## ----fig.width=6, fig.height=6-------------------------------------------
## The function DS uses multiple regression by default
ds.mu <- DS(pca.mu,eof.es,eofs=1:20)
plot(ds.mu,colbar1=list(breaks=seq(0,1,by=0.1)),
           colbar2=list(breaks=seq(-500,500,by=25)),new=FALSE)

## ----pca-of-residual-mu--------------------------------------------------
## Examine the PCA of residuals for the cold months to see if there are any structures left
  mu.2 <- as.station(as.residual(ds.mu))
  pca.mu.2 <- PCA(mu.2)
  plot(pca.mu.2,new=FALSE)
  grid()
  str(mu.2)

## ----fig.width=6, fig.height=6-------------------------------------------
ndhr.cal <- data.frame(y = coredata(subset(subset(ndhr100,is=2),it=c(1949,2004))),
                       X = coredata(subset(eof.es,it=c(1949,2004))))
ds.ndhr <- glm(y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7, data=ndhr.cal,family="poisson")
plot(1949:2004,ndhr.cal$y,pch=19,ylab=expression(N(X > 100*mm)),xlab='',main=loc(ndhr100)[2])
grid()
lines(1949:2004,exp(predict(ds.ndhr)),col='red')

## ----fig.width=6, fig.height=6-------------------------------------------
## Delhi
ndhr50.1.cal <- data.frame(y = coredata(subset(subset(ndhr50,is=1),it=c(1949,2004))),
                       X = coredata(subset(eof.es,it=c(1949,2004))))
ds.ndhr50.1 <- glm(y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7, data=ndhr50.1.cal,family="poisson")
plot(1949:2004,ndhr50.1.cal$y,pch=19,ylab=expression(N(X > 50*mm)),xlab='',main=loc(ndhr50)[1])
grid()
lines(1949:2004,exp(predict(ds.ndhr50.1)),col='red')

# Bangalore
ndhr50.3.cal <- data.frame(y = coredata(subset(subset(ndhr50,is=3),it=c(1949,2004))),
                       X = coredata(subset(eof.es,it=c(1949,2004))))
ds.ndhr50.3 <- glm(y ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7, data=ndhr50.3.cal,family="poisson")
plot(1949:2004,ndhr50.3.cal$y,pch=19,ylab=expression(N(X > 50*mm)),xlab='',main=loc(ndhr50)[3])
grid()
lines(1949:2004,exp(predict(ds.ndhr50.3)),col='red')

## ------------------------------------------------------------------------
print(paste('Mean number=',round(mean(ndhr100[,2],na.rm=TRUE),1),
            'Variance=',round(var(ndhr100[,2],na.rm=TRUE),1),
            'should equal',round(mean(ndhr100[,2],na.rm=TRUE),1),'for Poisson process'))
h <- hist(coredata(ndhr100[,2]),freq=FALSE,col='grey',
     main=loc(ndhr100)[2])
lines(h$mids-0.5,dpois(x=h$mids-0.5,lambda=mean(ndhr100[,2],na.rm=TRUE)),lwd=2,col='red')
grid()

## ------------------------------------------------------------------------
# Delhi
print(paste('Mean number=',round(mean(ndhr50[,1],na.rm=TRUE),1),
            'Variance=',round(var(ndhr50[,1],na.rm=TRUE),1),
            'should equal',round(mean(ndhr50[,1],na.rm=TRUE),1),'for Poisson process'))
h <- hist(coredata(ndhr50[,1]),freq=FALSE,col='grey',
     main=loc(ndhr50)[1])
lines(h$mids-0.5,dpois(x=h$mids-0.5,lambda=mean(ndhr50[,1],na.rm=TRUE)),lwd=2,col='red')
grid()

# Bangalore
print(paste('Mean number=',round(mean(ndhr50[,3],na.rm=TRUE),1),
            'Variance=',round(var(ndhr50[,3],na.rm=TRUE),1),
            'should equal',round(mean(ndhr50[,3],na.rm=TRUE),1),'for Poisson process'))
h <- hist(coredata(ndhr50[,3]),freq=FALSE,col='grey',
     main=loc(ndhr50)[3])
lines(h$mids-0.5,dpois(x=h$mids-0.5,lambda=mean(ndhr50[,3],na.rm=TRUE)),lwd=2,col='red')
grid()

## ------------------------------------------------------------------------
data(bjornholt)
ndhr30 <- aggregate(bjornholt,year,'count',threshold=30)
print(paste('Mean number=',round(mean(ndhr30,na.rm=TRUE),1),'Variance=',round(var(ndhr30,na.rm=TRUE),1),
            'should equal',round(mean(ndhr30,na.rm=TRUE),1),'for Poisson process'))
h <- hist(coredata(ndhr30),freq=FALSE,col='grey',
     main=loc(ndhr30))
lines(h$mids-0.5,dpois(x=h$mids,lambda=mean(ndhr30,na.rm=TRUE)),lwd=2,col='red')
grid()

## ---- fig.height=8,fig.width=8-------------------------------------------
ds.fw <- DS(pca.fw,eof.slp,eofs=1:20)
plot(ds.fw,new=FALSE)

## ----pca-of-residual-fw, fig.height=8,fig.width=8------------------------
## Examine the PCA of residuals for the cold monthsto see if there are any structures left
  fw.2 <- as.station(as.residual(ds.fw))
  pca.fw.2 <- PCA(fw.2)
  plot(pca.fw.2,new=FALSE)
  grid()

## ------------------------------------------------------------------------
if (!file.exists('dse.mu.climatrans.rda')) {
  index(pca.mu) <- as.Date(paste(year(pca.mu),'01-01',sep='-'))
  index(pca.fw) <- as.Date(paste(year(pca.fw),'01-01',sep='-'))
  dse.mu <- DSensemble(pca.mu,predictor=es,FUNX='tas2es',xfuns='tas2es',plot=TRUE)
  dse.fw <- DSensemble(pca.fw,predictor=slp,pattern="psl_Amon_ens_",plot=TRUE)
  save(file='dse.mu.climatrans.rda',dse.mu,dse.fw)
} else load('dse.mu.climatrans.rda')
mu.stations<-as.station(dse.mu)
fw.stations<-as.station(dse.fw)

## Plot the downscaled results for the three megacities:
for (i in 1:3) {
  plot(mu.stations[[i]],new=FALSE)
  plot(fw.stations[[i]],new=FALSE)
}

## ------------------------------------------------------------------------
## India:
Z2 <- precipcmip(lon=c(70,90),lat=c(7,30))
par(bty='n')
X <- zoo(x=apply(coredata(Z2$gcm),2,function(x) 100*(x/mean(x,na.rm=TRUE))),order.by=index(Z2$gcm))
plot(X,lwd=3,col=rgb(0,0.5,1,0.2),plot.type='single',ylab='Annual precip (%)')
lines(100*Z2$obs/mean(Z2$obs,na.rm=TRUE),lwd=3)
grid()

## ------------------------------------------------------------------------
  load('climatrans.tx.rda')
  plot(as.4seasons(anomaly(climatrans.tx)),new=FALSE)
  for(i in 1:3) lines(trend(subset(as.4seasons(anomaly(climatrans.tx)),is=i)))
  grid()
  print(summary(coredata(climatrans.tx)))

## ------------------------------------------------------------------------
## The predictor
  predictor <- retrieve('air.mon.mean.nc',param='air',lon=c(60,90),lat=c(10,35))
  predictor <- aggregate(subset(predictor,it='jja'),year,FUN='mean')
  index(predictor) <- year(predictor)
  eof.t2m <- EOF(predictor)

## ----Tx-predictand-------------------------------------------------------
## The predictand
  Tx <- subset(as.4seasons(climatrans.tx),it='jja')
  pca.tx <- PCA(Tx)
  class(pca.tx) <- c("pca","station", "annual", "zoo")
  index(pca.tx) <- year(pca.tx)

## ---- fig.height=8,fig.width=8-------------------------------------------
  ds.tx <- DS(pca.tx,eof.t2m)
  plot(ds.tx,new=FALSE)

## ---- fig.height=8,fig.width=8-------------------------------------------
## Check the residuals
  tx.2 <- as.station(as.residual(ds.tx))
  pca.tx.2 <- PCA(tx.2)
  plot(pca.tx.2,new=FALSE)

## ------------------------------------------------------------------------
  trh <- c(40,35,35)
  for (i in 1:3) {
    y <- subset(climatrans.tx,is=i)
    if (!file.exists(paste('dse.tx.climatrans.',loc(y),'.rda',sep=''))) {
      dse.tx <- DSensemble.t2m(y,biascorrect=TRUE,type='ncdf4',
                           predictor=predictor,nmin=60,verbose=FALSE)
      save(file=paste('dse.tx.climatrans.',loc(y),'.rda',sep=''),dse.tx)
    } else load(paste('dse.tx.climatrans.',loc(y),'.rda',sep=''))

    ## The function hotsummerdays estimates the number of hot days based on
    ## the assumption that the temperature is close to being normally distributed
    dse.tx[is.element(year(dse.tx),2100),] <- NA   # some of the results for 2100 are suspect 
    hw <- hotsummerdays(x=y,dse=dse.tx,threshold=trh[i],it=NULL,plot=TRUE,new=FALSE)
    #plot(hw)
  }

## ------------------------------------------------------------------------
  txq95 <- diagnose(subset(climatrans.tx,is=1))

## ------------------------------------------------------------------------
  tx.sd <- as.4seasons(anomaly(climatrans.tx),FUN='sd')
  plot(subset(tx.sd,it='jja'),new=FALSE)
  for(i in 1:3) lines(trend(subset(subset(tx.sd,it='jja'),is=i)))
  grid()
  print(trend(subset(tx.sd,it='jja'),result='pval'))
  print(colMeans(climatrans.tx))

## ------------------------------------------------------------------------
for (i in 1:3) {
  qqnorm(coredata(subset(subset(tx.sd,it='jja'),is=i)),main=loc(tx.sd)[i])
  qqline(coredata(subset(subset(tx.sd,it='jja'),is=i)),main=loc(tx.sd)[i])
}

## ------------------------------------------------------------------------
for(i in 1:3) {
  nhot1 <- annual(subset(climatrans.tx,is=i),FUN='count',threshold=trh[i])
  if (i==1) nhot <- nhot1 else nhot <- combine(nhot,nhot1)
}
plot(nhot,new=FALSE)
for(i in 1:3) lines(trend(subset(nhot,is=i)))
grid()
print(trend(nhot,result='pval')) 

## ------------------------------------------------------------------------
print(paste('Mean number=',mean(nhot[,1],na.rm=TRUE),'Variance=',var(nhot[,1],na.rm=TRUE),
            'should equal',mean(nhot[,1],na.rm=TRUE),'for Poisson process'))
hist(coredata(nhot[,1]),freq=FALSE,col='grey',
     main=loc(nhot)[1])
lines(dpois(x=seq(0,max(nhot[,1]),by=1),lambda=mean(nhot[,1],na.rm=TRUE)),lwd=2,col='red')
grid()

## ------------------------------------------------------------------------
N = round(10^seq(1,5,length=300))
qx <- rep(NA,length(N))
for (i in 1:length(N)) qx[i] <- max(rexp(N[i]))
plot(N,qx,log='x',main='Maximum value vs sample size',ylab='Maximum value')

## ------------------------------------------------------------------------
## Test weighted sum of daily rainfall.
## Simulate gridding 

N <- 100000
iw <- runif(4*N)
x <- rexp(4*N)
x[iw < 0.6] <- 0
W <- runif(4)
W <- W/sum(W)
dim(x) <- c(N,4)
y <- x %*% W

breaks <- seq(-1,20,by=0.1)
hx1 <- hist(x[,1],plot=FALSE,breaks=breaks)
hx2 <- hist(x[,2],plot=FALSE,breaks=breaks)
hx3 <- hist(x[,3],plot=FALSE,breaks=breaks)
hx4 <- hist(x[,4],plot=FALSE,breaks=breaks)
maxx <- max(hx1$counts,hx2$counts,hx3$counts,hx4$counts)
hy <- hist(y,breaks=breaks,xlim=c(0,5),ylim=c(0,maxx),col="grey",border="grey",
           main='Test gridded precip',xlab='amount',
           sub="'gridded' result = weighted sum of 4 stations")
lines(hx1$mids,hx1$counts,col='red',lwd=3)
lines(hx2$mids,hx2$counts,col='red',lwd=3)
lines(hx3$mids,hx3$counts,col='red',lwd=3)
lines(hx4$mids,hx4$counts,col='red',lwd=3)
legend(3,maxx,c('Individual simulation','Weighted sum'),
      col=c('red','grey'),lwd=c(3,10),bty='n')

## ----fig.width=12, fig.height=20-----------------------------------------
print('Compare common EOFs')
diag.slp <- biasdiag(slp,it=month.abb[season],verbose=FALSE)

## ----eval=FALSE----------------------------------------------------------
## 102197208 491010005000000000000000000000000000000000000000001090000000002820249 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 49102000000000000000000000000000300000000000000000000000000000000     1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 491110000000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 4911200000000000000000000000000000000000000000000000000000000         1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 491210000000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 49122000000000000000000000000000000000000000000000000000000000000     1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 500110000000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 50012000000000000000000000000000000000000000000000000000000000000     1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 500210000000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 50022000000000000000000000000000000000000000000000000                 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 500310000000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 50032000000000000000000000000000000000000000000000000000000000000     1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 500410000000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 5004200000000000000000000000000000000000000000000000000000000         1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 500510000000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 50052000000000000000000000000000000080003000500000000000000000000     1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 500610000000000000000000300050000000000000000004304750051001301470053 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 5006201220089005800790076000000200000000000000010001800380000         1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 500710005004100280000002000150513154901420030039108410729012200760051 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 50072004600080457074403560066006101070269046007320330001300710058     1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 500810028000300410147006100710081027701420152010700130015000000050005 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 50082000000250000000000030000000000000000000000000000005600030175     1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 500910470107207921283027400580208007900080000000003811201065800030000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 5009200030003002003350043000000380000000000000000000000000000         1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 501010028000000000000000000000000009900000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 50102000000000000000000000000000000000000000000000000000000000000     1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 501110000000000000000000000000000000000000000000000000000000000280000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 5011200000000000000000000000000000000000000000000000000000000         1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 501210000000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 50122000000000000000000000000000000000000000000000000000000000000     1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 51011    000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 51012000000000000000000000000000000000000000000000000000000000000     1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 51021    000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 51022000000000000000000000000000000000000000000000000                 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19
## 102197208 51031    000000000000000000000000000000000000000000000000000000000000 1207  BOMBAY CITY     BOMBAY SANTACRUZOBSY   19

## ----eval=FALSE----------------------------------------------------------
## Format Specifications
## Column                                           Data
## A                                                Catchment No. - 3 Digits
## B                                                Latitude - 2 Digits
## C                                                Longitudes - 2 Digits
## D                                                Station No. - 2 Digits
## E                                                Blank
## F                                                Year - 2 Digits
## G                                                Month - 2 Digits
## H                                                Card No. - 1 Digit
## I to W                                           Rainfall in 4 Digits (Without Decimal)
##                                                  If Col H has a value 1 - Rainfall from
## Dates 1 to 15; value 2 - Rainfall from Dates 17 to 31
## X                                                If Col H has value 1 - Rainfall for
## the date 16; value 2 - Rainfall for the month - (5 Digits)
## Y                                                State No. - 2 Digits
## Z                                                District No. - 2 Digits
## AA                                               STN No. - 2 Digits
## AB                                               District Name - 16 Characters
## AC                                               Station Name - 16 Characters
## AD                                               Class - 5 Characters
## AE                                               Secret Classification - 2 Digits
## AF                                               Century Year - 19 or 20

