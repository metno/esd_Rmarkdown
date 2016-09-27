## author: K Parding
## storm track analysis for paper 59

library(esd)

path.cmip <- "/vol/klimadata/work/abdelkaderm/CMIP5.monthly"
path.pre <- "/vol/fou/klima/kajsamp/data"

slp <- retrieve(file.path(path.pre,'slp.mon.mean.nc'),
                lon=c(-180,180),lat=c(0,90),type="ncdf4")
attr(slp,"unit") <- "hPa" ## millibar = hPa

#data("storms.barents")
is.barents <- list(lon=c(10,70),lat=c(65,80))
eofs <- 1:5

load("/vol/fou/klima/kajsamp/MIST2/cyclones.NH.rda")
barents <- subset(Y,it=Y$trajectory %in% subset(Y,is=is.barents)$trajectory)
#barents <- trackfilter(barents,param="trackcount",pmin=5)
#barents <- trackfilter(barents,param="tracklength",pmin=500)
barents <- subset(barents,it=c(1979,2011))

#save(file="storms.barents.rda",barents)
deep <- trackfilter(barents,param="pcent",pmax=980)

N <- as.station(subset(deep,is=is.barents))
attr(N,"param") <- "cyclones"
attr(N,"longname") <- "number of deep cyclones in the Barents region"
attr(N,"unit") <- "events/month"
N.year <- aggregate(N,by=year,FUN="sum")
attr(N.year,"unit") <- "events/year"

plot(aggregate(N,by=month,FUN="mean"),xlab="month")
dev.copy2pdf(file="barents.deepcyclones.seasoncycle.pdf")
dev.off()

## Downscale the number of deep cyclones in the Barents region

test.ds <- FALSE
if (test.ds) {
  is.predictor <- list(lon=c(-180,180),lat=c(60,90))
  predictor <- EOF(aggregate(subset(slp,is=is.predictor),by=year,FUN="mean"))
  ds.N <- DS(N.year,predictor,eofs=eofs)
  plot(ds.N)
  predictor2 <- EOF(aggregate(subset(slp,
   is=list(lon=c(-90,90),lat=c(30,80))),
   by=year,FUN="mean"))
  ds2.N <- DS(N.year,predictor2,eofs=eofs)
  plot(ds2.N)
  predictor3 <- EOF(aggregate(subset(slp,
   is=list(lon=c(-180,180),lat=c(30,90))),
   by=year,FUN="mean"))
  ds3.N <- DS(N.year,predictor3,eofs=eofs)
  plot(ds3.N)
}

force <- FALSE
eofs <- 1:5
biascorrect <- TRUE
select <- NULL
for (rcp in c("rcp26","rcp45","rcp85")) {
  fname <- paste("DSensemble",rcp,"cyclones.barents.annual.rda",sep=".")
  if (!file.exists(fname) | force) {
  dse.N <- DSensemble(N.year,path=path.cmip,rcp=rcp,
                  pattern="psl_Amon_ens_",
                  predictor=slp,
                  biascorrect=biascorrect,eofs=eofs,select=select,
                  plot=TRUE,verbose=TRUE,
                  lon=c(-180,180),
                  lat=c(60,90),rel.cord=FALSE)
    save(file=fname,dse.N)
  } else {
    load(fname)
  }

  fname <- paste("DSensemble",rcp,"cyclones.barents.annual.M09.rda",sep=".")
  if (!file.exists(fname) | force) {
    path.imilast <- "/vol/fou/klima/IMILAST"
    file.imilast <- "ERAinterim_.75_NH_M09_19790101_20111231_ST.txt"
    m <- read.imilast(file.imilast,path=path.imilast,verbose=TRUE)
    m <- Trackstats(m,verbose=TRUE)
    Z2 <- subset(m,it=m$trackcount>=5)
    Z2 <- subset(Z2,it=Z2$tracklength>=500)
    barents2 <- subset(Z2,is=is.barents,it=c(1979,2011))
    deep2 <- subset(barents2,it=barents2$slpmin<980)
    N2 <- as.station(deep2) # number of cyclones
    attr(N2,"variable") <- "cyclones"
    attr(N2,"longname") <- "number of deep cyclones in the Barents region"
    attr(N2,"unit") <- "events/month"
    N2.year <- aggregate(N2,by=year,FUN="sum")
    attr(N2.year,"unit") <- "events/year"
    dse.m09 <- DSensemble(N2.year,path=path.cmip,rcp=rcp,
                  pattern="psl_Amon_ens_",
                  predictor=slp,
                  biascorrect=biascorrect,eofs=eofs,select=select,
                  plot=TRUE,verbose=TRUE,
                  lon=c(-90,90),
                  lat=c(30,80))
    save(file=fname,dse.m09)
  } else {
    load(fname)
  }

  fname <- paste("DSensemble",rcp,"cyclones.barents.annual.M16.rda",sep=".")
  if (!file.exists(fname) | force) {
    file.imilast <- "ERAINTERIM_0.75_NH_M16_19790101_20091231_ST.txt"
    m <- read.imilast(file.imilast,path=path.imilast,verbose=TRUE)
    m <- Trackstats(m,verbose=TRUE)
    Z3 <- subset(m,it=m$trackcount>=5)
    Z3 <- subset(Z3,it=Z3$tracklength>=500)
    barents3 <- subset(Z3,is=is.barents,it=c(1979,2011))
    deep3 <- subset(barents3,it=barents3$depth>3000)
    N3 <- as.station(deep3) # number of cyclones
    attr(N3,"variable") <- "cyclones"
    attr(N3,"longname") <- "number of deep cyclones in the Barents region"
    attr(N3,"unit") <- "events/month"
    N3.year <- aggregate(N3,by=year,FUN="sum")
    attr(N3.year,"unit") <- "events/year"
    dse.m16 <- DSensemble(N3.year,path=path.cmip,rcp=rcp,
                  pattern="psl_Amon_ens_",
                  predictor=slp,
                  biascorrect=biascorrect,eofs=eofs,select=select,
                  plot=TRUE,verbose=TRUE,
                  lon=c(-90,90),
                  lat=c(30,80))
    save(file=fname,dse.m16)
  } else {
    load(fname)
  }

  plot(dse.N,map.show=TRUE,map.type="rectangle")
  title(rcp,cex.main=1.5,line=-1)
  dev.copy2pdf(file=paste("dse.barents.deepcyclones",rcp,"pdf",sep="."))
  dev.off()
  plot(dse.m09,map.show=TRUE,map.type="rectangle")
  title(paste(rcp,"M09"),cex.main=1.5,line=-1)
  dev.copy2pdf(file=paste("dse.barents.deepcyclones.m09",rcp,"pdf",sep="."))
  dev.off()
  plot(dse.m16,map.show=TRUE,map.type="rectangle")
  title(paste(rcp,"M16"),cex.main=1.5,line=-1)
  dev.copy2pdf(file=paste("dse.barents.deepcyclones.m16",rcp,"pdf",sep="."))
  dev.off()
  
}

stormdensity <- FALSE
if (stormdensity) {
  if (!file.exists("cdens.atlantic.rda")) {
    atlantic <- subset(Z,is=list(lon=c(-90,120),lat=c(40,90)),it=c(1979,2011))
    ok <- atlantic$pcent<980
    deeptracks <- unique(atlantic$trajectory[ok])
    it <- sapply(atlantic$trajectory,function(x) x %in% deeptracks)
    W <- subset(atlantic,it=it)
    cdens <- as.field(W,dx=1,dy=1,verbose=TRUE)
    save(file="cdens.atlantic.rda",cdens)
  } else {
    load("cdens.atlantic.rda")
  }
  X <- subset(cdens,is=is.barents)
  attr(X,"longname") <- "track density of deep cyclones"
  attr(X,"variable") <- "track~density"
  X.cycle <- aggregate(X,by=month,FUN="mean")
  X.4s <- as.4seasons(X)
  X.djf <- subset(X.4s,it="djf")
  X.mam <- subset(X.4s,it="mam")
  X.jja <- subset(X.4s,it="jja")
  X.son <- subset(X.4s,it="son")
  map(X,FUN="mean",colbar=list(pal="precip",breaks=seq(0,1.5,0.1)))
  dev.copy2pdf(file="barents.deepcyclones.trackdensity.pdf")
  dev.off()
  cb <- list(pal="precip",breaks=seq(0,2.6,0.2))
  par(fig=c(0,0.5,0.53,1),mar=c(2,2,2,1),xaxt="n",yaxt="n",bty="n")
  map(X.djf,colbar=cb,main="djf")
  par(fig=c(0.5,1,0.53,1),new=TRUE,mar=c(2,2,2,1),xaxt="n",yaxt="n",bty="n")
  map(X.mam,colbar=list(show=FALSE,pal=cb$pal,breaks=cb$breaks),main="mam")
  par(fig=c(0,0.5,0.03,0.5),new=TRUE,mar=c(2,2,2,1),xaxt="n",yaxt="n",bty="n")
  map(X.jja,colbar=list(show=FALSE,pal=cb$pal,breaks=cb$breaks),main="jja")
  par(fig=c(0.5,1,0.03,0.5),new=TRUE,mar=c(2,2,2,1),xaxt="n",yaxt="n",bty="n")
  map(X.son,colbar=list(show=FALSE,pal=cb$pal,breaks=cb$breaks),main="son")
  dev.copy2pdf(file="barents.deepcyclones.trackdensity.4s.pdf")
  dev.off()
  plot(EOF(X.cycle))
  dev.copy2pdf(file="barents.deepcyclones.trackdensity.eof.cycle.pdf")
  dev.off()
  slp.barents <- subset(slp,is=is.barents)#list(lon=c(-60,120),lat=c(50,80)))
  slp.cycle <- aggregate(slp.barents,by=month,FUN="mean")
  plot(EOF(slp.cycle),colbar=list(pal="t2m",breaks=seq(-0.2,0.2,0.05)))
  dev.copy2pdf(file="eof.cycle.slp.pdf")
  dev.off()
}


########################################################
## Seasons
########################################################
## N.4s <- as.4seasons(N)
## slp.4s <- as.4seasons(subset(slp,is=list(lon=c(-180,180),lat=c(30,90)))) 
## ds.djf <- DS(subset(N.4s,it="djf"),EOF(subset(slp.4s,it="djf")),eofs=eofs)
## ds.mam <- DS(subset(N.4s,it="mam"),EOF(subset(slp.4s,it="mam")),eofs=eofs)
## ds.jja <- DS(subset(N.4s,it="jja"),EOF(subset(slp.4s,it="jja")),eofs=eofs)
## ds.son <- DS(subset(N.4s,it="son"),EOF(subset(slp.4s,it="son")),eofs=eofs)     
## dse.djf <- DSensemble.season(N,season="djf",path=path.cmip,rcp=rcp,
##                   pattern="psl_Amon_ens_",
##                   predictor=slp,FUN="sum",nmin=3,
##                   biascorrect=TRUE,eofs=eofs,select=select,
##                   plot=TRUE,verbose=TRUE,
##                   lon=c(-180,180),
##                   lat=c(30,90))

## Ns <- as.seasons(N,start="10-01",end="03-31",FUN="sum")
## Ns <- Ns[attr(Ns,"n.valid")==6]
## slps <- as.seasons(slp,start="10-01",end="03-31",FUN="sum")
## slps <- subset(slps,it=attr(slps,"n.valid")==6)
##   EOF(subset(slp.4s,it="djf"))
## ds.djf <- DS(subset(N.4s,it="djf"),EOF(subset(slp.4s,it="djf")),eofs=eofs)
## dse.djf <- DSensemble.season(N,season="djf",path=path.cmip,rcp=rcp,
##                   pattern="psl_Amon_ens_",
##                   predictor=slp,FUN="sum",nmin=3,
##                   biascorrect=TRUE,eofs=eofs,select=select,
##                   plot=TRUE,verbose=TRUE,
##                   lon=c(-180,180),
##                   lat=c(30,90))

##plot.events <- function(x,it=NULL,is=NULL,param="count",FUN="mean",
##                        verbose=FALSE,...) {
##  if(verbose) print("plot.events")
##  x <- subset(x,it=it,is=is)
##  N <- as.station(x,param=param,FUN=FUN,verbose=verbose)
##  plot(N,verbose=verbose,...)
##}
