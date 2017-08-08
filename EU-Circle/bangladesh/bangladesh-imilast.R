## Rasmus Benestad
## Read IMILAST storms - select Bangladesh
## The IMILAST project excluded tropical cyclones.

library(esd)

find.file <- function(filename) {
  command <- paste("find -L $HOME -name",filename,sep=" ")
  fullpath <- system(command,intern=TRUE)
  if(length(fullpath)==0) return(FALSE)
  return(fullpath)
}

path <- sub("/bangladesh-imilast.R","",find.file("bangladesh-imilast.R")[1])

## HURDAT2 CANNOT BE USED FOR BANGLADESH.
## IT CONTAINS ONLY ATLANTIC (AND SOME PACIFIC?) STORMS. 
#hurdat2.storms <- find.file('hurdat2-1851-2016-apr2017.txt')[1]
#Y <- read.hurdat2(hurdat2.storms)

## IMILAST data for Bangladesh and the surrounding region
is <- list(lon=c(60,120),lat=c(0,40))
fname.imilast <- file.path(path,'imilast_bangladesh_ERAinterim_1.5.rda')
if(is.character(file.exists(fname.imilast))) {
  load(fname.imilast)
} else {
  fname <- find.file('ERAinterim_1.5_NH_*.txt')
  fname <- fname[!duplicated(gsub(".*/","",fname))]
  M <- gsub("_.*","",gsub(".*M","",fname))
  imilast <- list()
  for (i in seq(length(fname))) {
    print(paste("imilast M",M[i],sep=""))
    x <- subset(read.imilast(fname[i]),is=is)
    imilast[[paste("M",M[i],sep="")]] <- x
  }
  save(file=fname.imilast,imilast)
}

## Map cyclone trajectories
save.fig <- TRUE
is <- list(lon=c(75,105),lat=c(0,30))
## Cyclone 1BOB01, April 29, 1991
t1 <- seq(strptime("19910428 12",format="%Y%m%d %H"), 
          strptime("19910501 00",format="%Y%m%d %H"), by = "12 hours")
## May 1997 Bangladesh cyclone
t2 <- seq(strptime("19970515 12",format="%Y%m%d %H"), 
          strptime("19970518 00",format="%Y%m%d %H"), by = "12 hours")
## Cyclone Sidr, November 15, 2007
t3 <- seq(strptime("20071114 12",format="%Y%m%d %H"), 
          strptime("20071117 00",format="%Y%m%d %H"), by = "12 hours")
tvec <- list(t1,t2,t3)
M <- names(imilast)
for(m in M) {
  print(m)
  x <- imilast[[m]]
  param <- grep("pressure|slp|intensity|min",names(x),value=TRUE)[1]
  cb <- list(pal="bu",rev=TRUE,show=TRUE)
  if(grepl("slp|pressure",param)) {
    unit <- "hPa"
    if(grepl("depth",param)) {
      cb$breaks <- seq(0,10,2)
      cb$rev <- FALSE
    } else {
      cb$breaks <- seq(995,1025,5)
      cb$rev <- TRUE
    }
  } else if(grepl("z",param)) {
    unit <- "m^2/s^2"
    cb$breaks <- pretty(range(x[param]),n=10)
    cb$rev <- TRUE
  } else {
    unit <- NULL
    cb$breaks <- pretty(range(x[param]),n=10)
    cb$rev <- FALSE
  }
  if(m=="M14") cb$breaks <- pretty(c(0,60),n=10)
  for(times in tvec) {
    if(save.fig) {
      fname.fig <- file.path(path,paste("Figures/imilast.trajectories.bangladesh",
                             m,strftime(times[2],format="%Y%m%d"),"pdf",sep="."))
      pdf(file=fname.fig,width=10,height=6)
    } else {
      x11(width=10,height=6)
      #dev.new(width=10,height=6)
    }
    for(i in seq(1,min(6,length(t1)))) {
      if(i==1) {
        new <- FALSE
        main <- m
      } else {
        new <- TRUE
        main <- ""
      }
      if(i==4) {
        cb$show <- TRUE
        label.param <- param
        if(!is.null(unit)) label.param <- paste(label.param,"  (",unit,")",sep="")
      } else {
        cb$show <- FALSE
        label.param <- ""
      }
      yi <- ceiling(i/3)
      xi <- i-(yi-1)*3
      fi <- c((xi-1)/3,xi/3,(1-yi/2)*0.95+0.05,(1-(yi-1)/2)*0.95+0.05)
      par(fig=fi,new=new)
      map(x,it=times[i],xlim=is$lon,ylim=is$lat,projection="lonlat",
          main=main,param=param,label.param=label.param,
          lwd=3.5,alpha=0.5,border=TRUE,
          type=c("points","trajectory","colors","start"),
          mgp=c(1,0.5,0),mar=c(1.5,2,2,2),showaxis=TRUE,
          colbar=cb,verbose=FALSE,fig=fi,new=FALSE,add=new)
    }
    if(save.fig) dev.off()
  }
}
