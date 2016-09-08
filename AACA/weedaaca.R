## Need to iterate a couple of times before all duplicates are removed.

library(esd)

#load('t2m.aaca.rda')
load('pr.aaca.rda')
t2m.aaca <- pr.aaca

lonlat <- paste(round(lon(t2m.aaca),3),round(lat(t2m.aaca),3))
n <- length(lonlat)
dlonlat <- rownames(table(lonlat[duplicated(lonlat)]))

## Index for stations to keep 
keep <- rep(TRUE,n)

for (i in 1:length(dlonlat)) {
  is <- is.element(lonlat,dlonlat[i])
  y <- subset(t2m.aaca,is=is)
  nv <- apply(coredata(y),2,nv)
  print(c(sum(is),loc(y),nv))
  remove <- (1:n)[is]
  if (sum(nv==max(nv))>1) {
     remove <- remove[is.element(nv,max(nv))][-1]
   } else
     remove <- remove[nv==max(nv)]
  keep[remove] <- FALSE
  plot(y)
}

print(sum(keep))
t2m.aaca <- subset(t2m.aaca,is=keep)
t2m.aaca <- subset(t2m.aaca,it=c(1900,2015))
esd::map(t2m.aaca,FUN='mean',cex=-2)

#save(file='t2m.aaca.rda',t2m.aaca)
save(file='pr.aaca.rda',pr.aaca)


if (TRUE) {
## Remove duplicated records:
dy <- duplicated(stid(Y))
d2 <- (1:dim(Y)[2])[dy]
for (i in 1:length(d2))
  print(apply(subset(t2m.aaca,is=c(d2[i],d2[i]-1)),2,nv))
Y <- subset(Y,is=!dy)
eval(parse(text=paste('Y -> ',param,'.aaca',sep='')))
eval(parse(text=paste('save(file="',param,'.aaca.rda",',param,'.aaca)',sep='')))
}

param <- 'pr'
FUN <- 'mean'
it <- 'son'
reanalysis <- 'air.mon.mean.nc'
 
if (param=='precip') FUNX <- 'C.C.eq' else FUNX <- 'mean'


## Remove unrealistically high amounts.
## Get the predictand -> Y
if (!exists('Y')) {
   load(paste(param,'.aaca.rda',sep=''))
   Y <- eval(parse(text=paste(param,'.aaca',sep='')))
 }
yc <- coredata(Y)
yc[yc > 500] <- NA
yc -> coredata(Y)
