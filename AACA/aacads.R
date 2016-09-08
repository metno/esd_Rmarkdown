library(esd)

dayswithnosnow <- function(x) {
  x <- -x
  y <- annual(x,FUN='count',threshold=-1)
  attr(y,'unit') <- 'days < 1 cm'
  invisible(y)
}

load('sd.aaca.rda')
y <- dayswithnosnow(sd.aaca)
nok <- apply(coredata(y),2,nv)
y <- subset(y,is=nok > 40)
map(y,FUN='trend',cex=-2,colbar=list(breaks=seq(-22,22,by=2)))

