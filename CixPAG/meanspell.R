## Script that reads European temperature data and explores the connection between 
## the seasonal mean temperature and the mean length of the warm/cold spells 

it <- 'djf'
cold <- TRUE
if (cold) is <- 2 else is <- 1
threshold <- 0

library(esd)
ss <- select.station(src='ecad',param='tg',nmin=75)
d <- dim(ss)
x <- rep(NA,d[1]); y <- x
for (is in 1:d[1]) {
  Tx <- station(ss[is,])
  print(loc(Tx))
  z <- subset(Tx,is=is)
  ## Make sure that there are values above and below the given threshold - otherwise
  ## spell will not work.
  if ( (sum(z > threshold,na.rm=TRUE)>100) & (sum(z < threshold,na.rm=TRUE)> 100) ) {
    s <- spell(z,threshold=threshold)
    y[is] <- mean(subset(subset(s,is=2),it=it),na.rm=TRUE)
    x[is] <- mean(subset(z,it=it),na.rm=TRUE)
    plot(x,y,main='Mean temperature & mean length of freezing spells',
        sub='source: ECA&D',pch=19,col=rgb(0.5,0,0,0.3),
        xlab=expression(bar(T)),ylab=expression(bar(tau[T < T0])))
    calfit <- data.frame(x=max(x,na.rm=TRUE) - x[is.finite(x) & is.finite(y)],
                         y=y[is.finite(x) & is.finite(y)])
    attr(calfit,'max(x)') <- max(x,na.rm=TRUE)
    #fit <- lm(y ~ I(x) + I(x^2),data=calfit)
    #fit <- glm(y ~ x,data=calfit,family='poisson')
    #if (is > 10) lines(calfit$x[order(cal.fit$x)]-max(x,na.rm=TRUE),
    #                   predict(fit)[order(calfit$x)],col='red')
  }
  grid()
}

attr(x,'description') <- paste(it,'mean temperature (degC)')
if (cold) attr(y,'description') <- 'mean cold spell length (days)' else
          attr(y,'description') <- 'mean warm spell length (days)'
attr(y,'threshold') <- treshold
attr(x,'label') <- expression(bar(T))
attr(y,'label') <- xpression(bar(tau[T < T0]))
meanspell <- data.frame(meanT=x,meanL=y)
attr(meanspell,'fit') <- fit
save(meanspell,file=paste('meanspell',it,c('below','above')[c(cold,!cold)],threshold,'.rda',sep=''))


