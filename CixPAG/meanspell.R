## Script that reads European temperature data and explores the connection between 
## the seasonal mean temperature and the mean length of the warm/cold spells 

it <- 'djf'
if (it =='djf') {
  cold <- TRUE
  threshold <- 0
  is <- 2
} else {
  cold <- FALSE
  threshold <- 20
  is <- 1
}

AR <- function(n,mean=1,sd=1,a1=0.8) {
  rn <- rnorm(n,mean=mean,sd=sd)
  for (i in 2:n) rn[i] <- (a1*rn[i-1] + (1-a1)*rn[i])
  invisible(rn)
}

library(esd)
ss <- select.station(src='ecad',param='tg',nmin=75)
d <- dim(ss)
x <- rep(NA,d[1]); y <- x; std <- y
q.spell <- rep(NA,d[1]*10); dim(q.spell) <- c(d[1],10); q.geom <- q.spell

for (i in 1:d[1]) {
  z <- station(ss[i,])
  print(loc(z))
  ## Make sure that there are values above and below the given threshold - otherwise
  ## spell will not work.
  if ( (sum(z > threshold,na.rm=TRUE)>1000) & (sum(z < threshold,na.rm=TRUE)> 1000) ) {
    s <- spell(z,threshold=threshold)
    y[i] <- mean(subset(subset(s,is=is),it=it),na.rm=TRUE)
    
    ## Compare the spell-distribution with a geometric distribution
    q.spell[i,] <- quantile(subset(subset(s,is=is),it=it),probs=seq(0.05,0.95,by=0.1),na.rm=TRUE)
    q.geom[i,] <- qgeom(p=seq(0.05,0.95,by=0.1),prob=1/y[i])
    std[i] <- sd(subset(z,it=it),na.rm=TRUE)
    x[i] <- mean(subset(z,it=it),na.rm=TRUE)
    plot(x,y,main=paste('Mean',it,'temperature & mean length of intervals',
                        c('below','above')[c(cold,!cold)],threshold,'C'),
        sub='source: ECA&D',pch=19,col=rgb(0.5,0,0,0.3),
        xlab=expression(bar(T)),ylab=expression(bar(tau)))
  }
  grid()
}

## Plot results

plot(x,y,main=paste('Mean',it,'temperature & mean length of intervals',
                    c('below','above')[c(cold,!cold)],threshold,'C'),
     sub='source: ECA&D',pch=19,col=rgb(0.5,0,0,0.3),
     xlab=expression(bar(T)),ylab=expression(bar(tau)))
grid()

## Monte-Carlo simulations to compare spell length with 
mstd <- 2*max(std,na.rm=TRUE)
if (!cold) mx <- mean(x,na.rm=TRUE) else mx <- 0
nmc <- 300
ymc <- rep(NA,nmc); xmc <- ymc
for (m in seq(mx-mstd,mx+mstd,length=nmc)) {
  coredata(z) <- AR(length(z),mean=m,sd=mstd)
  s <- spell(z,threshold=threshold)
  ymc[i] <- mean(subset(subset(s,is=is),it=it),na.rm=TRUE)
  xmc[i] <- m
  points(xmc,ymc,pch=19,col='grey75')
}

points(x,y,pch=19,col=rgb(0.5,0,0,0.3))

ix <- order(x); x <- x[ix]; y <- y[ix]
ok <- is.finite(x) & is.finite(y)
x <- x[ok]; y <- y[ok]; std <- std[ok]
calfit <- data.frame(x=x,y=y)
attr(calfit,'max(x)') <- max(x,na.rm=TRUE)
#fit <- lm(y ~ I(x) + I(x^2),data=calfit)
fit <- glm(y ~ x,data=calfit,family='poisson')
lines(calfit$x,exp(predict(fit)),col='red')

attr(x,'description') <- paste(it,'mean temperature (degC)')
if (cold) attr(y,'description') <- 'mean cold spell length (days)' else
          attr(y,'description') <- 'mean warm spell length (days)'
attr(y,'threshold') <- threshold
attr(x,'label') <- expression(bar(T))
attr(x,'Monte-Carlo') <- xmc
attr(y,'label') <- expression(bar(tau[T < T0]))
attr(y,'Monte-Carlo)') <- ymc
meanspell <- data.frame(meanT=x,meanL=y,std=std)
attr(meanspell,'fit') <- fit
attr(meanspell,'geometric.fit') <- data.frame(q.spell=q.spell,q.geom = q.geom)
save(meanspell,file=paste('meanspell',it,c('below','above')[c(cold,!cold)],threshold,'.rda',sep=''))

## Test if the spell length statistics is close to geometric
plot(c(q.spell),c(q.geom),main='Spell length statistics & the geometric distribution',
     xlim=range(q.spell,q.geom,na.rm=TRUE),
     ylim=range(q.spell,q.geom,na.rm=TRUE),pch=19,col=rgb(0,0,0,0.2),
     xlab=expression(q[p]),ylab='qgeom(p,1/mean)')
grid()
lines(range(q.spell,q.geom,na.rm=TRUE),range(q.spell,q.geom,na.rm=TRUE),col='red')

## Pick Oslo
t2m <- subset(as.4seasons(station(ss[63,])),it=it)
sl <- subset(as.4seasons(spell(station(ss[63,]),threshold=threshold)),is=is,it=it)
t2msl <- zoo(exp(predict(fit,newdata=data.frame(x=coredata(t2m)))),order.by=index(t2m))

xy <- merge(t2msl,sl)
ok <- is.finite(xy$t2msl) & is.finite(xy$sl)
r <- round(cor(coredata(xy$t2msl)[ok],coredata(xy$sl)[ok]),2)

plot(coredata(xy$t2msl),coredata(xy$sl),main=paste('Correlation=',r))

plot(xy,plot.type='single',lwd=3,
     xlab='',ylab='days',sub=paste(loc(t2m),'correlation=',r),
     main=paste('Mean',toupper(it),'length of intervals',
                c('below','above')[c(cold,!cold)],threshold,'C'),
     col=c(rgb(0.7,0,0.2,0.5),rgb(0,0.2,0.7,0.5)))
grid()
legend(as.Date('1997-01-01'),5.75,c('Predicted','Observed'),
       col=c(rgb(0.7,0,0.2,0.5),rgb(0,0.2,0.7,0.5)),lty=1,lwd=3)

