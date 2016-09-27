library(esd)

data("geoborders",envir=environment())
ok <- (geoborders$y >= 65)

theta <- pi*geoborders$x[ok]/180; phi <- pi*geoborders$y[ok]/180
r <- 2.5*cos(phi)
x <- r*cos(theta)
y <- r*sin(theta)
x[abs(c(0,diff(theta))) > 0.1] <-  NA
y[abs(c(0,diff(theta))) > 0.1] <-  NA

par(bty="n",xaxt="n",yaxt="n",mar=rep(1,4))
plot(c(-1,1),c(-1,1),type='n',xlab="n",ylab="n")
lines(x,y,col="grey")

ss <- select.station(nmin=40,lat=c(65,90),param='t2m', it=c(1979,2012))
Theta <- pi*ss$longitude/180; Phi <- pi*ss$latitude/180

R <- 2.5*cos(Phi)
X <- R*cos(Theta)
Y <- R*sin(Theta)
points(X,Y,col='darkred',pch=19,cex=0.75)
