## Rasmus Benestad
## Read IMILAST storms - select Bangladesh
## The IMILAST project excluded tropical cyclones.

library(esd)
# 03 000001 001 1989010100 1989    1    1   0   -170   61.5    961.24
imilast.storms <- '~/Disk1/IMILAST/imilast.cyclones/ERAinterim_1.5_NH_M03_19890101_20090331_ST.txt'
Z <- events2trajectory(read.imilast(imilast.storms))
##Z <- subset(Z,is=list(lon=c(80,120),lat=c(10,30)))
map(Z)

## Use the hurdat2 data.