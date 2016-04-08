library('esd')
data("IPCC.AR5.Table.9.A.1")
data("global.t2m.gcm")
gcmnm <- attr(global.t2m.cmip5,"model")
gcmnm <- toupper(substr(gcmnm,1,nchar(gcmnm)-2))
gcmnm[is.element(substr(gcmnm,1,13),"CSIRO-MK3-6-0")] <- "CSIRO-MK3.6.0"
gcmnm[is.element(gcmnm,"FGOALS_G2")] <- "FGOALS-G2"
gcmnm[is.element(gcmnm,"FIO-ESM")] <- "FIO-ESM V1.0"
gcmnm[is.element(gcmnm,"BCC-CSM1-1")] <- "BCC-CSM1.1"
gcmnm[is.element(gcmnm,"BCC-CSM1-1-M")] <- "BCC-CSM1.1(M)"
gcmnm[is.element(gcmnm,"CESM1-BGC")] <- "CESM1(BGC)"
gcmnm[is.element(gcmnm,"CESM1-CAM5")] <- "CESM1(CAM5)"
gcmnm[is.element(gcmnm,"CNRM-CM5")] <- "CNRM-CM51"
gcmnm[is.element(gcmnm,"INMCM4")] <- "INM-CM4"
gcmnm[is.element(substr(gcmnm,1,7),"EC-EART")] <- "EC-EARTH"
n <- length(gcmnm); m <- 12
X <- matrix(rep(NA,n*m),n,m)
for (i in 1:n) {
  ii <- is.element(toupper(IPCC.AR5.Table.9.A.1$Model.Name),gcmnm[i])
  if (sum(ii)>0)
    X[i,] <- c(as.character(IPCC.AR5.Table.9.A.1$Atmosphere.Component.Name)[ii][1],
               as.character(IPCC.AR5.Table.9.A.1$Horizontal.Grid)[ii][1],
               IPCC.AR5.Table.9.A.1$Number.of.Vert.Levels[ii][1],
               IPCC.AR5.Table.9.A.1$Grid.Top.hPa[ii][1],
               as.character(IPCC.AR5.Table.9.A.1$Component.Name.or.type)[ii][1],
               as.character(IPCC.AR5.Table.9.A.1$Atmos.Chemistry.Component.Name)[ii][1],
               as.character(IPCC.AR5.Table.9.A.1$Land.Surface.Component.Name)[ii][1],
               as.character(IPCC.AR5.Table.9.A.1$Ocean.Component.Name)[ii][1],
               IPCC.AR5.Table.9.A.1$Number.of.Vertical.Levels[ii][1],
               as.character(IPCC.AR5.Table.9.A.1$z.Co.ord)[ii][1],
               as.character( IPCC.AR5.Table.9.A.1$Ocean.Biogeochemistry.Component.Name)[ii][1],
               as.character(IPCC.AR5.Table.9.A.1$Sea.Ice.Component.Name)[ii][1]) else
                 print(gcmnm[i])
}
colnames(X) <- c('Atm.model','atm.grid',
                 'atm.lev','atm.top',
                 'atm.type','Atm.chem',
                 'Land','Oce.model',
                 'oce.lev','oce.z',
                 'oce.bio','Seaice')
dT <- colMeans(window(zoo(global.t2m.cmip5),start=2070,end=2099))-
      colMeans(window(zoo(global.t2m.cmip5),start=1970,end=1999))
X <- as.data.frame(cbind(round(dT,2),X))

fit.all <- lm(dT ~ Atm.model + stm.grid + atm.lev + atm.top +
              atm.type + Atm.chem + Land + Oce.model + oce.lev + oce.z +
              oce.bio + Seaice,data=X)
print(summary(fit.all))

fit.atm <- lm(dT ~ Atm.model + stm.grid + atm.lev + atm.top,data=X)
print(summary(fit.atm))

fit.oce <- lm(dT ~ Oce.model + oce.lev + oce.z,data=X)
print(summary(fit.oce))
