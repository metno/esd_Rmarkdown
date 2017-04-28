## Organise decadal forecasts for SPECS
## Rasmus.Benestad@met.no
## Did not manage to install the R-package, so I source the R-code instead
rscripts <- list.files(path='~/Downloads/loadR/',pattern='.R',full.names = TRUE)
rscripts <- rscripts[-grep('.Rd',rscripts,fixed=TRUE)]
rscripts <- rscripts[-grep('.dic',rscripts,fixed=TRUE)]
rscripts <- rscripts[-grep('DESCRIPTION',rscripts,fixed=TRUE)]
for (i in 1:length(rscripts)) source(rscripts[i])

tasfiles <- list.files(path="~/Downloads/SPECS-decadal-fcst_hadcm3_i2p1/",pattern="tas",full.names = TRUE)

tasDECA <- loadDecadalForecast(dataset = tasfiles[1], # or  "mydata.nc"
                               latLim = c(0,80),
                               lonLim = c(-90,20),
                               var = "tas", 
                               dictionary = F,
                               years = 1971:2010,
                               season = 1:12)