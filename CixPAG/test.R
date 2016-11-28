## Test script:

projection <- list()
load('dse.aaca.t2m.rcp45.djf.eof.rda')
projection$djf <- Z
load('dse.aaca.t2m.rcp45.mam.eof.rda')
projection$mam <- Z
load('dse.aaca.t2m.rcp45.jja.eof.rda')
projection$jja <- Z
load('dse.aaca.t2m.rcp45.son.eof.rda')
projection$son <- Z

print('map')
map(Z[[1]])