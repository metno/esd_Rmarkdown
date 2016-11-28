## Evaluation of CMIP ensemble based on NCEP/NCAR 1, ERAINT, and MERRA reanalyses
## Surface temperature. Rank-test for grid-points values

testrank <- function(x,m,N,d1) {
  ## Find the rank number of m observations
  ## Need to unwrap the dimensions
  dim(x) <- c(N,d1)
  number <- apply(x,1,function(x) order(x)[1:m])
  return(c(number))
}

ECQ.ensemble <- function(obs=c('air.mon.mean.nc','ETAINT_t2m.mon.nc','MERRA'),
                         path='CMIP5.monthly',pattern='tas_',it=c(1980,2015),is=NULL,
                         anomaly=FALSE) {
  require(esd)
  
  ## Get all the GCM data
  gcms <- list.files(path=path,pattern=pattern,full.names=TRUE)
  n <- length(gcms)
  N <- n + length(obs)
  m <- length(obs)
  ## Use the grid of the first reanalysis and interploate the other fields onto this grid
  y1 <- subset(retrieve(obs[1]),it=it,is=is)
  if (anomaly) y1 <- anomaly(y1)
  d <- dim(y1)
  
  ## Set up a data matrix containing the data - start with the first reanalysis
  X <- rep(NA,N*d[1]*d[2]); dim(X) <- c(N,d)
  cnames <- rep("",N)
  X[1,,] <- coredata(y1)
  cnames[1] <- paste(attr(y1,'model_id'),attr(y1,'run'))
  
  if (m > 1) {
    for (i in 2:m) {
      y <- subset(retrieve(obs[i]),it=it,is=is)
      y <- regrid(y,is=y1)
      if (anomaly) y <- anomaly(y)
      X[i,,] <- coredata(y)
      cnames[i] <- paste(attr(y,'model_id'),attr(y,'run'))
    }
  }
  
  ## Do the GCMs:
  for (i in 1:n) {
    y <- subset(retrieve(gcms[i]),it=it,is=is)
    y <- regrid(y,is=y1)
    if (anomaly) y <- anomaly(y)
    X[i+m,,] <- coredata(y)
    cnames[i+m] <- paste(attr(y,'model_id'),attr(y,'run'))
  }
  
  ## convert the matrix to 2D
  dim(X) <- c(N*d[1],d[2])
  
  ## For each gridpoint, accumulate the ranks for all the time steps.
  ## The dimension of the results should be [m*d[1],d[2]], where te elements of the 
  ## first dimension is expected to follow a uniform distribution if the GCMs and reanalysis
  ## belong to the same population.
  testmap <- apply(X,1,testrank,m,N,d[1])
}