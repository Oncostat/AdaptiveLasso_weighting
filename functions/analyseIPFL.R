### Working directory
setwd(" ")

#####################################################################
# R packages installation                                           #
#####################################################################
library(dplyr)
library(survival) 
library(Matrix)
library(corpcor)
library(parallel)
library(ipflasso)

set.seed( )

pflist <-
  list(
    c(1, rep(2, 19)),
    c(2, 1, rep(2, 18)),
    c(rep(2, 10), 1, rep(2, 9)),
    c(1, 1, rep(2, 18)),
    c(1, rep(2, 9), 1, rep(2, 9)),
    c(2, 1, rep(2, 8), 1, rep(2, 9)),
    c(1, 1, rep(2, 8), 1, rep(2, 9))
  )

nfolds = 5

# Number of PC-cores used
pCores <- #80

# Definition of the cluster via 'makeCluster'
cl <- makeCluster(pCores, outfile = "outfilename.txt")
r <- 500

v <-
  paste0(
    './datas/scen', #path
    rep(1:8, each = r),
    '/simdata-',
    sprintf("%03d", c(1:r)),
    '.Rdata'
  )

set.seed(25845)
analyse <- function(character) {
  load(character)
  source("./functions/ipflasso.R")
  respath <- gsub('(datas|data-)', 'resultsIPFL', character)
  
  if (isTRUE(file.exists(respath))) {
    print('Found')
    flush.console()
  } else{
    res <- ipf.lasso(data[1:500, ], pflist, nfolds)
    flush.console()
    save(res, file = respath)
    print('File Saved')
    flush.console()
  }
}

# Execution of the program
start <- Sys.time()
z <- clusterApplyLB(cl = cl, v, analyse)
end <- Sys.time()
t.parallel <- difftime(end, start, units = "sec")
stopCluster(cl)
