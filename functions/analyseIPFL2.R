### Working directory ###
#setwd("...")

### R packages ###
library(dplyr)
library(survival) 
library(Matrix)
library(corpcor)
library(parallel)
library(ipflasso)
 

### Set a seed ###
#set.seed(...)


pflist <-
  list(
   c(2, rep(1, 19)),
   c(1, 2, rep(1, 18)),
   c(rep(1, 10), 2, rep(1, 9)),

   c(2, 2, rep(1, 18)),
   c(2, rep(1, 9), 2, rep(1, 9)),
   c(1, 2, rep(1, 8), 2, rep(1, 9)),
   c(2, 2, rep(1, 8), 2, rep(1, 9))
 )

nfolds <- 5

### Number of PC-cores used ###
pCores <- 80 

### Definition of the cluster via 'makeCluster' ###
cl <- makeCluster(pCores, outfile ="outfilename.txt")
r <- 500

v <- paste0('./datas/scen',rep(1:8,each=r),'/simdata-',sprintf("%03d",c(1:r)),'.Rdata')

set.seed(258452)
analyse <- function(character){
  load(character)
  source("./functions/ipflasso2.R")
  respath <- gsub('(datas|data-)', 'resultsIPFL2', character)
  
  if (isTRUE(file.exists(respath))){
    print('Found')
    flush.console()
  }else{
    res <- ipf.lasso2(data[1:500,], pflist, nfolds)          
    flush.console()
    save(res, file = respath)
    print('File Saved')
    flush.console()
  }  
}

### Execution of the program ###
start <- Sys.time()
z <- clusterApplyLB(cl = cl, v, analyse)
end <- Sys.time()
t.parallel <- difftime(end, start, units = "sec")
stopCluster(cl)

