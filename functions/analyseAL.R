########################################################################################
### Global function for the selection of biomarkers belonging to prespecified groups ###
########################################################################################

########################  PARAMETERS  ########################
############################################################################################
# - data        Dataset                                                                    #
# - method :    Variable weighting approaches                                              #   
#   - AC          The mean of the absolute values of the univariate regression coefficients#
#   - PCA         Principal component analysis by group of biomarkers                      #
#   - Lasso+PCA   Full Lasso + PCA                                                         #
#   - SW          Single Wald test                                                         #  
#   - MSW         Max single Wald test                                                     #  
#   - MSW*SW      Max single Wald test * Single Wald test                                  #
#   - ASW         Average single Wald test                                                 #
#   - ASW*SW      Average single Wald test * Single Wald test                              # 
# - paths       A vector indicating group membership for each biomaker                     #
# - nfolds      Number of folds. Default value is 5                                        #
# - thresh      Convergence threshold for coordinate descent. Default value is 1e-16       #
############################################################################################

### Working directory ###
#setwd("...")

### R packages ###
library(dplyr)
library(survival)
library(Matrix)
library(glmnet)
library(corpcor)
library(parallel) 

### Set a seed ###
#set.seed(...)

 
### Number of PC-cores used ###
pCores <- #80

### Definition of the cluster via 'makeCluster' ###
cl <- makeCluster(pCores, outfile ="OUTanalyseAL.txt")

### Number of replications ###
r <- 500

### Example of paths ###
paths = c(rep(1:20, each = 50))

v <- paste0('./datas/scen',rep(1:8, each =r),'/simdata-',sprintf("%03d",c(1:r)),'.Rdata')

analyse <- function(character){
  load(character)
  source("./functions/AL.R")
  respath <- gsub('(datas|data-)', 'resultsAL', character)
  
  if (isTRUE(file.exists(respath))){
    print('Found')
    flush.console()
  }else{
    methods = as.list(c("MC",
                        "PCA",
                        "Lasso+PCA",
                        "SW",
                        "MSW",
                        "MSW*SW",
                        "ASW",
                        "ASW*SW",
                        ))
    res <- lapply(methods,  FUN=function(X) AL(data[1:500,], method = X, paths = paths, nfolds = 5, thresh = 1e-16))          

    flush.console()
    save(res, file = respath)
    flush.console()
  } 
}

### Execution of the program ###
start <- Sys.time()
z <- clusterApplyLB(cl = cl, v, analyse)
end <- Sys.time()
t.parallel <- difftime(end, start, units = "sec")
stopCluster(cl)

