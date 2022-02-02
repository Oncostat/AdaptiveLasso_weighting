# @author Shaime Belhechmi

### Working directory ###
#setwd("...")

### R packages ###
library(grpreg)

### Set a seed ###
#set.seed(...)


### Biomarker groups ###
#group <- c(rep(1:9, each = 25), rep(10:20, each = 75))

### Number of PC-cores used ###
pCores <- 80

### Definition of the cluster via 'makeCluster' ###
cl <- makeCluster(pCores, outfile = "OUTcv.grpregGEL.txt")

### Number of iterations ###
r <- 500

v <-
  paste0('./datas/scen',
         rep(1:8, each = r),
         '/simdata-',
         sprintf("%03d", c(1:r)),
         '.Rdata')

analyse <- function(character) {
  load(character)
  source("./functions/gel.R")
  respath <- gsub('(datas|simdata-)', 'resultsGEL', character)
  
  if (isTRUE(file.exists(respath))) {
    #print('Found') #optional
    flush.console()
  } else{
    res <- cv.gel(data[1:500, ], group)
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
