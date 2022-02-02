### Working directory
setwd(" ")

#####################################################################
# R packages installation                                           #
#####################################################################
library(parallel)
library(SGL)

set.seed( )

#Index of biomarker groups
index <- c(rep(1:9, each = 25), rep(10:20, each = 75))

# Number of PC-cores used
pCores <- #80

# Definition of the cluster via 'makeCluster'
cl <- makeCluster(pCores, outfile = "OUTanalyseSGL.txt")

#Number of iterations
r <- 500

v <-
  paste0('./datas/scen',
         rep(1:8, each = r),
         '/simdata-',
         sprintf("%03d", c(1:r)),
         '.Rdata')


analyse <- function(character) {
  load(character)
  source("./functions/SGL.R")
  respath <- gsub('(datas|data-)', 'resultsSGL', character)
  
  if (isTRUE(file.exists(respath))) {
    #print('Found') #optional
    flush.console()
  } else{
    res <- cv.SGL(data[1:500, ], index)
    save(res, file = respath)
    flush.console()
  }
}

# Execution of the program
start <- Sys.time()
z <- clusterApplyLB(cl = cl, v, analyse)
end <- Sys.time()
t.parallel <- difftime(end, start, units = "sec")
stopCluster(cl)
