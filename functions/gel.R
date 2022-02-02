### Working directory ###
#setwd("...")

### R packages ###
library(grpreg)
 

### Set a seed ###
#set.seed(...)


### Biomarker groups ###
group <- c(rep(1:9, each = 25), rep(10:20, each = 75))

cv.gel <- function(data,
                   # seed,
                   group){
  
  biom <- grep("bm", names(data))
  nbiom <- length(biom)
  X = as.matrix(data[,nbiom])
  y = as.matrix(data[,c('time','status')])

  cv <-
    cv.grpsurv(
      X,
      y,
      # seed = seed,
      group = group,
      nfolds = 5,
      penalty = "gel",
      tau = 1 / 20,
      se = c('quick'),
      returnX = TRUE,
      returnY = TRUE,
      trace = TRUE
    )

  
  allresgel <- cv
  return(allresgel)

}
