#####################################################################
# R packages installation                                           #
#####################################################################
library(SGL)

set.seed( )

#Index of biomarker groups
# index <- c(rep(1:9, each = 25), rep(10:20, each = 75))

cv.SGL <- function(data, index) {
  biom <- grep("bm", names(data))
  
  X <- as.matrix(data[, biom])
  datas <- list(x = X,
                time = data$time,
                status = data$status)
  
  start <- Sys.time()
  cv <-
    cvSGL(
      data = datas,
      index = index,
      type = "cox",
      maxit = 1000,
      thresh = 0.001,
      min.frac = 0.05,
      nlam = 20,
      gamma = 0.8,
      standardize = FALSE,
      verbose = FALSE,
      step = 1,
      reset = 10,
      alpha = 0.95,
      lambdas = NULL
    )
  print(which((cv$fit$beta[, which(cv$lldiff == min(cv$lldiff))] != 0)))
  if (length(which((cv$fit$beta[, which(cv$lldiff == min(cv$lldiff))] != 0))) == 0) {
    coef.sgl <- "0"
    beta.sgl <- 0
  } else {
    coef.sgl <-
      sprintf("bm%03d", which((cv$fit$beta[, which(cv$lldiff == min(cv$lldiff))] != 0)))
    beta.sgl <-
      cv$fit$beta[, which(cv$lldiff == min(cv$lldiff))][which((cv$fit$beta[, which(cv$lldiff == min(cv$lldiff))] != 0))]
  }
  
  allressgl <- list(coef.sgl = coef.sgl,
                    beta.sgl = beta.sgl)
  print("Done")
  return(allressgl)
  
  
}