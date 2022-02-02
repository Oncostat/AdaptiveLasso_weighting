#####################################################################
# R packages installation                                           #
#####################################################################
library(ipflasso)

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

ipf.lasso <- function(data, pflist, nfolds) {
  biom <- grep("bm", names(data))
  nbiom <- length(biom)
  dataBiom <- data[, biom]
  Y <- cbind(time = data$time, status = data$status)
  X <- as.matrix(data[, biom])
  start <- Sys.time()
  cv.ipf <-
    cvr2.ipflasso(
      X = X,
      Y = Y,
      family = "cox",
      type.measure = "deviance",
      standardize = FALSE,
      blocks = list(
        block1 = 1:25,
        block2 = 26:50,
        block3 = 51:75,
        block4 = 76:100,
        block5 = 101:125,
        block6 = 126:150,
        block7 = 151:175,
        block8 = 176:200,
        block9 = 201:225,
        block10 = 226:250,
        block11 = 251:325,
        block12 = 326:400,
        block13 = 401:475,
        block14 = 476:550,
        block15 = 551:625,
        block16 = 626:700,
        block17 = 701:775,
        block18 = 776:850,
        block19 = 851:925,
        block20 = 926:1000
      ),
      pflist = pflist,
      nfolds = nfolds,
      ncv = 10
    )
  
  end <- Sys.time()
  t <- difftime(end, start, units = "sec")
  print(which(!(cv.ipf$coeff[-1, cv.ipf$ind.bestlambda - 1] == 0)))
  if (length(which(!(cv.ipf$coeff[-1, cv.ipf$ind.bestlambda - 1] == 0))) == 0) {
    coef.ipf <- "0"
    beta.ipf <- 0
  } else {
    coef.ipf <-
      attr(which(!(cv.ipf$coeff[-1, cv.ipf$ind.bestlambda - 1] == 0)), "names")
    beta.ipf <-
      unname(cv.ipf$coeff[-1, cv.ipf$ind.bestlambda - 1] [unname(which(!(cv.ipf$coeff[-1, cv.ipf$ind.bestlambda -
                                                                                        1] == 0)))])
  }
  
  allresipf <- list(coef.ipf = coef.ipf,
                    beta.ipf = beta.ipf)
  print("Done")
  return(allresipf)
  
}
