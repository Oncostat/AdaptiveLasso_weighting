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

### R packages ### 
# library(dplyr)
# library(survival)
# library(Matrix)
# library(glmnet)
# library(corpcor)
# library(parallel) 
 
AL <- function(data,
               paths,
               method = c(
                 "AC",        
                 "PCA",       
                 "Lasso+PCA", 
                 "SW",       
                 "MSW",       
                 "MSW*SW",    
                 "ASW",      
                 "ASW*SW",    
               ),
               nfolds = 5,
               thresh = 1e-16) {

  ### Control checking and data manipulation ###
  
  # Control of the data set
  if(missing(data))
    stop("\n 'data' is missing.")
  
  if (missing(paths))
    stop("Need to specify the pathways")
  else if (length(paths) <= 1)
    stop("Need to specify the 2 or more pathways")
  else
    path <- as.numeric(paths)
  
  # Control of the selected methods
  if (missing(method))
    stop("Need to specify the weight calculation method")
  else
    method <- match.arg(method)
  
  
  if (missing(nfolds))
    nfolds <- 5
  else
    nfolds <- as.numeric(nfolds)
  
  if (missing(thresh))
    thresh <- 1e-16
  else
    thresh <- as.numeric(thresh)
  
  
  biom <- grep("bm", names(data))
  nbiom <- length(biom)
  dataBiom <- data[, biom] 
  Y <- cbind(time = data$time, status = data$status)
  X <- as.matrix(data[, biom])
  nGp <- length(unique(paths))

  datai <- NULL
  for (i in 1:nGp) {
    datai[[i]] <-  dataBiom[, which(paths == i)]
  }
  
  ### Estimation of the weights for the adaptive lasso ###
  
  # Function for calculating the univariate coefficients beta
  Fcox <- function(i) {
    s <- summary(coxph(Surv(time, status) ~ data[, i], data = data))
    Beta <- log(s$conf.int[[1]])
    return(abs(Beta))
  }
  Betas = sapply(c(1:nbiom), Fcox)
  
  # Function for calculating the univariate coefficients of the 1st PCA
  FcoxPCA <- function(X) {  
    s <- summary(coxph(Surv(time, status) ~ X, data = data))
    Beta <- log(s$conf.int[[1]])
    return(abs(Beta))
  }
  
  # Function for calculating the Wald test 
  Fwald <- function(i) {
    Wtest <- coxph(Surv(time, status) ~ data[, i], data = data)$wald.test
    return(Wtest)
  }
  
  ### Weighting methods ###
  
  # AC
  if (method == "AC") {
    tab <- data.frame(x = factor(paths), y = Betas)
    w <- 1 / with(tab, tapply(y, x, mean)) 
    weights <- as.numeric(w[paths])
  }
  
  # PCA
  else if (method == "PCA") {
    ncomp <- 1
    UDV <- list(NULL)
    Fsvd <- function(X) {
      fast.svd(t(t(X) - colMeans(X)))
    }
    UDV <- lapply(datai, Fsvd)
    Fpca <- function(i) {
      pca <-
        as.matrix(UDV[[i]]$u[, 1:ncomp, drop = FALSE] %*% diag(UDV[[i]]$d[1:ncomp], ncomp))
      colnames(pca) <- gsub(" ", "", format(paste0("PCA1_", i)))
      return(pca)
    }
    PCA <- lapply(c(1:nGp), Fpca)
    w.pca <- 1 / sapply(PCA, FcoxPCA)
    #colnames(w.pca) <- paste0("W.PCA", 1:nGp) ### FACULTATIVE
    weights <- as.numeric(w.pca[paths])
  }

  # Lasso+PCA
  else if (method == "Lasso+PCA") {
    cvL <- cv.glmnet(
      X,
      Y,
      family = "cox",
      alpha = 1,
      nfolds = nfolds,
      grouped = TRUE,
      standardize = FALSE,
      thresh = thresh
    )
    
    fitL <- glmnet(
      X,
      Y,
      family = "cox",
      alpha = 1,
      lambda = cvL$lambda.min * ((nfolds - 1) / nfolds),
      standardize = FALSE,
      thresh = thresh
    )
    
    if (purrr::is_empty(coef(fitL)@Dimnames[[1]][which(as.vector(coef(fitL) != 0))])) {
      biom <- "NULL"
      weights <- rep(1, nbiom)
    } else{
      biom <- coef(fitL)@Dimnames[[1]][which(as.vector(coef(fitL) != 0))]
      DATA <- cbind(biom = colnames(dataBiom), path = paths)
      dataL <-
        merge(
          x = as.data.frame(DATA) ,
          y = as.data.frame(biom),
          by = "biom",
          all.y = TRUE
        )

      Fdata <- function(i) {
        j <- which(dataL$path == unique(dataL$path)[i])
        dataj <- data.frame(data[, as.character(dataL$biom[j])])
        colnames(dataj) <- as.character(dataL$biom[j])
        return(dataj)
      }
      DataL <- lapply(c(1:length(unique(dataL$path))), Fdata)
     
      ncomp <- 1
      UDV <- list(NULL)
      for(i in 1:length(DataL)){
        if(length(DataL[[i]])==1){ 
          UDV[[i]] <- DataL[[i]] 
        }else{
          UDV[[i]] <- fast.svd(t(t(DataL[[i]]) - colMeans(DataL[[i]])))
        }
      }
      
      Fpca <- function(i) {
        pca <-
          as.matrix(UDV[[i]]$u[, 1:ncomp, drop = FALSE] %*% diag(UDV[[i]]$d[1:ncomp], ncomp))
        colnames(pca) <- gsub(" ", "", format(paste0("PCA1_", i)))
        return(pca)
      }
      
      PCA <- list(NULL)
      for (i in 1:length(DataL)) {
        if (length(UDV[[i]]) == 1) {
          PCA[[i]] <- as.matrix(UDV[[i]]) 
          names(PCA[[i]]) <- paste0("PCA1_", i)
        } else{
          PCA[[i]] <-  Fpca(i)
        }
      }
      w.pca <- 1 / sapply(PCA, FcoxPCA)
      if (length(unique(dataL$path)) == length(unique(paths))) {
        tab <- 
          data.frame(path = sort(as.numeric(unique(dataL$path))), w. = w.pca)
      } else{
        tab1 <-
          data.frame(path = sort(as.numeric(unique(dataL$path))), w. = w.pca)
        tab2 <-
          data.frame(path = unique(paths)[-as.numeric(unique(dataL$path))], w. = max(w.pca))
        tab <- rbind(tab1, tab2)
        tab <- tab[order(tab$path), ]
      }
      weights <- as.numeric(tab$w.[paths])
    }
  }
  
  # SW
  else if (method == "SW") {
    tab <-
      data.frame(x = factor(paths), y = sapply(c(1:nbiom), Fwald))
    w <- 1 / tab$y
    weights <- as.numeric(w)
  }
  
  # MSW
  else if (method == "MSW") {
    tab <-
      data.frame(x = factor(paths), y = sapply(c(1:nbiom), Fwald))
    w <- 1 / with(tab, tapply(y, x, max))
    weights <- as.numeric(w[paths])
  }
  
  # MSW*SW
  else if (method == "MSW*SW") {
    tab <-
      data.frame(x = factor(paths), y = sapply(c(1:nbiom), Fwald))
    w <- 1 / with(tab, tapply(y, x, max))
    tab2 <- 1/tab[ ,2] 
    weights <- (as.numeric(w[paths]))*tab2
  }
  
  # ASW
  else if (method == "ASW") {
    tab <-
      data.frame(x = factor(paths), y = sapply(c(1:nbiom), Fwald))
    w <- 1 / with(tab, tapply(y, x, mean))
    weights <- (as.numeric(w[paths]))
  }
  
  # ASW*SW
  else if (method == "ASW*SW") {
    tab <-
      data.frame(x = factor(paths), y = sapply(c(1:nbiom), Fwald))
    w <- 1 / with(tab, tapply(y, x, mean))
    tab2 <- 1/tab[ ,2] 
    weights <- (as.numeric(w[paths]))*tab2
  }
  else {
    print(
      "the method must be one of (AC, PCA, Lasso+PCA, SW, MSW, MSW*SW, ASW, ASW*SW)"
    )
}
  
  
  ### Cross-validation to estimate the optimal lambda ###
  
  cv <- cv.glmnet(
    x = X,
    y = Y,
    family = "cox",
    alpha = 1,
    nfolds = nfolds,
    grouped = TRUE,
    standardize = FALSE,
    thresh = thresh,
    penalty.factor = weights
  )
  
  ### Fit of the final model ###
  
  fit.AL <- glmnet(
    x = X,
    y = Y,
    family = "cox",
    alpha = 1,
    lambda = cv$lambda.min * ((nfolds - 1) / nfolds),
    standardize = FALSE,
    thresh = thresh,
    penalty.factor = weights
  )

  
  if (length(which(coef(fit.AL) != 0)) == 0) {
    coef.AL <- "0"
    beta.AL <- 0
  } else{
    coef.AL <- coef(fit.AL)@Dimnames[[1]][which(coef(fit.AL) != 0)]
    beta.AL <- coef(fit.AL)@x
  }
  
  allresAL <- list(coef.AL = coef.AL,
                   beta.AL = beta.AL,
                   Betasi = Betas)
  return(allresAL)
  
}
