# @author Shaime Belhechmi

#########################################################################

###     Simulation function for survival data with biomarkers         ###

#########################################################################

########################  PARAMETERS  ########################
#########################################################################
# - N          the number of patients                                   #
# - p          the number of biomarkers                                 #
# - l          the number of groups/pathway of biomarkers               #  
# - p.active   the number of prognostic biomarkers                      #
# - l.active   the number of prognostic groups of biomarkers            #
# - w.act      the identifier of prognostic biomarkers                  #
# - t.prob     the treatment assignment probability                     #
# - m0         the baseline median survival times in years              #
# - alpha      the treatment effect                                     #
# - beta       the effect of prognostic biomarkers                      #
# - b.corr     the correlation within bm blocks                         #
# - b.corr.by  the size of the blocks of correlated bm's                #
# - wei.shape  the shape parameter of the Weibull distribution          #
# - recr       the recruitment period durations                         #
# - fu         the follow-up period duration                            #
# - timefac    the multiplicative factor for times                      #
#              for instance fix it to 365.25 if the mst are expressed   # 
#              in years, but times have to be expressed in days         #
#########################################################################

### R packages ###
library(MASS) 

simdata <- function(N, p, l, p.active, l.active, w.act, t.prob, m0, alpha, beta, b.corr,
                    b.corr.by, wei.shape, recr, fu, timefactor){   
  
###################################################################### 
### Parameters, checks and manipulation
  if(missing(p.active)) p.active <- 0
  if(missing(l.active)) l.active <- 0 
  if(missing(alpha)) alpha <- 0
  if(missing(beta)) beta <- NULL
  
### Sample size and number of biomarkers
  if(N <= 0 || p <= 0)
    stop("\nThe sample size 'n' and the number of biomarkers 'p' must be positive.")
  
### Active biomarkers
  if(p.active < 0 || l.active < 0)
    stop("\nThe number of active biomarkers 'p.active' and and the number of 
         active pathways 'l.active' must be positive or null.")
  
  if(p.active > p)
    stop("\nThe number of active biomarkers is larger than the overall number of biomarkers 'p'.")
  
  p.active <- round(p.active, 0)
  l.active <- round(l.active, 0)
  
### Effect of active biomarkers
  if(is.null(beta)){
    if(p.active > 0)
      stop("\nThe biomarkers' main effects 'beta' is missing, with no default.")
  }else{
    if(p.active > 0){
      if(length(beta) == 1) beta <- rep(beta, p.active)
      if(length(beta) > 1 && length(beta) < p.active){
        if(!length(beta) == 2)
          warning(paste0("\nThe length of the parameter 'beta' is not 1, 2 or ", p.active, " ('p.active'). ",
                         "Effects are randomly chosen within the range of the given list."))
        beta <- runif(n = p.active, min = min(beta), max = max(beta))
      }
      if(length(beta) > p.active)
        stop("\nThe biomarkers main effects 'beta' are larger than 'p.active'.")
    }else{
      warning("\nThe biomarkers'main effects 'beta' should be not taking into account as 'p.active' = 0.")
    }
  }
  
### Biomarkers' correlation
  if (!(b.corr >= 0 & b.corr < 1))
    stop("\nThe biomarkers correlation parameter 'b.corr' must be in [0, 1[.")
  
  if (p %% b.corr.by > 0)
    stop(paste0("\nThe size of the blocks of correlated biomarkers ",
                "'b.corr.by' must be a divisor of 'p', ",
                "the total number of biomarkers.\n"))

### Treatment assignment probability
  if(t.prob < 0 || t.prob > 1)
    stop("\nThe treatment assignment probability 't.prob' must be between 0 and 1.")

### Baseline median survival time
  if(m0 <= 0)
    stop("\nThe baseline median survival time 'm0' must be positive.")

### Weibull shape parameter
  if(wei.shape <= 0)
    stop("\nThe shape parameter 'wei.shape' must be positive.")
  
### Follow-up and recruitment period
  if(fu < 0 || recr < 0)
    stop("\nThe follow-up 'fu' and the recruitment period 'recr' must be positive.")

### Correlation
  if (length(as.vector(b.corr)) != 1)
    stop(paste0("\nThe biomarkers correlation parameter 'b.corr' ",
                "must be a scalar.\n"))
  
### Multiplicative factor for time
  if(timefactor <= 0)
    stop("\nThe multiplicate factor for time 'timefactor' must be positive.")
  m0 <- m0 * timefactor
######################################################################
  
### Data generation
### Generation of the treatment assignment
  treat = rbinom(N, 1, t.prob) - .5
  data <- data.frame(treat = rbinom(N, 1, t.prob) - .5)
  
### Correlation matrix
  n.blocks <- p %/% b.corr.by
  covMat <- diag(n.blocks) %x% matrix(b.corr^abs(matrix(1:b.corr.by, b.corr.by, b.corr.by, byrow = TRUE) - 
                                                   matrix(1:b.corr.by, b.corr.by, b.corr.by)), b.corr.by, b.corr.by)
  diag(covMat) <- 1
  
  
### Generation of the biomarkers
  data <- cbind(data, mvrnorm(N, rep(0, p), Sigma = covMat))
  rm(covMat)
  colnames(data)[2 :ncol(data)] <- sprintf("bm%03d",c(1:p))
  n.paths <- p/l
  
### Identifiant des pathways
  id.path <- sample(x = rep(1:l,n.paths), size = p, replace = FALSE)
  
### Effects of biomarkers
  HR <- rep(1, nrow(data))
  
  if(p.active == 0){
    wact <- NULL
  }else{
    wact <- sprintf("bm%03d", w.act)
    HR <- as.numeric(exp(as.matrix(data[, wact, drop=FALSE]) %*% matrix(beta, p.active)))
  }
 
### Derived Weibull scale parameters (parameter b = lambda^(-1/rho) in rweibull)
  t.effect <- alpha
  alpha <- rep(alpha, nrow(data))
  alpha[which(data$treat == -0.5)] <- 0
  wei.scale <- (m0 * (exp(alpha) * log(2) * HR) ^ (-1 / wei.shape))
 
### Event times
  data$time <- Vectorize(rweibull)(n = 1, shape = wei.shape, scale = wei.scale)

### Censoring times
  cTimes <- fu + recr * runif(N)
  data$status <- as.numeric(data$time <= cTimes)
  data$time <- pmin(data$time, cTimes)
  data <- data[,-1] #enlever la colonne traitement
  
### Appending attributes
  attributes(data) <- append(attributes(data), list(
    biomarkers = list(
      p.active = p.active,
      which.active = wact,
      beta.active = beta,
      l.active = l.active,
      correlation = b.corr,
      corr.block.size = b.corr.by),
    treatment = list(
      treatment.probability = t.prob,
      treatment.effect = t.effect,
      weibull.parameters = list(
        shape = wei.shape,
        scale = m0 ^ (-wei.shape) * log(2)),
      censoring = c("Recruitment time" = recr, 
                    "Minimum follow-up time" = fu),
      timefactor = timefactor)))
  return(data) 
}
