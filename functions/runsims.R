############################

### Generation of data   ###

############################ 

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

### Working directory ###
#setwd("./R code Simulation article/")

### R packages ### 
library(MASS)

### Set a seed ###
#set.seed()

### Number of replications ###
r <- 500

### Scenario 1: Complete null scenario
for (i in 1:r) {
  if (file.exists(paste0("./datas/scen1/simdata-", sprintf("%03d", i), ".Rdata"))) {
    print("Found")
  } else{
    data <-
      simdata(
        N = 1000,
        p = 1000,
        l = 20,
        t.prob = 0.5,
        m0 = 8,
        alpha = -0.7,
        beta = 0,
        b.corr = 0.8,
        b.corr.by = 5,
        wei.shape = 1,
        recr = 4,
        fu = 2,
        timefactor = 1
      )
    save(data, file = paste0("./datas/scen1/simdata-", sprintf("%03d", i), ".Rdata"))
  }
}

### Scenario 2: null scenario, without group effects 
for (i in 1:r) {
  if (file.exists(paste0("./datas/scen2/simdata-", sprintf("%03d", i), ".Rdata"))) {
    print("Found")
  } else{
    data <-
      simdata(
        N = 1000,
        p = 1000,
        l = 20,
        p.active = 100,
        l.active = 20,
        w.act = c(1:5, 26:30, 51:55, 76:80, 
                  101:105,126:130, 151:155,
                  176:180, 201:205,226:230, 
                  251:255,326:330, 401:405,
                  476:480, 551:555,626:630, 
                  701:705,776:780, 851:855,
                  926:930),
        t.prob = 0.5,
        m0 = 8,
        alpha = -0.7,
        beta = c(log(0.65), log(0.75)), 
        b.corr = 0.8,
        b.corr.by = 5,
        wei.shape = 1,
        recr = 4,
        fu = 2,
        timefactor = 1
      )
    save(data, file = paste0("./datas/scen2/simdata-", sprintf("%03d", i), ".Rdata"))
  }
}

### Scenario 3: l = 1; q = 8
for (i in 1:r) {
  if (file.exists(paste0("./datas/scen3/simdata-", sprintf("%03d", i), ".Rdata"))) {
    print("Found")
  } else{
    data <-
      simdata(
        N = 1000,
        p = 1000,
        l = 20,
        p.active = 8,
        l.active = 1,
        w.act = c(251:258),
        t.prob = 0.5,
        m0 = 8,
        alpha = -0.7,
        beta = c(log(0.65), log(0.75)), 
        b.corr = 0.8,
        b.corr.by = 5,
        wei.shape = 1,
        recr = 4,
        fu = 2,
        timefactor = 1
      )
    save(data, file = paste0("./datas/scen3/simdata-", sprintf("%03d", i), ".Rdata"))
  }
}

### Scenario 4: l = 2; q = 8
for (i in 1:r) {
  if (file.exists(paste0("./datas/scen4/simdata-", sprintf("%03d", i), ".Rdata"))) {
    print("Found")
  } else{
    data <-
      simdata(
        N = 1000,
        p = 1000,
        l = 20,
        p.active = 8,
        l.active = 2,
        w.act = c(1:4,251:254),
        t.prob = 0.5,
        m0 = 8,
        alpha = -0.7,
        beta = c(log(0.65), log(0.75)), 
        b.corr = 0.8,
        b.corr.by = 5,
        wei.shape = 1,
        recr = 4,
        fu = 2,
        timefactor = 1
      )
    save(data, file = paste0("./datas/scen4/simdata-", sprintf("%03d", i), ".Rdata"))
  }
}

### Scenario 5: l = 1; q = 32
for (i in 1:r) {
  if (file.exists(paste0("./datas/scen5/simdata-", sprintf("%03d", i), ".Rdata"))) {
    print("Found")
  } else{
    data <-
      simdata(
        N = 1000,
        p = 1000,
        l = 20,
        p.active = 32,
        l.active = 1,
        w.act = c(251:282),
        t.prob = 0.5,
        m0 = 8,
        alpha = -0.7,
        beta = c(log(0.85), log(0.95)), 
        b.corr = 0.8,
        b.corr.by = 5,
        wei.shape = 1,
        recr = 4,
        fu = 2,
        timefactor = 1
      )
    save(data, file = paste0("./datas/scen5/simdata-", sprintf("%03d", i), ".Rdata"))
  }
}

### Scenario 6: l = 2; q = 32
for (i in 1:r) {
  if (file.exists(paste0("./datas/scen6/simdata-", sprintf("%03d", i), ".Rdata"))) {
    print("Found")
  } else{
    data <-
      simdata(
        N = 1000,
        p = 1000,
        l = 20,
        p.active = 32,
        l.active = 2,
        w.act = c(1:16,251:266),
        t.prob = 0.5,
        m0 = 8,
        alpha = -0.7,
        beta = c(log(0.85), log(0.95)), 
        b.corr = 0.8,
        b.corr.by = 5,
        wei.shape = 1,
        recr = 4,
        fu = 2,
        timefactor = 1
      )
    save(data, file = paste0("./datas/scen6/simdata-", sprintf("%03d", i), ".Rdata"))
  }
}

### Scenario 7: l = 2; q = 16
for (i in 1:r) {
  if (file.exists(paste0("./datas/scen7/simdata-", sprintf("%03d", i), ".Rdata"))) {
    print("Found")
  } else{
    data <-
      simdata(
        N = 1000,
        p = 1000,
        l = 20,
        p.active = 16,
        l.active = 2,
        w.act = c(1:8,251:258),
        t.prob = 0.5,
        m0 = 8,
        alpha = -0.7,
        beta = c(log(0.85), log(0.95)), 
        b.corr = 0.8,
        b.corr.by = 5,
        wei.shape = 1,
        recr = 4,
        fu = 2,
        timefactor = 1
      )
    save(data, file = paste0("./datas/scen7/simdata-", sprintf("%03d", i), ".Rdata"))
  }
}

#### Scenario 8: l = 2; q = 8
for (i in 1:r) {
  if (file.exists(paste0("./datas/scen8/simdata-", sprintf("%03d", i), ".Rdata"))) {
    print("Found")
  } else{
    data <-
      simdata(
        N = 1000,
        p = 1000,
        l = 20,
        p.active = 8,
        l.active = 2,
        w.act = c(1:4,251:254),
        t.prob = 0.5,
        m0 = 8,
        alpha = -0.7,
        beta = c(log(0.85), log(0.95)), 
        b.corr = 0.8,
        b.corr.by = 5,
        wei.shape = 1,
        recr = 4,
        fu = 2,
        timefactor = 1
      )
    save(data, file = paste0("./datas/scen8/simdata-", sprintf("%03d", i), ".Rdata"))
  }
}


