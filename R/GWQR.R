library(GWmodel)
library(xgboost)
library(geosphere)
library(ggplot2)
library(spgwr)
library(dplyr)
library(quantreg)
library(tidyr)
library(mgcv)
library(purrr)
library(readr)
library(tmap)
library(MASS)
library(sf)
library(sp)
library(parallel) 
library(spData)
library(spdep)
library(shapviz)
library(mlspatial)
library(foreign)
library(ranger)
library(Metrics)
library(patchwork)


install.packages("remotes")

remotes::install_github(
  repo = "GWToolsClub/GWTools",
  ref = "main",
  dependencies = NA,
  upgrade = "never",
  build_vignettes = FALSE
)

library(GWTools)

install.packages("GWTools")

Income <- read_csv("C:/Users/chao6/OneDrive/桌面/yilingsim/greece.csv")
gvars <- c('UnemrT01', 'PrSect01' , 'Foreig01')
gcoords <- Income[, c("X", "Y")]

gf_string <- paste("Income01 ~", paste(gvars, collapse = " + "))
gf <- as.formula(gf_string) 

boston <- read_csv("C:/Users/chao6/OneDrive/桌面/yilingsim/Data.csv")
bco <- read_csv("C:/Users/chao6/OneDrive/桌面/yilingsim/Coords.csv")

bcoordx <- bco$Xcoord
bcoordy <- bco$Ycoord

bovars <- c('CRIM', 'ZN' , 'INDUS', 'CHAS', 'NOX', 'RM', 'AGE', 'DIS', 
            'RAD', 'TAX', 'PTRATIO', 'B', 'LSTAT')


bof_string <- paste("CMEDV ~", paste(bovars, collapse = " + "))
bof <- as.formula(bof_string) 

atlantic <- read_csv("C:/Users/chao6/OneDrive/桌面/yilingsim/data_atlantic_1998_2012.csv")
avars <- c('POV', 'SMOK', 'PM25', 'NO2', 'SO2')

af_string <- paste("Rate ~", paste(avars, collapse = " + "))
af <- as.formula(af_string) 

obesity <- read_csv("C:/Users/chao6/OneDrive/桌面/yilingsim/Obesity.csv")
ocoordx <- obesity$Lon
ocoordy <- obesity$Lat
colnames(obesity) <- make.names(colnames(obesity))
ovars <- c("X..Black", "X..Ame.Indi.and.AK.Native", "X..Asian", 
           "X..Nati.Hawa.and.Paci.Island", "X..Hispanic.or.Latino", "X..male", "X..married",
           "X..age.18.29", "X..age.30.39", "X..age.40.49", "X..age.50.59", "X..age...60", 
           "X...highschool", "median.income", "X..unemployment", "X..below.poverty.line",            
           "X..food.stamp.SNAP", "median.value.units.built", "median.year.units.built",         
           "X..renter.occupied.housing.units", "population.density" )


of_string <- paste("obesity_rate ~", paste(ovars, collapse = " + "))
of <- as.formula(of_string) 

coords_income <- cbind(X = Income$X, Y = Income$Y)
coords_bco <- cbind(X = bco$Xcoord, Y = bco$Ycoord)
coords_atlantic <- cbind(X = atlantic$x, Y = atlantic$y)
coords_obesity <- cbind(Lon = obesity$Lon, Lat = obesity$Lat)

#######GWQR
bw_gwqr_cv <- function(formula, tau = 0.5, cordxy, distmethod = "euclidean",
                       kernel, adaptive = TRUE, data, locallin = FALSE,
                       bw_min = NA, bw_max = NA, parallel = TRUE, core = NULL) {
  start <- proc.time()
  
  # ----- kernel function ----
  kernel_weights <- function(d, h, nr, kernel, adaptive) {
    if (adaptive == TRUE) {
      rh <- round(h)
      sd <- sort(d)
      bw <- sd[rh]
    } else {
      bw <- h
    }
    
    # w <- rep(0, length(d))
    if (kernel == "global") {
      w <- rep(1, nr)
    } else if (kernel == "gaussian") { # /*gaussian kernel*/
      w <- exp((-0.5) * ((d / bw)**2))
    } else if (kernel == "exponential") {
      w <- exp(-d / bw) # /*exponential kernel*/
    } else if (kernel == "bisquare") {
      w <- (1 - (d / bw)^2)^2
      index <- which(d > bw)
      w[index] <- 0 # /*bisquare nearest neighbor*/
    } else if (kernel == "tricube") {
      w <- (1 - abs(d / bw)^3)^3
      index <- which(d > bw)
      w[index] <- 0
    } else if (kernel == "boxcar") {
      w <- ifelse(d <= bw, 1, 0)
    } else {
      stop("Unsupported kernel type. Choose from
           'global', 'gaussian', 'exponential',
           'bisquare', 'tricube' or 'boxcar'.")
    }
    w <- as.numeric(w)
    return(w)
  }
  
  # ---- distance matrix ----
  dist_mat <- gwdist(cordxy, method = distmethod)
  
  all_dist <- sort(dist_mat)[which(sort(dist_mat) != 0)]
  all_dist <- all_dist[!duplicated(all_dist)]
  
  if (adaptive == TRUE) {
    auto_bw_min <- floor(0.1 * nrow(data))
    auto_bw_max <- ceiling(0.9 * nrow(data))
  } else {
    auto_bw_min <- all_dist[max(1, floor(length(all_dist) * 0.1))]
    auto_bw_max <- all_dist[ceiling(length(all_dist) * 0.9)]
  }
  
  if (!is.na(bw_min) && !is.na(bw_max)) {
    if (bw_min > bw_max)
      stop("Error: bw_min cannot be greater than bw_max")
  } else if (!is.na(bw_min) && is.na(bw_max)) {
    bw_max <- auto_bw_max
    if (bw_min > bw_max)
      stop("Error: The computed bw_max is smaller than the provided bw_min.")
  } else if (is.na(bw_min) && !is.na(bw_max)) {
    bw_min <- auto_bw_min
    if (bw_min > bw_max)
      stop("Error: The computed bw_min is larger than the provided bw_max.")
  } else {
    bw_min <- auto_bw_min
    bw_max <- auto_bw_max
  }
  
  # ---- CV function ----
  cvh <- function(h) {
    # ---- Prepare data ----
    mat <- model.frame(formula, data = data)
    y <- mat[, 1]
    x <- as.matrix(model.matrix(formula, data = data))
    if (attr(terms(formula), "intercept")) x <- x[, -1] else x <- x
    nr <- length(y)
    xx <- cbind(rep(1, nr), x)
    nvar <- ncol(xx)
    cordx <- cordxy[, 1]
    cordy <- cordxy[, 2]
    
    # ---- estimate function ----
    estimate <- function(times) {
      beta <- c()
      muhat <- c()
      for (i in times) {
        d <- gwdist(cordxy, method = distmethod, target_idx = i)
        # leave-one-out
        cid <- which(d != 0)
        dd <- d[cid]
        
        w <- kernel_weights(dd, h, nr - 1, kernel, adaptive)
        
        regdata <- data.frame(data[cid, ], w)
        # for local linear estimation
        if (locallin == TRUE) {
          xstart <- cordx[i]
          ystart <- cordy[i]
          
          xd <- cordx - xstart
          yd <- cordy - ystart
          
          xdz <- xd * x
          ydz <- yd * x
          
          allx <- cbind(x, xd, yd, xdz, ydz)
          callx <- allx[cid, ]
          
          regdata <- data.frame(callx)
        }
        
        fit <- rq(formula, tau = tau, weights = w,
                  method = "fn", data = regdata)
        fit1 <- summary(fit,se="nid")
        bmat <- coef(fit1)[1:nvar,1]
        #bmat <- fit$coefficients[1:nvar]
        beta <- rbind(beta, bmat)
        
        yhat <- xx[i, ] %*% as.matrix(bmat[1:nvar])
        muhat <- rbind(muhat, yhat)
      }
      output <- list(beta = beta, muhat = muhat)
      return(output)
    }
    
    # ---- parallel ----
    if (parallel == TRUE) {
      if (is.null(core)) {
        core <- detectCores() - 2
      }
      
      cl <- makeCluster(core)
      clusterExport(cl, c(
        "y", "nr", "cordxy", "cordx", "cordy",
        "formula", "xx", "x", "nvar", "kernel", "tau",
        "kernel_weights", "gwdist", "data", "bw_min",
        "bw_max", "h", "adaptive", "locallin"
      ),
      envir = environment()
      )
      clusterEvalQ(cl, c(library(quantreg), library(MASS)))
      sp <- clusterSplit(cl, 1:nr)
      
      res <- clusterApply(cl = cl, sp, fun = estimate)
      # str(res)
      on.exit(stopCluster(cl), add = TRUE)
    } else {
      res <- lapply(split(1:nr, 1:nr), estimate)
    }
    
    beta <- c()
    muhat <- c()
    for (l in 1:length(res)) {
      beta <- rbind(beta, res[[l]]$beta)
      muhat <- rbind(muhat, res[[l]]$muhat)
    }
    
    error <- y - muhat
    wsum <- sum((error * (tau - ifelse(error < 0, 1, 0))))
    print(sprintf("bandwidth: %.4f, wsum value: %.4f", h, wsum))
    out <- list(y = y, x = x, beta = beta, muhat = muhat,
                error = error, wsum = wsum)
    return(out)
  }
  
  # /*Golden Section Search*/
  eps <- 1
  r <- (sqrt(5) - 1) / 2
  a0 <- bw_min
  b0 <- bw_max
  x1 <- r * a0 + (1 - r) * b0
  x2 <- (1 - r) * a0 + r * b0
  fx1 <- cvh(x1)$wsum
  # print(cbind(fx1))
  fx2 <- cvh(x2)$wsum
  
  it <- 0
  while ((b0 - a0) > eps) {
    it <- it + 1
    if (fx1 < fx2) {
      b0 <- x2
      x2 <- x1
      fx2 <- fx1
      x1 <- r * a0 + (1 - r) * b0
      fx1 <- cvh(x1)$wsum
    } else {
      a0 <- x1
      x1 <- x2
      fx1 <- fx2
      x2 <- (1 - r) * a0 + r * b0
      fx2 <- cvh(x2)$wsum
    }
  }
  h <- ifelse(fx1 <= fx2, a0, b0) #* 將fx1小於fx2的b0改成a0;
  # print(paste(time1,Sys.time()))
  # print(paste(a0,b0))
  # print(cvh(h)$wsum)
  bw.time <- proc.time() - start
  # print(bw.time)
  names(h) <- c("lambda")
  return(h)
}

gwqr <- function(formula, tau = 0.5, cordxy, distmethod = "euclidean",
                 h, kernel, adaptive = TRUE, data, locallin = FALSE,
                 parallel = TRUE, core = NULL) {
  start <- proc.time()
  
  if (!is.matrix(cordxy)) {
    stop("cordxy must be a matrix.")
  }
  
  # ----- kernel function ----
  kernel_weights <- function(d, h, nr, kernel, adaptive) {
    if (adaptive == TRUE) {
      rh <- round(h)
      sd <- sort(d)
      bw <- sd[rh]
    } else {
      bw <- h
    }
    
    # w <- rep(0, length(d))
    if (kernel == "global") {
      w <- rep(1, nr)
    } else if (kernel == "gaussian") { # /*gaussian kernel*/
      w <- exp((-0.5) * ((d / bw)**2))
      w <- as.numeric(w)
    } else if (kernel == "exponential") {
      w <- exp(-d / bw)
      w <- as.numeric(w) # /*exponential kernel*/
    } else if (kernel == "bisquare") {
      w <- (1 - (d / bw)^2)^2
      index <- which(d > bw)
      w[index] <- 0 # /*bisquare nearest neighbor*/
    } else if (kernel == "tricube") {
      w <- (1 - abs(d / bw)^3)^3
      index <- which(d > bw)
      w[index] <- 0
    } else if (kernel == "boxcar") {
      w <- ifelse(d <= bw, 1, 0)
    } else {
      stop("Unsupported kernel type. Choose from
           'global', 'gaussian', 'exponential',
           'bisquare', 'tricube' or 'boxcar'.")
    }
    w <- as.numeric(w)
    return(w)
  }
  
  # ---- Prepare data ----
  mat <- model.frame(formula, data = data)
  y <- mat[, 1]
  x <- as.matrix(model.matrix(formula, data = data)) # designed x matrix including intercept (1's)
  ox <- x[, -1] # designed x matrix excluding intercept (1's)
  nvar <- ncol(x)
  nr <- length(y)
  cordx <- cordxy[, 1]
  cordy <- cordxy[, 2]
  
  
  # ---- estimate function ----
  estimate <- function(times) {
    beta <- c()
    muhat <- c()
    std <- c()
    tvalue <- c()
    for (i in times) {
      d <- gwdist(cordxy, method = distmethod, target_idx = i)
      w <- kernel_weights(d, h, nr, kernel, adaptive)
      w[i] <- 0
      regdata <- data.frame(ox)
      
      # for local linear estimation
      if (locallin == TRUE) {
        xstart <- cordx[i]
        ystart <- cordy[i]
        
        xd <- cordx - xstart
        yd <- cordy - ystart
        
        xdz <- xd * ox
        ydz <- yd * ox
        
        allx <- cbind(ox, xd, yd, xdz, ydz)
        
        regdata <- data.frame(allx)
      }
      fit <- quantreg::rq( y ~ ., tau = tau, weights = w,
                           method = "fn", data = regdata)
      fit1 <- quantreg::summary.rq(fit, se = "nid")
      bmat <- coef(fit1)[1:nvar, 1]
      semat <- coef(fit1)[1:nvar, 2]
      tmat <- coef(fit1)[1:nvar, 3]
      
      std <- rbind(std, semat)
      tvalue <- rbind(tvalue, tmat)
      yhat <- x[i, ] %*% as.matrix(bmat[1:nvar])
      beta <- rbind(beta, bmat)
      muhat <- rbind(muhat, yhat)
      # print(muhat)
    }
    # colnames(beta)[ncol(beta)]<- c("lambda")
    output <- list(beta = beta, std = std, muhat = muhat, tvalue = tvalue)
    return(output)
  }
  
  if (parallel == TRUE) {
    if (is.null(core)) {
      core <- parallel::detectCores() - 2
    }
    
    cl <- parallel::makeCluster(core)
    parallel::clusterExport(cl, c(
      "nr", "nvar", "tau", "x", "formula", "data", "y", "cordxy", "kernel", "h",
      "kernel_weights", "gwdist"
    ), envir = environment())
    parallel::clusterEvalQ(cl, c(library(quantreg), library(MASS)))
    sp <- parallel::clusterSplit(cl, 1:nr)
    
    res <- parallel::clusterApply(cl = cl, sp, fun = estimate)
    # str(res)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    res <- lapply(split(1:nr, 1:nr), estimate)
  }
  
  beta <- c()
  muhat <- c()
  std <- c()
  tvalue <- c()
  for (l in 1:length(res)) {
    beta <- rbind(beta, res[[l]]$beta)
    std <- rbind(std, res[[l]]$std)
    tvalue <- rbind(tvalue, res[[l]]$tvalue)
    muhat <- rbind(muhat, res[[l]]$muhat)
  }
  
  error <- y - muhat
  # loss function
  wsum <- sum((error * (tau - ifelse(error < 0, 1, 0))))
  out <- list(
    formula = formula, beta = beta, std = std, muhat = muhat,
    tvalue = tvalue, error = error, wsum = wsum,
    call = match.call()
  )
  class(out) <- "gwqr_model"
  return(out)
}

gwqr_parLOO <- function(formula,tau=0.5,geoid,cordx,cordy,h,kernel,indata,locallin=F,LOO=F,corenum){
  mat<-model.frame(formula, data = indata)
  y=mat[, 1]
  x=as.matrix(model.matrix(formula, data = indata))  #designed x matrix including intercept (1's)
  ox=x[,-1]  #designed x matrix excluding intercept (1's)
  nvar=ncol(x)
  nobs=length(y)
  
  beta=c(); muhat=c(); std=c(); tvalue=c()
  #parallel settings
  cl<-makeCluster(corenum)
  clusterExport(cl,c("nobs","nvar","tau","x","ox","formula","indata","locallin","geoid",
                     "y","cordx","cordy","kernel","h"), envir=environment())
  clusterEvalQ(cl,c(library(quantreg),library(MASS)))
  sp = clusterSplit(cl,1:nobs)
  
  #估計函數
  estimate<-function(times){
    beta=c(); muhat=c(); std=c(); tvalue=c()
    for (i in times){
      xstart=cordx[i]
      ystart=cordy[i]
      xd=cordx-xstart
      yd=cordy-ystart
      d<-sqrt(xd**2+yd**2)
      rh=round(h)
      sd=sort(d)
      disth=sd[rh]
      # kernel function
      if (kernel=="global") w=rep(1,nobs)
      if (kernel=="gaussian"){
        w=exp((-0.5)*((d/disth)**2))
        w=as.numeric(w)}# /*gaussian kernel*/
      
      if (kernel=="exponential"){
        w=exp(-d/disth)
        w=as.numeric(w)} #/*exponential kernel*/
      if (kernel=="bisquare"){
        #w=rep(0,nobs)
        w=(1-(d/disth)^2)^2
        index=which(d>disth)
        w[index]=0} #/*bisquare nearest neighbor*/
      
      if (!LOO) {w[i]=0}
      newdata = data.frame(ox)
      if (locallin==T) {  #for local linear estimation
        nvar1=3*nvar
        xdz=xd*ox; ydz=yd*ox
        allx = cbind(ox,xd,yd,xdz,ydz)
        newdata = data.frame(allx) }
      fit <-rq(y~.,tau=tau,weights=w,method="fn",data=newdata)
      fit1 <- summary(fit,se="nid")
      bmat <- coef(fit1)[1:nvar,1]
      #print(bmat)
      semat <- coef(fit1)[1:nvar,2]
      #print(semat)
      tmat <- coef(fit1)[1:nvar,3]
      
      std=rbind(std,semat)
      tvalue=rbind(tvalue,tmat)
      yhat = x[i,]%*%as.matrix(bmat[1:nvar])
      beta=rbind(beta,bmat)
      muhat=rbind(muhat,yhat)
      #print(muhat)
    }
    output=list(beta=beta,std=std,tvalue=tvalue,muhat=muhat)
    return(output)
  }
  
  #執行平行運算
  res<-clusterApply(cl=cl,sp,fun=estimate)
  #str(res)   #Investigate the structure of res
  #停止平行運算
  stopCluster(cl)
  #合併平行運算的output
  beta=c()
  muhat=c()
  std=c()
  tvalue=c()
  for (l in 1:length(res)){
    beta=rbind(beta,res[[l]]$beta)
    std=rbind(std,res[[l]]$std)
    tvalue=rbind(tvalue,res[[l]]$tvalue)
    muhat=rbind(muhat,res[[l]]$muhat)
  }
  error=y-muhat
  #loss function
  wsum=sum((error*(tau-ifelse(error<0,1,0))))
  out=list(y=y,x=x,geoid=geoid,cordx=cordx,cordy=cordy,beta=beta,std=std,tvalue=tvalue,muhat=muhat,error=error,wsum=wsum)  
  return(out)
}

#Greece
gqr10 <- bw_gwqr_cv(
  formula = gf, tau = 0.1,bw_min = 30, bw_max = 300,
  cordxy = coords_income, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = Income, locallin = FALSE,
  parallel = TRUE, core = 4)

ggwqr10 <- gwqr(
  formula = gf,
  tau = 0.1,
  cordxy = coords_income,
  distmethod = "euclidean",
  h = 200,
  kernel = "bisquare",
  adaptive = TRUE,
  data = Income,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

ggwqr10LOO <- gwqr_parLOO(
  formula = gf,
  tau = 0.1,
  geoid=1:325,
  cordx=Income$X,cordy=Income$Y,
  h = 200,
  kernel = "bisquare",
  indata = Income,
  locallin = FALSE,
#  LOO=T,
  core = 4
)
View(ggwqr10LOO)
# ggwqr10 <- gwqr(
#   formula = gf,
#   tau = 0.1,
#   cordxy = coords_income,
#   distmethod = "euclidean",
#   h = 174,
#   kernel = "bisquare",
#   adaptive = TRUE,
#   data = Income,
#   locallin = FALSE,
#   parallel = TRUE,
#   core = 4
# )

gqr25 <- bw_gwqr_cv(
  formula = gf, tau = 0.25,bw_min = 30, bw_max = 300,
  cordxy = coords_income, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = Income, locallin = FALSE,
  parallel = TRUE, core = 4)

ggwqr25 <- gwqr(
  formula = gf,
  tau = 0.25,
  cordxy = coords_income,
  distmethod = "euclidean",
  h = 130,
  kernel = "bisquare",
  adaptive = TRUE,
  data = Income,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)
# ggwqr25 <- gwqr(
#   formula = gf,
#   tau = 0.25,
#   cordxy = coords_income,
#   distmethod = "euclidean",
#   h = 127,
#   kernel = "bisquare",
#   adaptive = TRUE,
#   data = Income,
#   locallin = FALSE,
#   parallel = TRUE,
#   core = 4
# )

gqr50 <- bw_gwqr_cv(
  formula = gf, tau = 0.5,bw_min = 30, bw_max = 300,
  cordxy = coords_income, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = Income, locallin = FALSE,
  parallel = TRUE, core = 4)


ggwqr50 <- gwqr(
  formula = gf,
  tau = 0.5,
  cordxy = coords_income,
  distmethod = "euclidean",
  h = 124,
  kernel = "bisquare",
  adaptive = TRUE,
  data = Income,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)
# ggwqr50 <- gwqr(
#   formula = gf,
#   tau = 0.5,
#   cordxy = coords_income,
#   distmethod = "euclidean",
#   h = 134,
#   kernel = "bisquare",
#   adaptive = TRUE,
#   data = Income,
#   locallin = FALSE,
#   parallel = TRUE,
#   core = 4
# )

gqr75 <- bw_gwqr_cv(
  formula = gf, tau = 0.75,bw_min = 30, bw_max = 300,
  cordxy = coords_income, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = Income, locallin = FALSE,
  parallel = TRUE, core = 4)


ggwqr75 <- gwqr(
  formula = gf,
  tau = 0.75,
  cordxy = coords_income,
  distmethod = "euclidean",
  h = 111,
  kernel = "bisquare",
  adaptive = TRUE,
  data = Income,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)
# ggwqr75 <- gwqr(
#   formula = gf,
#   tau = 0.75,
#   cordxy = coords_income,
#   distmethod = "euclidean",
#   h = 44,
#   kernel = "bisquare",
#   adaptive = TRUE,
#   data = Income,
#   locallin = FALSE,
#   parallel = TRUE,
#   core = 4
# )

gqr90 <- bw_gwqr_cv(
  formula = gf, tau = 0.9,bw_min = 30, bw_max = 300,
  cordxy = coords_income, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = Income, locallin = FALSE,
  parallel = TRUE, core = 4)

ggwqr90 <- gwqr(
  formula = gf,
  tau = 0.9,
  cordxy = coords_income,
  distmethod = "euclidean",
  h = 156,
  kernel = "bisquare",
  adaptive = TRUE,
  data = Income,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)
# ggwqr90 <- gwqr(
#   formula = gf,
#   tau = 0.9,
#   cordxy = coords_income,
#   distmethod = "euclidean",
#   h = 229,
#   kernel = "bisquare",
#   adaptive = TRUE,
#   data = Income,
#   locallin = FALSE,
#   parallel = TRUE,
#   core = 4
# )


#Boston
bqr10 <- bw_gwqr_cv(
  formula = bof, tau = 0.1, bw_min = 350, bw_max = 460,
  cordxy = coords_bco, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = boston, locallin = FALSE,
  parallel = TRUE, core = 4)

bgwqr10loo <- gwqr(
  formula = bof,
  tau = 0.1,
  cordxy = coords_bco,
  distmethod = "euclidean",
  h = bqr10,
  kernel = "bisquare",
  adaptive = TRUE,
  data = boston,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

bqr25 <- bw_gwqr_cv(
  formula = bof, tau = 0.25,bw_min = 350, bw_max = 460,
  cordxy = coords_bco, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = boston, locallin = FALSE,
  parallel = TRUE, core = 4)

bgwqr25loo <- gwqr(
  formula = bof,
  tau = 0.25,
  cordxy = coords_bco,
  distmethod = "euclidean",
  h = bqr25,
  kernel = "bisquare",
  adaptive = TRUE,
  data = boston,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

bqr50 <- bw_gwqr_cv(
  formula = bof, tau = 0.5,bw_min = 350, bw_max = 460,
  cordxy = coords_bco, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = boston, locallin = FALSE,
  parallel = TRUE, core = 4)

bgwqr50loo <- gwqr(
  formula = bof,
  tau = 0.5,
  cordxy = coords_bco,
  distmethod = "euclidean",
  h = bqr50,
  kernel = "bisquare",
  adaptive = TRUE,
  data = boston,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

bqr75 <- bw_gwqr_cv(
  formula = bof, tau = 0.75,bw_min = 350, bw_max = 460,
  cordxy = coords_bco, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = boston, locallin = FALSE,
  parallel = TRUE, core = 4)

bgwqr75loo <- gwqr(
  formula = bof,
  tau = 0.75,
  cordxy = coords_bco,
  distmethod = "euclidean",
  h = bqr75,
  kernel = "bisquare",
  adaptive = TRUE,
  data = boston,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

bqr90 <- bw_gwqr_cv(
  formula = bof, tau = 0.9,bw_min = 350, bw_max = 460,
  cordxy = coords_bco, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = boston, locallin = FALSE,
  parallel = TRUE, core = 4)

bgwqr90loo <- gwqr(
  formula = bof,
  tau = 0.9,
  cordxy = coords_bco,
  distmethod = "euclidean",
  h = bqr90,
  kernel = "bisquare",
  adaptive = TRUE,
  data = boston,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

#Atlantic
aqr10 <- bw_gwqr_cv(
  formula = af, tau = 0.1,bw_min = 60, bw_max = 600,
  cordxy = coords_atlantic, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = atlantic, locallin = FALSE,
  parallel = TRUE, core = 4)

agwqr10loo <- gwqr(
  formula = af,
  tau = 0.1,
  cordxy = coords_atlantic,
  distmethod = "euclidean",
  h = aqr10,
  kernel = "bisquare",
  adaptive = TRUE,
  data = atlantic,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

aqr25 <- bw_gwqr_cv(
  formula = af, tau = 0.25,bw_min = 60, bw_max = 600,
  cordxy = coords_atlantic, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = atlantic, locallin = FALSE,
  parallel = TRUE, core = 4)

agwqr25loo <- gwqr(
  formula = af,
  tau = 0.25,
  cordxy = coords_atlantic,
  distmethod = "euclidean",
  h = aqr25,
  kernel = "bisquare",
  adaptive = TRUE,
  data = atlantic,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

aqr50 <- bw_gwqr_cv(
  formula = af, tau = 0.5,bw_min = 60, bw_max = 600,
  cordxy = coords_atlantic, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = atlantic, locallin = FALSE,
  parallel = TRUE, core = 4)

agwqr50loo <- gwqr(
  formula = af,
  tau = 0.5,
  cordxy = coords_atlantic,
  distmethod = "euclidean",
  h = aqr50,
  kernel = "bisquare",
  adaptive = TRUE,
  data = atlantic,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

aqr75 <- bw_gwqr_cv(
  formula = af, tau = 0.75,bw_min = 60, bw_max = 600,
  cordxy = coords_atlantic, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = atlantic, locallin = FALSE,
  parallel = TRUE, core = 4)

agwqr75loo <- gwqr(
  formula = af,
  tau = 0.75,
  cordxy = coords_atlantic,
  distmethod = "euclidean",
  h = aqr75,
  kernel = "bisquare",
  adaptive = TRUE,
  data = atlantic,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

aqr90 <- bw_gwqr_cv(
  formula = af, tau = 0.9,bw_min = 60, bw_max = 600,
  cordxy = coords_atlantic, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = atlantic, locallin = FALSE,
  parallel = TRUE, core = 4)

agwqr90loo <- gwqr(
  formula = af,
  tau = 0.9,
  cordxy = coords_atlantic,
  distmethod = "euclidean",
  h = aqr90,
  kernel = "bisquare",
  adaptive = TRUE,
  data = atlantic,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

#obesity
oqr10 <- bw_gwqr_cv(
  formula = of, tau = 0.1,bw_min = 350, bw_max = 1800,
  cordxy = coords_obesity, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = obesity, locallin = FALSE,
  parallel = TRUE, core = 4)

ogwqr10loo <- gwqr(
  formula = of,
  tau = 0.1,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 673,
  kernel = "bisquare",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

# ogwqr10LOO <- gwqr_parLOO(
#   formula = of,
#   tau = 0.1,
#   geoid=1:1995,
#   cordx=obesity$Lon,cordy=obesity$Lat,
#   h = 672,
#   kernel = "bisquare",
#   indata = obesity,
#   locallin = FALSE,
#   #  LOO=T,
#   core = 4
# )
# View(ogwqr10LOO)

oqr25 <- bw_gwqr_cv(
  formula = of, tau = 0.25,bw_min = 350, bw_max = 1800,
  cordxy = coords_obesity, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = obesity, locallin = FALSE,
  parallel = TRUE, core = 4)

ogwqr25loo <- gwqr(
  formula = of,
  tau = 0.25,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 610,
  kernel = "bisquare",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

oqr50 <- bw_gwqr_cv(
  formula = of, tau = 0.5,bw_min = 350, bw_max = 1800,
  cordxy = coords_obesity, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = obesity, locallin = FALSE,
  parallel = TRUE, core = 4)

ogwqr50loo <- gwqr(
  formula = of,
  tau = 0.5,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 606,
  kernel = "bisquare",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

ogwqr50t <- gwqr(
  formula = of,
  tau = 0.5,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 554,
  kernel = "bisquare",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

oqr75 <- bw_gwqr_cv(
  formula = of, tau = 0.75,bw_min = 350, bw_max = 1800,
  cordxy = coords_obesity, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = obesity, locallin = FALSE,
  parallel = TRUE, core = 4)

ogwqr75loo <- gwqr(
  formula = of,
  tau = 0.75,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 617,
  kernel = "bisquare",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

oqr90 <- bw_gwqr_cv(
  formula = of, tau = 0.9,bw_min = 350, bw_max = 1800,
  cordxy = coords_obesity, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = obesity, locallin = FALSE,
  parallel = TRUE, core = 4)

ogwqr90loo <- gwqr(
  formula = of,
  tau = 0.9,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 610,
  kernel = "bisquare",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

calc_gwqr_metrics <- function(gwqr_obj, tau, p = NULL, dataset_name = "Dataset") {
  yhat <- as.numeric(gwqr_obj$muhat)
  error <- as.numeric(gwqr_obj$error)
  y <- error + yhat
  n <- length(y)

  if (is.null(p)) {
    p <- ncol(as.matrix(gwqr_obj$beta))
  }
  
  pinball_point <- function(y, qhat, tau) {
    u <- y - qhat
    pmax(tau * u, (tau - 1) * u)
  }
  
  loss_gwqr <- pinball_point(y, yhat, tau)
  
  pinball_loss <- mean(loss_gwqr)
  sum_loss <- sum(loss_gwqr)
  
  wsum_from_obj <- as.numeric(gwqr_obj$wsum)
  
  baseline_pred <- as.numeric(quantile(y, probs = tau, names = FALSE))
  
  loss_base <- pinball_point(y, baseline_pred, tau)
  baseline_sum <- sum(loss_base)
  baseline_pinball_loss <- mean(loss_base)
  
  pseudo_r2 <- 1 - (sum_loss / baseline_sum)
  adj_pseudo_r2 <- 1 - ((sum_loss / (n - p)) / (baseline_sum / (n - 1)))
  
  out <- data.frame(
    Dataset = dataset_name,
    Method = "GWQR",
    tau = tau,
    n = n,
    p = p,
    Pinball_Loss = pinball_loss,
    Sum_Pinball_Loss = sum_loss,
    #Baseline_Pinball_Loss = baseline_pinball_loss,
    #Baseline_Sum_Loss = baseline_sum,
    Pseudo_R2 = pseudo_r2,
    #Adj_Pseudo_R2 = adj_pseudo_r2,
    wsum_from_gwqr_object = wsum_from_obj
    #check_wsum_difference = sum_loss - wsum_from_obj
  )
  
  return(out)
}
calc_qrtree_metrics <- function(gwqr_obj, tau, p = NULL, dataset_name = "Dataset") {
  yhat <- as.numeric(gwqr_obj$yhat_local_oob)
  error <- as.numeric(gwqr_obj$error_local_oob)
  y <- error + yhat
  n <- length(y)
  
  if (is.null(p)) {
    p <- 18
  }
  
  pinball_point <- function(y, qhat, tau) {
    u <- y - qhat
    pmax(tau * u, (tau - 1) * u)
  }
  
  loss_gwqr <- pinball_point(y, yhat, tau)
  
  pinball_loss <- mean(loss_gwqr)
  sum_loss <- sum(loss_gwqr)
  
  wsum_from_obj <- as.numeric(gwqr_obj$wsum_oob)
  
  baseline_pred <- as.numeric(quantile(y, probs = tau, names = FALSE))
  
  loss_base <- pinball_point(y, baseline_pred, tau)
  baseline_sum <- sum(loss_base)
  baseline_pinball_loss <- mean(loss_base)
  
  pseudo_r2 <- 1 - (sum_loss / baseline_sum)
  adj_pseudo_r2 <- 1 - ((sum_loss / (n - p)) / (baseline_sum / (n - 1)))
  
  out <- data.frame(
    Dataset = dataset_name,
    Method = "QRTree",
    tau = tau,
    n = n,
    p = p,
    Pinball_Loss = pinball_loss,
    Sum_Pinball_Loss = sum_loss,
    #Baseline_Pinball_Loss = baseline_pinball_loss,
    #Baseline_Sum_Loss = baseline_sum,
    Pseudo_R2 = pseudo_r2,
    #Adj_Pseudo_R2 = adj_pseudo_r2,
    wsum_from_gwqr_object = wsum_from_obj
    #check_wsum_difference = sum_loss - wsum_from_obj
  )
  
  return(out)
}

gwqr_metric <- calc_qrtree_metrics(
  gwqr_obj = rqxgb10,
  tau = 0.1,
  dataset_name = "R"
)
gwqr_metric


gwqr_metric <- calc_gwqr_metrics(
  gwqr_obj = rgwqr90loo,
  tau = 0.9,
  dataset_name = "R"
)
gwqr_metric


##算Moran's I
auto_pick_knn <- function(coords, y, k_grid = 4:12, longlat=TRUE) {
  out <- data.frame(k = c(), ncomp = c(), I = c(), Z = c(), p = c())
  for (k in k_grid) {
    nb0 <- knn2nb(knearneigh(coords, k = k, longlat = longlat))
    nb <- make.sym.nb(nb0) # 對稱化
    ncomp <- n.comp.nb(nb)$nc
    lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
    mt <- suppressWarnings(moran.test(y, lw, zero.policy = TRUE))
    out <- rbind(out, data.frame(
      k = k, ncomp = ncomp, I = mt$estimate[["Moran I statistic"]],Z = mt$statistic,p = mt$p.value))
  }
  out
}

res_knnf <- auto_pick_knn(coords_real, rgwqr90loo$error, k_grid = 1:300, longlat=FALSE)
res_knnf[res_knnf$ncomp == 1, ][which.max(res_knnf[res_knnf$ncomp == 1, ]$Z), ]
res_knnf[which.max(res_knnf$Z),]


#Greece
gqr10 <- gwqr(
  formula = gf,
  tau = 0.1,
  cordxy = coords_income,
  distmethod = "euclidean",
  h = 325,
  kernel = "global",
  adaptive = TRUE,
  data = Income,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

gqr25 <- gwqr(
  formula = gf,
  tau = 0.25,
  cordxy = coords_income,
  distmethod = "euclidean",
  h = 325,
  kernel = "global",
  adaptive = TRUE,
  data = Income,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

gqr50 <- gwqr(
  formula = gf,
  tau = 0.5,
  cordxy = coords_income,
  distmethod = "euclidean",
  h = 325,
  kernel = "global",
  adaptive = TRUE,
  data = Income,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

gqr75 <- gwqr(
  formula = gf,
  tau = 0.75,
  cordxy = coords_income,
  distmethod = "euclidean",
  h = 325,
  kernel = "global",
  adaptive = TRUE,
  data = Income,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

gqr90 <- gwqr(
  formula = gf,
  tau = 0.9,
  cordxy = coords_income,
  distmethod = "euclidean",
  h = 325,
  kernel = "global",
  adaptive = TRUE,
  data = Income,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)


#Boston

bqr10 <- gwqr(
  formula = bof,
  tau = 0.1,
  cordxy = coords_bco,
  distmethod = "euclidean",
  h = 506,
  kernel = "global",
  adaptive = TRUE,
  data = boston,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)
bqr25 <- gwqr(
  formula = bof,
  tau = 0.25,
  cordxy = coords_bco,
  distmethod = "euclidean",
  h = 506,
  kernel = "global",
  adaptive = TRUE,
  data = boston,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

bqr50 <- gwqr(
  formula = bof,
  tau = 0.5,
  cordxy = coords_bco,
  distmethod = "euclidean",
  h = 506,
  kernel = "global",
  adaptive = TRUE,
  data = boston,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)
bqr75 <- gwqr(
  formula = bof,
  tau = 0.75,
  cordxy = coords_bco,
  distmethod = "euclidean",
  h = 506,
  kernel = "global",
  adaptive = TRUE,
  data = boston,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

bqr90 <- gwqr(
  formula = bof,
  tau = 0.9,
  cordxy = coords_bco,
  distmethod = "euclidean",
  h = 506,
  kernel = "global",
  adaptive = TRUE,
  data = boston,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

#Atlantic
aqr10 <- gwqr(
  formula = af,
  tau = 0.1,
  cordxy = coords_atlantic,
  distmethod = "euclidean",
  h = 666,
  kernel = "global",
  adaptive = TRUE,
  data = atlantic,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

aqr25 <- gwqr(
  formula = af,
  tau = 0.25,
  cordxy = coords_atlantic,
  distmethod = "euclidean",
  h = 666,
  kernel = "global",
  adaptive = TRUE,
  data = atlantic,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

aqr50 <- gwqr(
  formula = af,
  tau = 0.5,
  cordxy = coords_atlantic,
  distmethod = "euclidean",
  h = 666,
  kernel = "global",
  adaptive = TRUE,
  data = atlantic,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

aqr75 <- gwqr(
  formula = af,
  tau = 0.75,
  cordxy = coords_atlantic,
  distmethod = "euclidean",
  h = 666,
  kernel = "global",
  adaptive = TRUE,
  data = atlantic,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

aqr90 <- gwqr(
  formula = af,
  tau = 0.9,
  cordxy = coords_atlantic,
  distmethod = "euclidean",
  h = 666,
  kernel = "global",
  adaptive = TRUE,
  data = atlantic,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

#obesity
oqr10 <- gwqr(
  formula = of,
  tau = 0.1,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 1995,
  kernel = "global",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

oqr25 <- gwqr(
  formula = of,
  tau = 0.25,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 1995,
  kernel = "global",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

oqr50 <- gwqr(
  formula = of,
  tau = 0.5,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 1995,
  kernel = "global",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

oqr75 <- gwqr(
  formula = of,
  tau = 0.75,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 1995,
  kernel = "global",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

oqr90 <- gwqr(
  formula = of,
  tau = 0.9,
  cordxy = coords_obesity,
  distmethod = "euclidean",
  h = 1995,
  kernel = "global",
  adaptive = TRUE,
  data = obesity,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)



data_dir <- "C:/Users/chao6/OneDrive/桌面/yilingsim/sim2_datasets" 
n_rep    <- 100       
tau_list <- c(0.1, 0.25, 0.5, 0.75, 0.9)
core_num <- 4

all_results <- list()

for (t_idx in seq_along(tau_list)) {
  tau_value <- tau_list[t_idx]
  qtrue_col <- paste0("V", t_idx)  
  
  metrics_rep <- data.frame(
    rep      = 1:n_rep,
    pinball  = NA_real_,
    r2       = NA_real_,
    adj_r2   = NA_real_,
    q_mse    = NA_real_
  )
  
  for (r in 1:n_rep) {
    rep_str    <- sprintf("%03d", r)
    sim_file   <- file.path(data_dir, paste0("sim2_rep", rep_str, ".csv"))
    qtrue_file <- file.path(data_dir, paste0("qtrue2_rep", rep_str, ".csv"))
    
    if (!file.exists(sim_file) || !file.exists(qtrue_file)) {
      warning(sprintf("tau=%.2f, rep %s 找不到檔案，略過", tau_value, rep_str))
      next
    }
    
    sim_data   <- read.csv(sim_file)
    qtrue_data <- read.csv(qtrue_file)
    
    cordxy <- as.matrix(sim_data[, c("u", "v")])
    y <- sim_data$y2
    n <- length(y)
    p_val  <- 6 
    
    fit <- gwqr(
      formula = y2 ~ x1 + x2 + x3 + x4 + x5,
      tau = tau_value,
      cordxy = cordxy,
      distmethod = "euclidean",
      h = n,        
      kernel = "global",
      adaptive = TRUE,
      data = sim_data,
      locallin = FALSE,
      parallel = TRUE,
      core = core_num
    )
    
    yhat_fix <- as.numeric(fit$muhat)
    qtrue    <- qtrue_data[[qtrue_col]]
    
    # ---- 指標計算 ----
    err_fix   <- y - yhat_fix
    pinball_fix <- mean(pmax(tau_value * err_fix, (tau_value - 1) * err_fix))
    sum_fix     <- sum(pmax(tau_value * err_fix, (tau_value - 1) * err_fix))
    
    baseline_pred <- as.numeric(quantile(y, probs = tau_value, names = FALSE))
    err_base      <- y - baseline_pred
    baseline_sum  <- sum(pmax(tau_value * err_base, (tau_value - 1) * err_base))
    
    r2_fix     <- 1 - (sum_fix / baseline_sum)
    adj_r2_fix <- 1 - ((sum_fix / (n - p_val)) / (baseline_sum / (n - 1)))
    q_mse_fix  <- mean((qtrue - yhat_fix)^2)
    
    metrics_rep$pinball[r] <- pinball_fix
    metrics_rep$r2[r]      <- r2_fix
    metrics_rep$adj_r2[r]  <- adj_r2_fix
    metrics_rep$q_mse[r]   <- q_mse_fix
    
    cat(sprintf("tau = %.2f | rep %s 完成\n", tau_value, rep_str))
  }
  
  all_results[[t_idx]] <- metrics_rep
}

summary_table <- do.call(rbind, lapply(seq_along(tau_list), function(t_idx) {
  m <- all_results[[t_idx]]
  data.frame(
    tau          = tau_list[t_idx],
    pinball_mean = mean(m$pinball, na.rm = TRUE),
    r2_mean = mean(m$r2, na.rm = TRUE),
    adj_r2_mean  = mean(m$adj_r2, na.rm = TRUE),
    q_mse_mean   = mean(m$q_mse, na.rm = TRUE),
    n_valid_rep  = sum(!is.na(m$pinball))
  )
}))

print(summary_table)

write.csv(summary_table, "gwqr_loo_summary.csv", row.names = FALSE)

##Polygon算Moran's I
moran_queen <- function(poly, y, style = "W", zero.policy = TRUE) {
  
  nb <- poly2nb(poly, queen = TRUE)
  ncomp <- n.comp.nb(nb)$nc
  n_no_neighbor <- sum(card(nb) == 0)
  
  lw <- nb2listw(nb, style = style, zero.policy = zero.policy)
  mt <- moran.test(y, lw, zero.policy = zero.policy)
  
  cat("=====================================================\n")
  cat(sprintf(" Queen 鄰接結果\n"))
  cat("=====================================================\n")
  cat(sprintf(" 連通元件數 (ncomp)     : %d\n", ncomp))
  cat(sprintf(" 無鄰居的孤立多邊形數量  : %d\n", n_no_neighbor))
  cat(sprintf(" Moran's I : %.4f\n", mt$estimate[["Moran I statistic"]]))
  cat(sprintf(" Z value : %.4f\n", mt$statistic))
  cat(sprintf(" p-value : %.4f\n", mt$p.value))
  cat("-----------------------------------------------------\n\n")
  
  list(nb = nb, listw = lw, moran = mt, ncomp = ncomp)
}


result_df <- data.frame(
  GEOID = as.character(obesity$Census.tract.code),
  error = ogwqr10loo$error
)
poly <- st_read("C:/Users/chao6/OneDrive/桌面/yilingsim/tl_2018_36_tract/tl_2018_36_tract.shp")
poly$GEOID <- as.character(poly$GEOID)

poly_merged <- poly %>%
  left_join(result_df %>% dplyr::select(GEOID, error), by = "GEOID")

cat("未合併到資料的 tract 數量:", sum(is.na(poly_merged$error)), "\n")
cat("原始 poly 筆數:", nrow(poly), " / mor 筆數:", nrow(result_df), "\n")
poly_clean <- poly_merged %>% filter(!is.na(error))
res_queen <- moran_queen(poly_clean, poly_clean$error)

#Greece Moran's I
result_df <- data.frame(
  X = as.numeric(Income$X),
  Y = as.numeric(Income$Y),
  error = ggwqr90$error
) %>%
  mutate(
    X_key = round(X, 6),
    Y_key = round(Y, 6)
  )

poly <- st_read("C:/Users/chao6/OneDrive/桌面/yilingsim/GR_Municipalities.shp",quiet = TRUE)
poly <- poly %>%
  st_make_valid() %>%
  mutate(
    X_key = round(as.numeric(X), 6),
    Y_key = round(as.numeric(Y), 6)
  )

poly_merged <- poly %>%
  left_join(
    result_df %>%
      dplyr::select(
        X_key,
        Y_key,
        error
      ),
    by = c("X_key", "Y_key")
  )

cat("未合併到資料的 municipality 數量：",sum(is.na(poly_merged$error)), "\n")
cat("原始 polygon 筆數：",nrow(poly)," / Moran 資料筆數：",nrow(result_df),"\n")
res_queen <- moran_queen(poly_merged,poly_merged$error)


##實價登錄資料
real113 <- read_csv("C:/Users/chao6/OneDrive/桌面/yilingsim/realprice_combined_113_allclean.csv")
real113$realy <- log(real113$總價元)
real113$liv <- real113$`建物現況格局-廳`
real113$bath <- real113$`建物現況格局-衛`
real113$ele <- ifelse(real113$電梯 == "有", 1, 0)
real113$area <- real113$建物移轉總面積平方公尺
real113$lon <- real113$經度
real113$lat <- real113$緯度
realvars <- c( 'ele', 'liv', 'bath', 'area', 'count_floor', 'house_age',
               'Dependency_ratio', 'BACHELOR_UP_rate', 'Old_ratio', 'popdense',
               'aqi_idw', 'daily_rain_idw', 't_max_idw',
               'count_med_500m', 'count_phar_500m', 'dist_to_mrt_meter', 'dist_to_rail_meter')

realcoords <- real113[, c("lon", "lat")]

rf_string <- paste("realy ~", paste(realvars, collapse = " + "))
rf <- as.formula(rf_string)
coords_real <- cbind(X = real113$lon, Y = real113$lat)

##QR
rqr10loo <- gwqr(
  formula = rf,
  tau = 0.1,
  cordxy = coords_real,
  distmethod = "euclidean",
  h = 2509,
  kernel = "global",
  adaptive = TRUE,
  data = real113,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

rqr25loo <- gwqr(
  formula = rf,
  tau = 0.25,
  cordxy = coords_real,
  distmethod = "euclidean",
  h = 2509,
  kernel = "global",
  adaptive = TRUE,
  data = real113,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

rqr50loo <- gwqr(
  formula = rf,
  tau = 0.5,
  cordxy = coords_real,
  distmethod = "euclidean",
  h = 2509,
  kernel = "global",
  adaptive = TRUE,
  data = real113,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

rqr75loo <- gwqr(
  formula = rf,
  tau = 0.75,
  cordxy = coords_real,
  distmethod = "euclidean",
  h = 2509,
  kernel = "global",
  adaptive = TRUE,
  data = real113,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

rqr90loo <- gwqr(
  formula = rf,
  tau = 0.9,
  cordxy = coords_real,
  distmethod = "euclidean",
  h = 2509,
  kernel = "global",
  adaptive = TRUE,
  data = real113,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)


##GWQR
rqr10 <- bw_gwqr_cv(
  formula = rf, tau = 0.1,bw_min = 250, bw_max = 2260,
  cordxy = coords_real, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = real113, locallin = FALSE,
  parallel = TRUE, core = 4)

rgwqr10loo <- gwqr(
  formula = rf,
  tau = 0.1,
  cordxy = coords_real,
  distmethod = "euclidean",
  h = rqr10,
  kernel = "bisquare",
  adaptive = TRUE,
  data = real113,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

rqr25 <- bw_gwqr_cv(
  formula = rf, tau = 0.25,bw_min = 250, bw_max = 2260,
  cordxy = coords_real, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = real113, locallin = FALSE,
  parallel = TRUE, core = 4)

rgwqr25loo <- gwqr(
  formula = rf,
  tau = 0.25,
  cordxy = coords_real,
  distmethod = "euclidean",
  h = rqr25,
  kernel = "bisquare",
  adaptive = TRUE,
  data = real113,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

rqr50 <- bw_gwqr_cv(
  formula = rf, tau = 0.5,bw_min = 250, bw_max = 2260,
  cordxy = coords_real, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = real113, locallin = FALSE,
  parallel = TRUE, core = 4)

rgwqr50loo <- gwqr(
  formula = rf,
  tau = 0.5,
  cordxy = coords_real,
  distmethod = "euclidean",
  h = rqr50,
  kernel = "bisquare",
  adaptive = TRUE,
  data = real113,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

rqr75 <- bw_gwqr_cv(
  formula = rf, tau = 0.75,bw_min = 250, bw_max = 2260,
  cordxy = coords_real, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = real113, locallin = FALSE,
  parallel = TRUE, core = 4)

rgwqr75loo <- gwqr(
  formula = rf,
  tau = 0.75,
  cordxy = coords_real,
  distmethod = "euclidean",
  h = rqr75,
  kernel = "bisquare",
  adaptive = TRUE,
  data = real113,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)

rqr90 <- bw_gwqr_cv(
  formula = rf, tau = 0.9,bw_min = 250, bw_max = 2260,
  cordxy = coords_real, distmethod = "euclidean",
  kernel = "bisquare", adaptive = TRUE,
  data = real113, locallin = FALSE,
  parallel = TRUE, core = 4)

rgwqr90loo <- gwqr(
  formula = rf,
  tau = 0.9,
  cordxy = coords_real,
  distmethod = "euclidean",
  h = rqr90,
  kernel = "bisquare",
  adaptive = TRUE,
  data = real113,
  locallin = FALSE,
  parallel = TRUE,
  core = 4
)
# saveRDS(rqr10loo, file="C:/Users/chao6/OneDrive/桌面/yilingsim/real113/qr10realprice.rds")
# saveRDS(rqr25loo, file="C:/Users/chao6/OneDrive/桌面/yilingsim/real113/qr25realprice.rds")
# saveRDS(rqr50loo, file="C:/Users/chao6/OneDrive/桌面/yilingsim/real113/qr50realprice.rds")
# saveRDS(rqr75loo, file="C:/Users/chao6/OneDrive/桌面/yilingsim/real113/qr75realprice.rds")
# saveRDS(rqr90loo, file="C:/Users/chao6/OneDrive/桌面/yilingsim/real113/qr90realprice.rds")
# saveRDS(rgwqr10loo, file="C:/Users/chao6/OneDrive/桌面/yilingsim/real113/gwqr10realprice.rds")
# saveRDS(rgwqr25loo, file="C:/Users/chao6/OneDrive/桌面/yilingsim/real113/gwqr25realprice.rds")
# saveRDS(rgwqr50loo, file="C:/Users/chao6/OneDrive/桌面/yilingsim/real113/gwqr50realprice.rds")
# saveRDS(rgwqr75loo, file="C:/Users/chao6/OneDrive/桌面/yilingsim/real113/gwqr75realprice.rds")
# saveRDS(rgwqr90loo, file="C:/Users/chao6/OneDrive/桌面/yilingsim/real113/gwqr90realprice.rds")
