#### GWQRF Estimation
gwqrf_par <- function(formula, cordx, cordy, h, kernel, indata, distkm=F, corenum, 
                      model_type = "quantile",tau=0.5, mtry=NULL, num_trees=500, max_depth=0, min_node_size=5) {  
  time1 = Sys.time()
  cat("執行估計:", format(time1, "%H:%M:%S"), "\n")
  mat <- model.frame(formula, data = indata)
  y = mat[, 1]
  x = as.matrix(model.matrix(formula, data = indata))  
  ox = x[,-1]  
  nvar = ncol(ox)
  nobs = length(y) 
  if (is.null(mtry)) {mtry = max(floor(nvar/3), 2)}
  
  cl <- makeCluster(corenum)
  clusterExport(cl, c("nobs","nvar","x","ox","formula","indata","mtry","formula",
                      "y","cordx","cordy","h","kernel","distkm", "tau","model_type",
                      "num_trees", "max_depth", "min_node_size"), envir=environment())    
  clusterEvalQ(cl,c(library(ranger),library(geosphere)))#,library(fastshap))) 
  sp = clusterSplit(cl, 1:nobs)
  
  # 估計函數
  estimate <- function(times){
    yhat_rf_oob=c(); yhat_rf_full=c(); #shapv=c()
    #forest_full=list()
    muhat <- numeric(length(times))
    for(j in seq_along(times)){
      i <- times[j]
      if (distkm){
        d <- geosphere::distCosine(cbind(cordx, cordy), cbind(cordx[i], cordy[i]), r=6371)} 
      else { d <- sqrt((cordx - cordx[i])^2 + (cordy - cordy[i])^2) }
      rh <- min(max(round(h), 2L), nobs)  # rh 限制在合法範圍內
      disth <- sort(d)[rh]
      
      # kernal
      if (kernel == "global") {
        w <- rep(1, nobs) }
      else if (kernel == "gaussian") {
        w <- exp(-0.5 * (d / disth)^2)}
      else if (kernel == "exponential") {
        w <- exp(-d / disth) }
      else { 
        w1 = (1-(d/disth)**2)**2    #adaptive bisquare
        w = rep(0,length(w1))
        index=which(d<=disth)
        w[index] = w1[index]
      }
      
      #LOO
      w_oob <- w
      w_oob[i] <- 0
      valid_oob <- which(w_oob > 0) 
      #regdata_oob <- data.frame(indata[valid_oob,,drop=FALSE],w = w_oob[valid_oob]) 
      #regdata_oob <- data.frame(indata[,,drop=FALSE],w = w_oob)         
      rstoob = ranger::ranger(formula = formula, data = indata[,,drop=FALSE], case.weights = w_oob, 
                              quantreg = (model_type == "quantile"),
                              num.trees = num_trees, mtry = mtry, num.threads = 1, 
                              max.depth = max_depth, min.node.size = min_node_size, #keep.inbag = TRUE,
                              #importance = "impurity",
                              seed = 123 
      )
      newdata = indata[i, , drop = FALSE]        
      if(model_type == "quantile") {
        yhat_rf_oob[j] = predict(rstoob, data = newdata, type = "quantiles", quantiles = tau)$predictions[1,1] }
      else {yhat_rf_oob[j] = predict(rstoob, data = newdata)$predictions[1,1] }
      
      #非LOO
      #w_full  <- w
      #valid_full <- which(w_full > 0)
      #regdata_full <- data.frame(indata[valid_full,,drop=FALSE],w = w_full[valid_full])  
      #regdata_full <- data.frame(indata[,,drop=FALSE],w = w)       
      rstfull = ranger::ranger(formula = formula, data =indata[,,drop=FALSE], case.weights = w, 
                               quantreg = (model_type == "quantile"),
                               num.trees = num_trees, mtry = mtry, num.threads = 1, 
                               max.depth = max_depth, min.node.size = min_node_size, #keep.inbag = TRUE,
                               #importance = "impurity",
                               seed = 123 
      )
      newdata = indata[i, , drop = FALSE]        
      if(model_type == "quantile") {
        yhat_rf_full[j] = predict(rstfull, data = newdata, type = "quantiles", quantiles = tau)$predictions[1,1]  }
      else {yhat_rf_full[j] = predict(rstfull, data = newdata)$predictions[1,1] }
      #forest_full[[i]] = rstfull
      
      # #SHAP (LOO)
      # shap_X = regdata_full[, all.vars(formula)[-1], drop = FALSE]  #Shap values
      # shap_test = newdata[, all.vars(formula)[-1], drop = FALSE]
      # set.seed(123)
      # if(model_type == "quantile") {
      #   shapv1 = fastshap::explain(rstfull, X = shap_X, 
      #                              pred_wrapper = function(object, newdata){ as.numeric(predict(object,data = newdata,type = "quantiles",quantiles = tau)$predictions)},
      #                              newdata = shap_test,nsim=50)
      # }
      # else { shapv1 = fastshap::explain(rstfull, X = shap_X, 
      #                                   pred_wrapper = function(object, newdata){ as.numeric(predict(object,data = newdata)$predictions)},
      #                                   newdata = shap_test,nsim=50)
      # }
      # shapv = rbind(shapv, shapv1)
    }
    output = list(yhat_rf_oob=yhat_rf_oob,yhat_rf_full=yhat_rf_full)#,shapv=shapv)
    return(output)
  }
  
  res <- clusterApply(cl, sp, fun = estimate)
  stopCluster(cl)
  yhat_rf_oob  <- unlist(lapply(res, function(x) x$yhat_rf_oob))
  yhat_rf_full <- unlist(lapply(res, function(x) x$yhat_rf_full))
  #shapv <- do.call(rbind, lapply(res, function(x) x$shapv))
  
  error_rf_oob=y-yhat_rf_oob
  error_rf_full=y-yhat_rf_full
  if (model_type == "quantile") {
    wsum_rf_oob=sum((error_rf_oob*(tau-ifelse(error_rf_oob<0,1,0))))
    wsum_rf_full=sum((error_rf_full*(tau-ifelse(error_rf_full<0,1,0)))) }
  else { 
    wsum_rf_oob=sum(error_rf_oob^2)
    wsum_rf_full=sum(error_rf_full^2)
  } 
  
  time2 <- Sys.time()
  run_time <- time2 - time1
  cat("運算完成。總耗時:", format(run_time), "\n")
  out=list(y=y,ox=ox,yhat_rf_oob=yhat_rf_oob, yhat_rf_full=yhat_rf_full,#shaprf=shapv,
           error_rf_oob=error_rf_oob, error_rf_full=error_rf_full,
           wsum_rf_oob=wsum_rf_oob,wsum_rf_full=wsum_rf_full)  
  return(out)
}

tryqrf50 = gwqrf_par(formula = rf, cordx = realcoords$lon, cordy = realcoords$lat, h=2509, 
                     kernel = "global", indata = real113, model_type="quantile",tau = 0.5, distkm = FALSE,
                     corenum = 4, mtry=15, num_trees=159, max_depth=27, min_node_size=2)


###########
###GWQRF###
###########
gwqrf05 <- readRDS("C:/Users/vivia/Downloads/gwqrf50.rds")


####
dat <- read.csv("C:/Users/vivia/Desktop/Paper/real113/real113q05.csv")
tau <- 0.5

pinball_point <- function(y, qhat, tau) {
  u <- y - qhat
  pmax(tau * u, (tau - 1) * u)
}

error_global_full <- dat$y - dat$yhat_global_full
error_global_oob  <- dat$y - dat$yhat_global_loo

pinball_global_full <- pinball_point(dat$y, dat$yhat_global_full, tau)
pinball_global_oob  <- pinball_point(dat$y, dat$yhat_global_loo, tau)

# 計算總 pinball loss
wsum_global_full <- sum(pinball_global_full, na.rm = TRUE)
wsum_global_oob  <- sum(pinball_global_oob, na.rm = TRUE)

# 計算平均 pinball loss
mean_global_full <- mean(pinball_global_full, na.rm = TRUE)
mean_global_oob  <- mean(pinball_global_oob, na.rm = TRUE)

# 輸出結果
cat("=== Global Quantile Model Pinball Loss ===\n")
cat("Tau =", tau, "\n\n")

cat("Full model:\n")
cat("  Total pinball loss =", wsum_global_full, "\n")
cat("  Mean  pinball loss =", mean_global_full, "\n\n")

cat("OOB model:\n")
cat("  Total pinball loss =", wsum_global_oob, "\n")
cat("  Mean  pinball loss =", mean_global_oob, "\n")

