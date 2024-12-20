bootsfitRI_AUC2 <- function(i, seed, N, increment, beta, tau, gamma1, gamma2,
                              alpha1, alpha2, vee1, vee2, lambda1, lambda2, CL,
                              CU, covbw, quadpoint, maxiter, 
                              n.cv, landmark.time, horizon.time,
                              method = "GH", silent = TRUE, metric = c("AUC", "Cindex"), model = c("JMH", "FastJM")) {
  
  
  data <- simJMdataRI(seed = seed, N = N, increment = increment,
                      beta = beta, tau = tau, gamma1 = gamma1, gamma2 = gamma2,
                      alpha1 = alpha1, alpha2 = alpha2, vee1 = vee1, vee2 = vee2,
                      lambda1 = lambda1, lambda2 = lambda2,
                      CL = CL, CU = CU, covbw = covbw)
  ydata <- data$ydata
  cdata <- data$cdata
  
  set.seed(seed + i)
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  
  if (!metric %in% c("AUC", "Cindex"))
    stop("Please specify a metric for prediction accuracy assessment: AUC or Cindex.")
  if (!model %in% c("JMH", "FastJM"))
    stop("Please specify a joint model. JMH or FastJM")
  
  AUC.cv <- list()
  
  for (t in 1:n.cv) {
    train.ydata <- ydata[ydata$ID %in% folds[[t]], ]
    train.cdata <- cdata[cdata$ID %in% folds[[t]], ]
    
    if (model == "JMH") {
      fit <- try(JMH::JMMLSM(cdata = train.cdata, ydata = train.ydata, 
                        long.formula = Y ~ Z1 + Z2 + Z3 + time,
                        surv.formula = Surv(survtime, cmprsk) ~ var1 + var2 + var3,
                        variance.formula = ~ Z1 + Z2 + Z3 + time, 
                        quadpoint = quadpoint, random = ~ 1|ID), silent = silent)
    } else {
      fit <- try(FastJM::jmcs(cdata = train.cdata, ydata = train.ydata, 
                              long.formula = Y ~ Z1 + Z2 + Z3 + time,
                              surv.formula = Surv(survtime, cmprsk) ~ var1 + var2 + var3,
                              quadpoint = quadpoint, random = ~ 1|ID, opt = "optim"), silent = silent)
    }
    
    
    if ('try-error' %in% class(fit)) {
      AUC.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      AUC.cv[[t]] <- NULL
    } else {
      
      surv.formula <- fit$SurvivalSubmodel
      surv.var <- all.vars(surv.formula)
      
      val.ydata <- ydata[!ydata$ID %in% folds[[t]], ]
      val.cdata <- cdata[!cdata$ID %in% folds[[t]], ]
      
      val.cdata <- val.cdata[val.cdata$survtime > landmark.time, ]
      val.ydata <- val.ydata[val.ydata$ID %in% val.cdata$ID, ]
      val.ydata <- val.ydata[val.ydata$time <= landmark.time/increment, ]
      
      if (model == "JMH") {
        survfit <- try(JMH::survfitJMMLSM(fit, ynewdata = val.ydata, cnewdata = val.cdata, 
                                     u = horizon.time, method = method, 
                                     Last.time = landmark.time,
                                     obs.time = "time", quadpoint = quadpoint), silent = silent)
      } else {
        survfit <- try(FastJM::survfitjmcs(fit, ynewdata = val.ydata, cnewdata = val.cdata, 
                                           u = horizon.time, method = method, 
                                           Last.time = landmark.time,
                                           obs.time = "time", quadpoint = quadpoint), silent = silent)
      }
      
      if ('try-error' %in% class(survfit)) {
        AUC.cv[[t]] <- NULL
      } else {
        
        CIF <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 3))
        colnames(CIF) <- c("ID", "CIF1", "CIF2")
        CIF$ID <- val.cdata$ID
        CIF$time <- val.cdata[, surv.var[1]]
        CIF$status <- val.cdata[, surv.var[2]]
        mean.AUC <- matrix(NA, nrow = length(horizon.time), ncol = 2)
        for (j in 1:length(horizon.time)) {
          
          ## extract estimated CIF
          for (k in 1:nrow(CIF)) {
            CIF[k, 2] <- survfit$Pred[[k]][j, 2]
            CIF[k, 3] <- survfit$Pred[[k]][j, 3]
          }
          
          if (metric == "AUC") {
            ROC <- timeROC::timeROC(T = CIF$time, delta = CIF$status,
                                    weighting = "marginal",
                                    marker = CIF$CIF1, cause = 1,
                                    times = horizon.time[j])
            
            mean.AUC[j, 1] <- ROC$AUC_1[2]
            
            CIF$status2 <- ifelse(CIF$status == 2, 1,
                                  ifelse(CIF$status == 1, 2, 0)
            )
            
            ROC <- timeROC::timeROC(T = CIF$time, delta = CIF$status2,
                                    weighting = "marginal",
                                    marker = CIF$CIF2, cause = 1,
                                    times = horizon.time[j])
            
            mean.AUC[j, 2] <- ROC$AUC_1[2]
          } else {
            
            mean.AUC[j, 1] <- CindexCR(CIF$time, CIF$status, CIF$CIF1, 1)
            mean.AUC[j, 2] <- CindexCR(CIF$time, CIF$status, CIF$CIF2, 2)
            
          }
          
          
        }
        
        AUC.cv[[t]] <- mean.AUC
        
        writeLines(paste0("The ", t, " th validation is done!")) 
      }
    }
  }
  
  return(AUC.cv)
  
}