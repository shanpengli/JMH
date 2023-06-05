bootsfitRISF_Brier <- function(i, seed, N, increment, beta, tau, gamma1,
                             alpha1, vee1, lambda1, CL,
                             CU, covbw, quadpoint, maxiter, 
                             n.cv, landmark.time, horizon.time, method = "GH") {
  
  data <- simJMdataRISF(seed = seed + i, N = N, increment = increment,
                      beta = beta, tau = tau, gamma1 = gamma1, 
                      alpha1 = alpha1, vee1 = vee1, 
                      lambda1 = lambda1, 
                      CL = CL, CU = CU, covbw = covbw)
  ydata <- data$ydata
  cdata <- data$cdata
  
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  Brier.cv <- list()
  
  for (t in 1:n.cv) {
    train.ydata <- ydata[ydata$ID %in% folds[[t]], ]
    train.cdata <- cdata[cdata$ID %in% folds[[t]], ]
    
    fit <- try(JMMLSM(cdata = train.cdata, ydata = train.ydata, 
                      long.formula = Y ~ Z1 + Z2 + Z3 + time,
                      surv.formula = Surv(survtime, cmprsk) ~ var1 + var2 + var3,
                      variance.formula = ~ Z1 + Z2 + Z3 + time, 
                      quadpoint = quadpoint, random = ~ 1|ID), silent = TRUE)
    
    if ('try-error' %in% class(fit)) {
      Brier.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      Brier.cv[[t]] <- NULL
    } else {
      
      val.ydata <- ydata[!ydata$ID %in% folds[[t]], ]
      val.cdata <- cdata[!cdata$ID %in% folds[[t]], ]
      
      surv.formula <- fit$SurvivalSubmodel
      surv.var <- all.vars(surv.formula)
      ## fit a Kalplan-Meier estimator
      fitKM <- survfit(Surv(survtime, cmprsk == 0) ~ 1, data = val.cdata)
      
      val.cdata <- val.cdata[val.cdata$survtime > landmark.time, ]
      val.ydata <- val.ydata[val.ydata$ID %in% val.cdata$ID, ]
      val.ydata <- val.ydata[val.ydata$time <= landmark.time/increment, ]
      
      survfit <- try(survfitJMMLSM(fit, seed = seed, ynewdata = val.ydata, cnewdata = val.cdata, 
                                   u = horizon.time, method = method, 
                                   Last.time = rep(landmark.time, nrow(val.cdata)),
                                   obs.time = "time", quadpoint = quadpoint), silent = TRUE)
      
      if ('try-error' %in% class(survfit)) {
        Brier.cv[[t]] <- NULL
      } else {
        Surv <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 2))
        colnames(Surv) <- c("ID", "Surv")
        Surv$ID <- val.cdata$ID
        Gs <- summary(fitKM, times = landmark.time)$surv
        mean.Brier <- matrix(NA, nrow = length(horizon.time), ncol = 1)
        mean.MAE <- matrix(NA, nrow = length(horizon.time), ncol = 1)
        for (j in 1:length(horizon.time)) {
          fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
          if ('try-error' %in% class(fitKM.horizon)) {
            mean.Brier[j, 1] <- NA
          } else {
            Gu <- fitKM.horizon$surv
            ## true counting process
            N1 <- vector()
            Gt <- vector()
            W.IPCW <- vector()
            for (i in 1:nrow(Surv)) {
              if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] == 1) {
                N1[i] <- 1
                Gt[i] <- summary(fitKM, times = val.cdata[i, surv.var[1]])$surv
              } else {
                N1[i] <- 0
                Gt[i] <- NA
              }
              
              if (val.cdata[i, surv.var[1]] > horizon.time[j]) {
                W.IPCW[i] <- 1/(Gu/Gs)
              } else if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] == 1) {
                W.IPCW[i] <- 1/(Gt[i]/Gs)
              } else {
                W.IPCW[i] <- NA
              }
            }
            ## extract estimated Survival probability
            for (k in 1:nrow(Surv)) {
              Surv[k, 2] <- survfit$Pred[[k]][j, 2]
            }
            
            RAWData.Brier <- data.frame(Surv, N1, W.IPCW)
            colnames(RAWData.Brier)[1:2] <- c("ID", "Surv")
            RAWData.Brier$Brier <- RAWData.Brier$W.IPCW*
              abs(1 - RAWData.Brier$Surv - RAWData.Brier$N1)^2
            mean.Brier[j, 1] <- sum(RAWData.Brier$Brier, na.rm = TRUE)/nrow(RAWData.Brier)
          }
        }
        Brier.cv[[t]] <- mean.Brier
      }
      
    }
    
  }
  
  return(Brier.cv)
  
}