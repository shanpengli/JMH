bootsfitMESA_Brier2 <- function(i, seed = 100, N = 200,
                              beta = c(80, -3, -3, 0.5, 5),
                              tau = c(2, 1, 1, 0.05, 0.3),
                              gamma1 = c(0.1, -0.3),
                              gamma2 = c(0.1, -0.3),
                              alpha1 = c(0, -0.03, 0),
                              alpha2 = c(-0.02, 0, 0.02),
                              vee1 = 0.7,
                              vee2 = 0.8,
                              covbw = matrix(c(200, -70, -70, 10, 
                                               -70, 150, 70, -2, 
                                               -70, 70, 150, -2, 
                                               10, -2, -2, 0.7), 
                                             nrow = 4, 
                                             ncol = 4),
                              lambda1 = 1e-5,
                              lambda2 = 3e-5, 
                              lambdaC = 1e-2, 
                              incremin = 1, 
                              incremax = 2.5,
                              Cmax = 18,
                              maxiter = 1000,
                              quadpoint = 6, 
                              n.cv, landmark.time, horizon.time,
                              method = "GH", silent = TRUE, model = c("JMH", "FastJM")) {
  
  
  data <- simJMMESA(seed = seed + i, N = N, beta = beta,
                    tau = tau,
                    gamma1 = gamma1,
                    gamma2 = gamma2,
                    alpha1 = alpha1,
                    alpha2 = alpha2,
                    vee1 = vee1,
                    vee2 = vee2,
                    lambda1 = lambda1,
                    lambda2 = lambda2,
                    lambdaC = lambdaC,
                    incremin = incremin,
                    incremax = incremax,
                    Cmax = Cmax,
                    covbw = covbw)
  ydata <- data$ydata
  cdata <- data$cdata
  rate <- data$rate
  
  set.seed(seed + i)
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  if (!model %in% c("JMH", "FastJM"))
    stop("Please specify a joint model. JMH or FastJM")
  
  Brier.cv <- list()
  
  for (t in 1:n.cv) {
    train.ydata <- ydata[ydata$ID %in% folds[[t]], ]
    train.cdata <- cdata[cdata$ID %in% folds[[t]], ]
    
    if (model == "JMH") {
      fit <- try(JMH::JMMLSM(cdata = train.cdata, ydata = train.ydata, 
                             long.formula = Y ~ timens1 + timens2 + Z1 + Z2,
                             surv.formula = Surv(survtime, cmprsk) ~ var1 + var2,
                             variance.formula = ~ timens1 + timens2 + Z1 + Z2, 
                             quadpoint = quadpoint, random = ~ timens1 + timens2|ID,
                             maxiter = maxiter,
                             epsilon = 5e-3,
                             initial.para = list(beta = beta,
                                                 tau = tau,
                                                 gamma1 = gamma1,
                                                 gamma2 = gamma2,
                                                 alpha1 = alpha1,
                                                 alpha2 = alpha2,
                                                 vee1 = vee1,
                                                 vee2 = vee2,
                                                 Sig = covbw)), silent = silent)
    } else {
      fit <- try(FastJM::jmcs(cdata = train.cdata, ydata = train.ydata, 
                              long.formula = Y ~ timens1 + timens2 + Z1 + Z2,
                              surv.formula = Surv(survtime, cmprsk) ~ var1 + var2,
                              quadpoint = quadpoint, random = ~ timens1 + timens2|ID,
                              maxiter = maxiter,
                              tol = 5e-3), silent = silent)
    }
    
    
    if ('try-error' %in% class(fit)) {
      Brier.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      Brier.cv[[t]] <- NULL
    } else {
      
      surv.formula <- fit$SurvivalSubmodel
      surv.var <- all.vars(surv.formula)
      
      ## fit a Kalplan-Meier estimator
      
      val.ydata <- ydata[!ydata$ID %in% folds[[t]], ]
      val.cdata <- cdata[!cdata$ID %in% folds[[t]], ]
      
      ## fit a Kalplan-Meier estimator
      fitKM <- survival::survfit(Surv(survtime, cmprsk == 0) ~ 1, data = val.cdata)
      
      val.cdata <- val.cdata[val.cdata$survtime > landmark.time, ]
      val.ydata <- val.ydata[val.ydata$ID %in% val.cdata$ID, ]
      val.ydata <- val.ydata[val.ydata$time <= landmark.time, ]
      
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
        Brier.cv[[t]] <- NULL
      } else {
        
        CIF <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 3))
        colnames(CIF) <- c("ID", "CIF1", "CIF2")
        CIF$ID <- val.cdata$ID
        Gs <- summary(fitKM, times = landmark.time)$surv
        mean.Brier <- matrix(NA, nrow = length(horizon.time), ncol = 2)
        for (j in 1:length(horizon.time)) {
          fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
          if ('try-error' %in% class(fitKM.horizon)) {
            mean.Brier[j, 1] <- NA
            mean.Brier[j, 2] <- NA
            mean.MAE[j, 1] <- NA
            mean.MAE[j, 2] <- NA
          } else {
            Gu <- fitKM.horizon$surv
            ## true counting process
            N1 <- vector()
            N2 <- vector()
            Gt <- vector()
            W.IPCW <- vector()
            for (i in 1:nrow(CIF)) {
              if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] == 1) {
                N1[i] <- 1
              } else {
                N1[i] <- 0
              }
              
              if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] == 2) {
                N2[i] <- 1
              } else {
                N2[i] <- 0
              }
              
              if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] != 0) {
                Gt[i] <- summary(fitKM, times = val.cdata[i, surv.var[1]])$surv
              } else {
                Gt[i] <- NA
              }
              
              if (val.cdata[i, surv.var[1]] > horizon.time[j]) {
                W.IPCW[i] <- 1/(Gu/Gs)
              } else if (val.cdata[i, surv.var[1]] <= horizon.time[j] && val.cdata[i, surv.var[2]] != 0) {
                W.IPCW[i] <- 1/(Gt[i]/Gs)
              } else {
                W.IPCW[i] <- NA
              }
            }
            ## extract estimated CIF
            for (k in 1:nrow(CIF)) {
              CIF[k, 2] <- survfit$Pred[[k]][j, 2]
              CIF[k, 3] <- survfit$Pred[[k]][j, 3]
            }
            
            RAWData.Brier <- data.frame(CIF, N1, N2, W.IPCW)
            colnames(RAWData.Brier)[1:3] <- c("ID", "CIF1", "CIF2")
            RAWData.Brier$Brier1 <- RAWData.Brier$W.IPCW*
              abs(RAWData.Brier$CIF1 - RAWData.Brier$N1)^2
            RAWData.Brier$Brier2 <- RAWData.Brier$W.IPCW*
              abs(RAWData.Brier$CIF2 - RAWData.Brier$N2)^2
            mean.Brier1 <- sum(RAWData.Brier$Brier1, na.rm = TRUE)/nrow(RAWData.Brier)
            mean.Brier2 <- sum(RAWData.Brier$Brier2, na.rm = TRUE)/nrow(RAWData.Brier)
            mean.Brier[j, 1] <- mean.Brier1
            mean.Brier[j, 2] <- mean.Brier2
          }
          
          
        }
        
        Brier.cv[[t]] <- mean.Brier
        
        writeLines(paste0("The ", t, " th validation is done!"))
        
      }
    }
  }
  
  return(Brier.cv)
  
}