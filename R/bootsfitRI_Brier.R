##' @export
##'

bootsfitRI_Brier <- function(i, seed, N, increment, beta, tau, gamma1, gamma2,
                             alpha1, alpha2, vee1, vee2, lambda1, lambda2, CL,
                             CU, covbw, quadpoint, maxiter, 
                             n.cv, landmark.time, horizon.time, method = "GH") {
  
  data <- simJMdataRI(seed = seed + i, N = N, increment = increment,
                      beta = beta, tau = tau, gamma1 = gamma1, gamma2 = gamma2,
                      alpha1 = alpha1, alpha2 = alpha2, vee1 = vee1, vee2 = vee2,
                      lambda1 = lambda1, lambda2 = lambda2,
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
                      quadpoint = 15, random = ~ 1|ID), silent = TRUE)
    
    if ('try-error' %in% class(fit)) {
      Brier.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      Brier.cv[[t]] <- NULL
    } else {
      
      val.ydata <- ydata[!ydata$ID %in% folds[[t]], ]
      val.cdata <- cdata[!cdata$ID %in% folds[[t]], ]
      
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
          } else {
            Gu <- fitKM.horizon$surv
            ## true counting process
            N1 <- vector()
            N2 <- vector()
            Gt <- vector()
            W.IPCW <- vector()
            for (i in 1:nrow(CIF)) {
              if (val.cdata$survtime[i] <= horizon.time[j] && val.cdata$cmprsk[i] == 1) {
                N1[i] <- 1
              } else {
                N1[i] <- 0
              }
              
              if (val.cdata$survtime[i] <= horizon.time[j] && val.cdata$cmprsk[i] == 2) {
                N2[i] <- 1
              } else {
                N2[i] <- 0
              }
              
              if (val.cdata$survtime[i] <= horizon.time[j] && val.cdata$cmprsk[i] != 0) {
                Gt[i] <- summary(fitKM, times = val.cdata$survtime[i])$surv
              } else {
                Gt[i] <- NA
              }
              
              if (val.cdata$survtime[i] > horizon.time[j]) {
                W.IPCW[i] <- 1/(Gu/Gs)
              } else if (val.cdata$survtime[i] <= horizon.time[j] && val.cdata$cmprsk[i] != 0) {
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
            RAWData.Brier$Brier1 <- RAWData.Brier$W.IPCW*(RAWData.Brier$CIF1 - RAWData.Brier$N1)^2
            RAWData.Brier$Brier2 <- RAWData.Brier$W.IPCW*(RAWData.Brier$CIF2 - RAWData.Brier$N2)^2
            
            mean.Brier1 <- sum(RAWData.Brier$Brier1, na.rm = TRUE)/nrow(RAWData.Brier)
            mean.Brier2 <- sum(RAWData.Brier$Brier2, na.rm = TRUE)/nrow(RAWData.Brier)
            mean.Brier[j, 1] <- mean.Brier1
            mean.Brier[j, 2] <- mean.Brier2
          }
          
          
        }
        
        Brier.cv[[t]] <- mean.Brier
      }
      
    }
    
  }
  
  return(Brier.cv)
  
}