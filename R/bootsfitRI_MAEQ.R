##' @export
##'

bootsfitRI_MAEQ <- function(i, seed, N, increment, beta, tau, gamma1, gamma2,
                            alpha1, alpha2, vee1, vee2, lambda1, lambda2, CL,
                            CU, covbw, quadpoint, maxiter, 
                            n.cv, landmark.time, horizon.time, quintile.width,
                            method = "GH") {
  

  data <- simJMdataRI(seed = seed + i, N = N, increment = increment,
                      beta = beta, tau = tau, gamma1 = gamma1, gamma2 = gamma2,
                      alpha1 = alpha1, alpha2 = alpha2, vee1 = vee1, vee2 = vee2,
                      lambda1 = lambda1, lambda2 = lambda2,
                      CL = CL, CU = CU, covbw = covbw)
  ydata <- data$ydata
  cdata <- data$cdata
  
  folds <- caret::groupKFold(c(1:nrow(cdata)), k = n.cv)
  groups <- 1/quintile.width
  if (floor(groups) != groups)
    stop("The reciprocal of quintile.width must be an integer.")
  
  MAEQ.cv <- list()

  for (t in 1:n.cv) {
    train.ydata <- ydata[ydata$ID %in% folds[[t]], ]
    train.cdata <- cdata[cdata$ID %in% folds[[t]], ]
    
    fit <- try(JMMLSM(cdata = train.cdata, ydata = train.ydata, 
                      long.formula = Y ~ Z1 + Z2 + Z3 + time,
                      surv.formula = Surv(survtime, cmprsk) ~ var1 + var2 + var3,
                      variance.formula = ~ Z1 + Z2 + Z3 + time, 
                      quadpoint = quadpoint, random = ~ 1|ID), silent = TRUE)
    
    if ('try-error' %in% class(fit)) {
      MAEQ.cv[[t]] <- NULL
    } else if (fit$iter == maxiter) {
      MAEQ.cv[[t]] <- NULL
    } else {
      
      surv.formula <- fit$SurvivalSubmodel
      surv.var <- all.vars(surv.formula)
      
      val.ydata <- ydata[!ydata$ID %in% folds[[t]], ]
      val.cdata <- cdata[!cdata$ID %in% folds[[t]], ]
      
      val.cdata <- val.cdata[val.cdata$survtime > landmark.time, ]
      val.ydata <- val.ydata[val.ydata$ID %in% val.cdata$ID, ]
      val.ydata <- val.ydata[val.ydata$time <= landmark.time/increment, ]
      
      survfit <- try(survfitJMMLSM(fit, seed = seed, ynewdata = val.ydata, cnewdata = val.cdata, 
                                   u = horizon.time, method = method, 
                                   Last.time = rep(landmark.time, nrow(val.cdata)),
                                   obs.time = "time", quadpoint = quadpoint), silent = TRUE)
      
      if ('try-error' %in% class(survfit)) {
        MAEQ.cv[[t]] <- NULL
      } else { 
        
        AllCIF1 <- list()
        AllCIF2 <- list()
        for (j in 1:length(horizon.time)) {
          CIF <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 3))
          colnames(CIF) <- c("ID", "CIF1", "CIF2")
          CIF$ID <- val.cdata$ID
          ## extract estimated CIF
          for (k in 1:nrow(CIF)) {
            CIF[k, 2] <- survfit$Pred[[k]][j, 2]
            CIF[k, 3] <- survfit$Pred[[k]][j, 3]
          }
          ## group subjects based on CIF
          quant1 <- quantile(CIF$CIF1, probs = seq(0, 1, by = quintile.width))
          EmpiricalCIF1 <- rep(NA, groups)
          PredictedCIF1 <- rep(NA, groups)
          
          for (i in 1:groups) {
            subquant <- CIF[CIF$CIF1 > quant1[i] &
                              CIF$CIF1 <= quant1[i+1], 1:2]
            quantsubdata <- val.cdata[val.cdata$ID %in% subquant$ID, surv.var]
            
            quantsubCIF <- GetEmpiricalCIF(data = quantsubdata, 
                                           time = surv.var[1],
                                           status = surv.var[2])
            
            quantsubRisk1 <- quantsubCIF$H1
            ii <- 1
            while (ii <= nrow(quantsubRisk1)) {
              if (quantsubRisk1[ii, 1] > horizon.time[j]) {
                if (ii >= 2) {
                  EmpiricalCIF1[i] <- quantsubRisk1[ii-1, 4]
                } else {
                  EmpiricalCIF1[i] <- 0
                }
                break
              } else {
                ii <- ii + 1
              }
            }
            if (is.na(EmpiricalCIF1[i])) {
              if (nrow(quantsubRisk1) == 0) {
                EmpiricalCIF1[i] <- 0
              } else {
                EmpiricalCIF1[i] <- quantsubRisk1[nrow(quantsubRisk1), 4] 
              }
            }
            PredictedCIF1[i] <- mean(subquant$CIF1)
          }
          AllCIF1[[j]] <- data.frame(EmpiricalCIF1, PredictedCIF1)
          
          quant2 <- quantile(CIF$CIF2, probs = seq(0, 1, by = quintile.width))
          EmpiricalCIF2 <- rep(NA, groups)
          PredictedCIF2 <- rep(NA, groups)
          for (i in 1:(1/quintile.width)) {
            subquant <- CIF[CIF$CIF2 > quant2[i] &
                              CIF$CIF2 <= quant2[i+1], c(1, 3)]
            quantsubdata <- cdata[cdata[, ID] %in% subquant$ID, surv.var]
            
            quantsubCIF <- GetEmpiricalCIF(data = quantsubdata, 
                                           time = surv.var[1],
                                           status = surv.var[2])
            
            quantsubRisk2 <- quantsubCIF$H2
            ii <- 1
            while (ii <= nrow(quantsubRisk2)) {
              if (quantsubRisk2[ii, 1] > horizon.time[j]) {
                if (ii >= 2) {
                  EmpiricalCIF2[i] <- quantsubRisk2[ii-1, 4]
                } else {
                  EmpiricalCIF2[i] <- 0
                }
                break
              } else {
                ii <- ii + 1
              }
            }
            if (is.na(EmpiricalCIF2[i])) {
              if (nrow(quantsubRisk2) == 0) {
                EmpiricalCIF2[i] <- 0
              } else {
                EmpiricalCIF2[i] <- quantsubRisk2[nrow(quantsubRisk2), 4] 
              }
            }
            PredictedCIF2[i] <- mean(subquant$CIF2)
          }
          AllCIF2[[j]] <- data.frame(EmpiricalCIF2, PredictedCIF2)
          
        }
        names(AllCIF1) <- names(AllCIF2) <- horizon.time
        writeLines(paste0("The ", t, " th validation is done!"))
        
        
        
      }
      result <- list(AllCIF1 = AllCIF1, AllCIF2 = AllCIF2)
    }
    MAEQ.cv[[t]] <- result
  }
  
  return(MAEQ.cv)
  
}