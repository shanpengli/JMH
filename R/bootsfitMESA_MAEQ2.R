bootsfitMESA_MAEQ2 <- function(i, seed = 100, N = 200,
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
                               n.cv, landmark.time, horizon.time, quantile.width = 0.25,
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
  groups <- 1/quantile.width
  if (floor(groups) != groups)
    stop("The reciprocal of quantile.width must be an integer.")
  if (!model %in% c("JMH", "FastJM"))
    stop("Please specify a joint model. JMH or FastJM")
  
  MAEQ.cv <- list()
  
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
          quant1 <- quantile(CIF$CIF1, probs = seq(0, 1, by = quantile.width))
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
          
          quant2 <- quantile(CIF$CIF2, probs = seq(0, 1, by = quantile.width))
          EmpiricalCIF2 <- rep(NA, groups)
          PredictedCIF2 <- rep(NA, groups)
          for (i in 1:(1/quantile.width)) {
            subquant <- CIF[CIF$CIF2 > quant2[i] &
                              CIF$CIF2 <= quant2[i+1], c(1, 3)]
            quantsubdata <- cdata[cdata$ID %in% subquant$ID, surv.var]
            
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
        result <- list(AllCIF1 = AllCIF1, AllCIF2 = AllCIF2)
        MAEQ.cv[[t]] <- result
        writeLines(paste0("The ", t, " th validation is done!"))
      }
    }
  }
  
  return(MAEQ.cv)
  
}