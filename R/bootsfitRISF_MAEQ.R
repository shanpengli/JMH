##' @export
##'

bootsfitRISF_MAEQ <- function(i, seed, N, increment, beta, tau, gamma1,
                            alpha1, vee1, lambda1, CL, CU, covbw, quadpoint, maxiter, 
                            n.cv, landmark.time, horizon.time, quintile.width,
                            method = "GH") {
  
  
  data <- simJMdataRISF(seed = seed + i, N = N, increment = increment,
                        beta = beta, tau = tau, gamma1 = gamma1,
                        alpha1 = alpha1, vee1 = vee1,
                        lambda1 = lambda1, CL = CL, CU = CU, covbw = covbw)
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
        AllSurv <- list()
        for (j in 1:length(horizon.time)) {
          Surv <- as.data.frame(matrix(0, nrow = nrow(val.cdata), ncol = 2))
          colnames(Surv) <- c("ID", "Surv")
          Surv$ID <- val.cdata$ID
          ## extract estimated survival prob
          for (k in 1:nrow(Surv)) {
            Surv[k, 2] <- survfit$Pred[[k]][j, 2]
          }
          ## group subjects based on survival prob
          quant <- quantile(Surv$Surv, probs = seq(0, 1, by = quintile.width))
          EmpiricalSurv <- rep(NA, groups)
          PredictedSurv <- rep(NA, groups)
          for (i in 1:groups) {
            subquant <- Surv[Surv$Surv > quant[i] &
                               Surv$Surv <= quant[i+1], c(1, 2)]
            quantsubdata <- cdata[cdata$ID %in% subquant$ID, surv.var]
            colnames(quantsubdata) <- c("time", "status")
            fitKM <- survfit(Surv(time, status) ~ 1, data = quantsubdata)
            fitKM.horizon <- try(summary(fitKM, times = horizon.time[j]), silent = TRUE)
            if ('try-error' %in% class(fitKM.horizon)) {
              EmpiricalSurv[i] <- summary(fitKM, times = max(quantsubdata$time))$surv
            } else {
              EmpiricalSurv[i] <- summary(fitKM, times = horizon.time[j])$surv
            }
            PredictedSurv[i] <-mean(subquant$Surv)
          }
          AllSurv[[j]] <- data.frame(EmpiricalSurv, PredictedSurv)
        }
        names(AllSurv) <- horizon.time
        result <- list(AllSurv = AllSurv)
        MAEQ.cv[[t]] <- result
      }
    }
  }
  
  return(MAEQ.cv)
  
}