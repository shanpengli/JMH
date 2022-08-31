##' @export
##'

bootsfitRI_DP <- function(i, seed, N, increment, beta, tau, gamma1, gamma2,
                          alpha1, alpha2, vee1, vee2, lambda1, lambda2, CL,
                          CU, covbw, quadpoint, maxiter, u, 
                          ynewdata, cnewdata, M) {
  
  data <- simJMdataRI(seed = seed + i, N = N, increment = increment,
                      beta = beta, tau = tau, gamma1 = gamma1, gamma2 = gamma2,
                      alpha1 = alpha1, alpha2 = alpha2, vee1 = vee1, vee2 = vee2,
                      lambda1 = lambda1, lambda2 = lambda2,
                      CL = CL, CU = CU, covbw = covbw)
  ydata <- data$ydata
  cdata <- data$cdata
  
  fit <-   try(JMMLSM(cdata = cdata, ydata = ydata, 
                      long.formula = Y ~ Z1 + Z2 + Z3 + time,
                      surv.formula = Surv(survtime, cmprsk) ~ var1 + var2 + var3,
                      variance.formula = ~ Z1 + Z2 + Z3 + time, maxiter = maxiter, epsilon = 1e-04, 
                      quadpoint = quadpoint, random = ~ 1|ID), silent = TRUE)
  

  if ('try-error' %in% class(fit)) {
    pred1 <- NULL
    pred2 <- NULL
    coef <- list(pred1 = pred1, pred2 = pred2)
    return(coef)
  } else if (fit$iter == maxiter) {
    pred1 <- NULL
    pred2 <- NULL
    coef <- list(pred1 = pred1, pred2 = pred2)
    return(coef)
  } else {
    
    survfit <- try(survfit3JMMLSM(fit, seed = seed + i, ynewdata = ynewdata, cnewdata = cnewdata, 
                                  u = u, M = M, burn.in = 0.2, 
                                  n.nested = 400, simulate = TRUE, quadpoint = quadpoint), silent = TRUE)
    
    if ('try-error' %in% class(survfit)) {
      pred1 <- NULL
      pred2 <- NULL
      coef <- list(pred1 = pred1, pred2 = pred2)
      return(coef)
    } else {
      
      pred1 <- list()
      pred2 <- list()
      for (i in 1:nrow(cnewdata)) {
        pred1[[i]] <- survfit$Pred[[i]]$`Cumulative incidence probabilities for type 1 failure`[, c(2, 4, 5)]
        pred2[[i]] <- survfit$Pred[[i]]$`Cumulative incidence probabilities for type 2 failure`[, c(2, 4, 5)]
        colnames(pred1[[i]]) <- c("mean", "95%Lower", "95%Upper")
        colnames(pred2[[i]]) <- c("mean", "95%Lower", "95%Upper")
        rownames(pred1[[i]]) <- u
        rownames(pred2[[i]]) <- u
      }
      names(pred1) <- cnewdata$ID
      names(pred2) <- cnewdata$ID
      
      coef <- list(pred1 = pred1, pred2 = pred2)
      return(coef)
    }
  }
  
  
  
  
}