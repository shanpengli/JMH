##' Joint modeling of longitudinal continuous data and competing risks
##' @title Joint Modelling for Continuous outcomes
##' @param ydata a longitudinal data frame in long format.
##' @param cdata a survival data frame with competing risks or single failure.
##' Each subject has one data entry.
##' @param long.formula a formula object with the response variable and fixed effects covariates
##' to be included in the longitudinal sub-model.
##' @param surv.formula a formula object with the survival time, event indicator, and the covariates
##' to be included in the survival sub-model.
##' @param variance.formula an one-sided formula object with the fixed effects covariates to model the variance of longituidnal sub-model.
##' @param random a one-sided formula object describing the random effects part of the longitudinal sub-model.
##' For example, fitting a random intercept model takes the form ~ 1|ID.
##' Alternatively. Fitting a random intercept and slope model takes the form ~ x1 + ... + xn|ID.
##' @param maxiter the maximum number of iterations of the EM algorithm that the function will perform. Default is 10000.
##' @param epsilon Tolerance parameter. Default is 0.0001.
##' @param quadpoint the number of pseudo-adaptive Gauss-Hermite quadrature points
##' to be chosen for numerical integration. Default is 6 which produces stable estimates in most dataframes.
##' @param print.para Print detailed information of each iteration. Default is FALSE, i.e., not to print the iteration details.
##' @export
##' 

JMMLSM <- function(cdata, ydata,
                long.formula,
                surv.formula,
                variance.formula, 
                random,
                maxiter = 1000, epsilon = 1e-04, 
                quadpoint = 10, print.para = FALSE) {
  
  random.form <- all.vars(random)
  ID <- random.form[length(random.form)]
  if (length(random.form) == 1) {
    RE <- NULL
    model <- "intercept"
  } else {
    RE <- random.form[-length(random.form)]
    model <- "interslope"
  }
  
  # getinit <- Getinit(cdata = cdata, ydata = ydata, long.formula = long.formula,
  #                    surv.formula = surv.formula, variance.formula = variance.formula,
  #                    model = model, ID = ID, RE = RE)

  getinit <- GetinitFake(cdata, ydata, long.formula, surv.formula, variance.formula,
                         model, ID, RE)

  cdata <- getinit$cdata
  
  survival <- all.vars(surv.formula)
  status <- as.vector(cdata[, survival[2]])
  if (prod(c(0, 1, 2) %in% unique(status))) {
    ## initialize parameters
    
    getriskset <- Getriskset(cdata = cdata, surv.formula = surv.formula)
    
    ## number of distinct survival time
    H01 <- getriskset$tablerisk1
    H02 <- getriskset$tablerisk2
    
    ## initialize parameters
    beta <- getinit$beta
    tau <- getinit$tau
    gamma1 <- getinit$gamma1
    gamma2 <- getinit$gamma2
    alpha1 <- getinit$alpha1
    alpha2 <- getinit$alpha2
    vee1 <- getinit$vee1
    vee2 <- getinit$vee2
    Sig <- getinit$Sig
    p1a <- ncol(Sig) - 1
    if (p1a == 2) Sigb <- Sig[1:2, 1:2]
    if (p1a == 1) Sigb <- as.matrix(Sig[1, 1])
    
    CompetingRisk <- TRUE
  } else {
    ## initialize parameters
    
    getriskset <- Getriskset(cdata = cdata, surv.formula = surv.formula)
    
    ## number of distinct survival time
    H01 <- as.matrix(getriskset$tablerisk1)
    
    ## initialize parameters
    beta <- getinit$beta
    tau <- getinit$tau
    gamma1 <- getinit$gamma1
    alpha1 <- getinit$alpha1
    vee1 <- getinit$vee1
    Sig <- getinit$Sig
    p1a <- ncol(Sig) - 1
    if (p1a == 2) Sigb <- Sig[1:2, 1:2]
    if (p1a == 1) Sigb <- as.matrix(Sig[1, 1])
    CompetingRisk <- FALSE
  }
  
  # beta <- getinit$Beta$betahat
  # tau <- getinit$Tau$betahat
  # gamma1 <- getinit$gamma1
  # gamma2 <- getinit$gamma2
  # alpha1 <- as.vector(rep(0, ncol(getinit$Z)))
  # alpha2 <- alpha1
  # vee1 <- 0
  # vee2 <- vee1
  # 
  # Sigb <- matrix(c(1200, -21, -21, 0.38), nrow = 2, ncol = 2)
  # Sigw <- 0.1
  # Sig <- matrix(0, nrow = (nrow(Sigb)+1), ncol = (nrow(Sigb)+1))
  # for (i in 1:(nrow(Sigb))) {
  #   for (j in 1:(nrow(Sigb))) {
  #     Sig[i, j] <- Sigb[i, j]
  #   }
  # }
  # Sig[(nrow(Sigb)+1), (nrow(Sigb)+1)] <- Sigw
  
  gq_vals <- statmod::gauss.quad(n = quadpoint, kind = "hermite")
  xs <- gq_vals$nodes
  ws <- gq_vals$weights
  
  if (nrow(Sigb) == 2) {
    xsmatrix <- matrix(0, nrow = 3, ncol = quadpoint^3)
    wsmatrix <- xsmatrix
    xsmatrix[3, ] <- rep(xs, quadpoint^2)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(xs[i], quadpoint)
      Total <- c(Total, sub)
    }
    xsmatrix[2, ] <- rep(Total, quadpoint)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(xs[i], quadpoint^2)
      Total <- c(Total, sub)
    }
    xsmatrix[1, ] <- Total
    xsmatrix <- t(xsmatrix)
    
    wsmatrix[3, ] <- rep(ws, quadpoint^2)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(ws[i], quadpoint)
      Total <- c(Total, sub)
    }
    wsmatrix[2, ] <- rep(Total, quadpoint)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(ws[i], quadpoint^2)
      Total <- c(Total, sub)
    }
    wsmatrix[1, ] <- Total
    wsmatrix <- t(wsmatrix)
  } else {
    xsmatrix <- matrix(0, nrow = 2, ncol = quadpoint^2)
    wsmatrix <- xsmatrix
    xsmatrix[2, ] <- rep(xs, quadpoint)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(xs[i], quadpoint)
      Total <- c(Total, sub)
    }
    xsmatrix[1, ] <- Total
    xsmatrix <- t(xsmatrix)
    
    wsmatrix[2, ] <- rep(ws, quadpoint)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(ws[i], quadpoint)
      Total <- c(Total, sub)
    }
    wsmatrix[1, ] <- Total
    wsmatrix <- t(wsmatrix)
  }
  
  iter=0
  
  Z <- getinit$Z
  X1 <- getinit$X1
  W <- getinit$W
  Y <- getinit$Y
  X2 <- getinit$X2
  survtime <- getinit$survtime
  cmprsk <- getinit$cmprsk
  mdata <- getinit$mdata
  n <- nrow(mdata)
  mdata <- as.data.frame(mdata)
  mdata <- as.vector(mdata$ni)
  mdataS <- rep(0, n) 
  mdataS[1] <- 1
  mdataCum <- cumsum(mdata)
  mdata2 <- mdata - 1
  mdataS[2:n] <- mdataCum[2:n] - mdata2[2:n]
  
  if (CompetingRisk == TRUE) {
    repeat
    {
      iter <- iter + 1
      prebeta <- beta
      pretau <- tau
      pregamma1 <- gamma1
      pregamma2 <- gamma2
      prealpha1 <- alpha1
      prealpha2 <- alpha2
      prevee1 <- vee1
      prevee2 <- vee2
      preSig <- Sig
      preH01 <- H01
      preH02 <- H02
      if (print.para) {
        writeLines("iter is:")
        print(iter)
        writeLines("beta is:")
        print(beta)
        writeLines("tau is:")
        print(tau)
        writeLines("gamma1 is:")
        print(gamma1)
        writeLines("gamma2 is:")
        print(gamma2)
        writeLines("alpha1 is:")
        print(alpha1)
        writeLines("alpha2 is:")
        print(alpha2)
        writeLines("vee1 is:")
        print(vee1)
        writeLines("vee2 is:")
        print(vee2)
        writeLines("Sig is:")
        print(Sig)
      }
      
      GetEfun <- GetE(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, 
                      Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix)
      
      GetMpara <- GetM(GetEfun, beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, 
                       Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS)
      
      beta <- GetMpara$beta
      tau <- GetMpara$tau
      gamma1 <- GetMpara$gamma1
      gamma2 <- GetMpara$gamma2
      alpha1 <- GetMpara$alpha1
      alpha2 <- GetMpara$alpha2
      vee1 <- GetMpara$vee1
      vee2 <- GetMpara$vee2
      Sig <- GetMpara$Sig
      H01 <- GetMpara$H01
      H02 <- GetMpara$H02
      
      if((Diff(beta, prebeta, tau, pretau, gamma1, pregamma1, gamma2, pregamma2,
               alpha1, prealpha1, alpha2, prealpha2, vee1, prevee1, vee2, prevee2,
               Sig, preSig, H01, preH01, H02, preH02, epsilon) == 0) || (iter == maxiter) || (!is.list(GetEfun))
         || (!is.list(GetMpara))) {
        break
      }
      
    }
    
    if (iter == maxiter) {
      writeLines("program stops because of nonconvergence")
      convergence = 0
      result <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, 
                     H02, Sig, iter, convergence)
      
      names(result) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1",
                         "vee2", "H01", "H02", "Sig", "iter", "convergence")
      
      return(result)
    } else if (!is.list(GetEfun)) {
      writeLines("Something wrong in the E steps")
      beta <- NULL
      tau <- NULL
      gamma1 <- NULL
      gamma2 <- NULL
      alpha1 <- NULL
      alpha2 <- NULL
      vee1 <- NULL
      vee2 <- NULL
      H01 <- NULL
      H02 <- NULL
      Sig <- NULL
      iter <- NULL
      result <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, 
                     H02, Sig, iter)
      names(result) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1",
                         "vee2", "H01", "H02", "Sig", "iter")
      
      return(result)
    } else if (!is.list(GetMpara)) {
      writeLines("Something wrong in the M steps")
      beta <- NULL
      tau <- NULL
      gamma1 <- NULL
      gamma2 <- NULL
      alpha1 <- NULL
      alpha2 <- NULL
      vee1 <- NULL
      vee2 <- NULL
      H01 <- NULL
      H02 <- NULL
      Sig <- NULL
      iter <- NULL
      result <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, 
                     H02, Sig, iter)
      names(result) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1",
                         "vee2", "H01", "H02", "Sig", "iter")
      
      return(result)
    } else {
      
      convergence = 1
      
      GetEfun <- GetE(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, 
                      Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix)
      
      FUNENW <- as.vector(GetEfun$FUNENW)
      FUNBENW <- as.matrix(GetEfun$FUNBENW)
      FUNBS <- as.matrix(GetEfun$FUNBS)
      FUNBW <- as.matrix(GetEfun$FUNBW)
      FUNWS <- as.vector(GetEfun$FUNWS)
      FUNBSENW <- as.matrix(GetEfun$FUNBSENW) 
      FUNEC <- as.matrix(GetEfun$FUNEC)
      FUNBEC <- as.matrix(GetEfun$FUNBEC)
      FUNBSEC <- as.matrix(GetEfun$FUNBSEC)
      FUNWEC <- as.matrix(GetEfun$FUNWEC)
      FUNWSEC <- as.matrix(GetEfun$FUNWSEC)
      FUNB <- as.matrix(GetEfun$FUNB)
      FUNW <- as.vector(GetEfun$FUNW)
      
      getcov <- getCov(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2,  H01, 
                       H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS,
                       FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC,
                       FUNBSEC, FUNWEC, FUNWSEC,FUNB, FUNW)
      
      vcov <- getcov$vcov
      sebeta <- getcov$sebeta
      setau <- getcov$setau
      segamma1 <- getcov$segamma1
      segamma2 <- getcov$segamma2
      sealpha1 <- getcov$sealpha1
      sealpha2 <- getcov$sealpha2
      sevee1 <- getcov$sevee1
      sevee2 <- getcov$sevee2
      seSig <- getcov$seSig
      
      long <- all.vars(long.formula)
      survival <- all.vars(surv.formula)
      
      names(beta) <- c("intercept", long[-1])
      variance.form <- all.vars(variance.formula)
      names(tau) <- paste0(c("intercept", variance.form), "_var")
      
      names(gamma1) <- paste0(survival[-(1:2)], "_1")
      names(gamma2) <- paste0(survival[-(1:2)], "_2")
      
      PropComp <- as.data.frame(table(cdata[, survival[2]]))
      
      LongOut <- long[1]
      LongX <- paste0(long[-1], collapse = "+")
      FunCall_long <- as.formula(paste(LongOut, LongX, sep = "~"))
      
      LongVarOut <- "log(sigma^2)"
      LongVarX <- paste0(variance.form, collapse = "+") 
      FunCall_longVar <- as.formula(paste(LongVarOut, LongVarX, sep = "~"))
      
      SurvOut <- paste0("Surv(", survival[1], ",", survival[2], ")")
      SurvX <- paste0(survival[-(1:2)], collapse = "+")
      FunCall_survival <- as.formula(paste(SurvOut, SurvX, sep = "~"))
      
      ## return the joint modelling result
      mycall <- match.call()
      
      result <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, 
                     H02, Sig, iter, convergence, vcov, sebeta, setau, segamma1,
                     segamma2, sealpha1, sealpha2, sevee1, sevee2, seSig, CompetingRisk,
                     quadpoint, ydata, cdata, PropComp, FunCall_long, FunCall_longVar,
                     FunCall_survival, random, mycall)
      
      names(result) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1",
                         "vee2", "H01", "H02", "Sig", "iter", "convergence", "vcov",
                         "sebeta", "setau", "segamma1", "segamma2", "sealpha1", "sealpha2", 
                         "sevee1", "sevee2", "seSig", "CompetingRisk", "quadpoint",
                         "ydata", "cdata", "PropEventType", "LongitudinalSubmodelmean",
                         "LongitudinalSubmodelvariance", "Survivalsubmodel", "random",
                         "call")
      
      class(result) <- "JMH"
      
      return(result)
    }
  } else {
    repeat
    {
      iter <- iter + 1
      prebeta <- beta
      pretau <- tau
      pregamma1 <- gamma1
      prealpha1 <- alpha1
      prevee1 <- vee1
      preSig <- Sig
      preH01 <- H01
      if (print.para) {
        writeLines("iter is:")
        print(iter)
        writeLines("beta is:")
        print(beta)
        writeLines("tau is:")
        print(tau)
        writeLines("gamma1 is:")
        print(gamma1)
        writeLines("alpha1 is:")
        print(alpha1)
        writeLines("vee1 is:")
        print(vee1)
        writeLines("Sig is:")
        print(Sig)
      }
      
      GetEfun <- GetESF(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, 
                        X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix)
      
      
      GetMpara <- GetMSF(GetEfun, beta, tau, gamma1, alpha1, vee1, H01,
                       Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS)
      
      beta <- GetMpara$beta
      tau <- GetMpara$tau
      gamma1 <- GetMpara$gamma1
      alpha1 <- GetMpara$alpha1
      vee1 <- GetMpara$vee1
      Sig <- GetMpara$Sig
      H01 <- GetMpara$H01
      
      if((DiffSF(beta, prebeta, tau, pretau, gamma1, pregamma1,
               alpha1, prealpha1, vee1, prevee1, Sig, preSig, H01, preH01, epsilon) == 0) 
         || (iter == maxiter) || (!is.list(GetEfun)) || (!is.list(GetMpara))) {
        break
      }
      
    }
    
    if (iter == maxiter) {
      writeLines("program stops because of nonconvergence")
      convergence = 0
      result <- list(beta, tau, gamma1, alpha1, vee1, H01, Sig, iter, convergence)
      
      names(result) <- c("beta", "tau", "gamma1", "alpha1", "vee1", "H01", "Sig", 
                         "iter", "convergence")
      
      return(result)
    } else if (!is.list(GetEfun)) {
      writeLines("Something wrong in the E steps")
      beta <- NULL
      tau <- NULL
      gamma1 <- NULL
      alpha1 <- NULL
      vee1 <- NULL
      H01 <- NULL
      Sig <- NULL
      iter <- NULL
      result <- list(beta, tau, gamma1, alpha1, vee1, H01, Sig, iter)
      names(result) <- c("beta", "tau", "gamma1", "alpha1", "vee1", "H01", "Sig", "iter")
      return(result)
    } else if (!is.list(GetMpara)) {
      writeLines("Something wrong in the M steps")
      beta <- NULL
      tau <- NULL
      gamma1 <- NULL
      alpha1 <- NULL
      vee1 <- NULL
      H01 <- NULL
      Sig <- NULL
      iter <- NULL
      result <- list(beta, tau, gamma1, alpha1, vee1, H01, Sig, iter)
      names(result) <- c("beta", "tau", "gamma1", "alpha1", "vee1", "H01", "Sig", "iter")
      return(result)
    } else {
      
      convergence = 1
      
      GetEfun <- GetESF(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, 
                        X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix)
      
      FUNENW <- as.vector(GetEfun$FUNENW)
      FUNBENW <- as.matrix(GetEfun$FUNBENW)
      FUNBS <- as.matrix(GetEfun$FUNBS)
      FUNBW <- as.matrix(GetEfun$FUNBW)
      FUNWS <- as.vector(GetEfun$FUNWS)
      FUNBSENW <- as.matrix(GetEfun$FUNBSENW) 
      FUNEC <- as.matrix(GetEfun$FUNEC)
      FUNBEC <- as.matrix(GetEfun$FUNBEC)
      FUNBSEC <- as.matrix(GetEfun$FUNBSEC)
      FUNWEC <- as.matrix(GetEfun$FUNWEC)
      FUNWSEC <- as.matrix(GetEfun$FUNWSEC)
      FUNB <- as.matrix(GetEfun$FUNB)
      FUNW <- as.vector(GetEfun$FUNW)
      
      getcov <- getCovSF(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, 
                         X2, survtime, cmprsk, mdata, mdataS,
                         FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC,
                         FUNBSEC, FUNWEC, FUNWSEC,FUNB, FUNW)
      
      vcov <- getcov$vcov
      sebeta <- getcov$sebeta
      setau <- getcov$setau
      segamma1 <- getcov$segamma1
      sealpha1 <- getcov$sealpha1
      sevee1 <- getcov$sevee1
      seSig <- getcov$seSig
      
      long <- all.vars(long.formula)
      survival <- all.vars(surv.formula)
      
      names(beta) <- c("intercept", long[-1])
      variance.form <- all.vars(variance.formula)
      names(tau) <- paste0(c("intercept", variance.form), "_var")
      
      names(gamma1) <- paste0(survival[-(1:2)], "_1")
      
      PropComp <- as.data.frame(table(cdata[, survival[2]]))
      
      LongOut <- long[1]
      LongX <- paste0(long[-1], collapse = "+")
      FunCall_long <- as.formula(paste(LongOut, LongX, sep = "~"))
      
      LongVarOut <- "log(sigma^2)"
      LongVarX <- paste0(variance.form, collapse = "+") 
      FunCall_longVar <- as.formula(paste(LongVarOut, LongVarX, sep = "~"))
      
      SurvOut <- paste0("Surv(", survival[1], ",", survival[2], ")")
      SurvX <- paste0(survival[-(1:2)], collapse = "+")
      FunCall_survival <- as.formula(paste(SurvOut, SurvX, sep = "~"))
      
      ## return the joint modelling result
      mycall <- match.call()
      
      result <- list(beta, tau, gamma1, alpha1, vee1, H01, Sig, iter, convergence, 
                     vcov, sebeta, setau, segamma1, sealpha1, sevee1, seSig, CompetingRisk,
                     quadpoint, ydata, cdata, PropComp, FunCall_long, FunCall_longVar,
                     FunCall_survival, random, mycall)
      
      names(result) <- c("beta", "tau", "gamma1", "alpha1", "vee1", "H01", "Sig", 
                         "iter", "convergence", "vcov", "sebeta", "setau", "segamma1", 
                         "sealpha1", "sevee1", "seSig", "CompetingRisk", "quadpoint",
                         "ydata", "cdata", "PropEventType", "LongitudinalSubmodelmean",
                         "LongitudinalSubmodelvariance", "Survivalsubmodel", "random",
                         "call")
      
      class(result) <- "JMH"
      
      return(result)
    }
    
    
    
    
  }

  

    
}