##' @export
##'

JMMLSM <- function(cdata, ydata,
                long.formula,
                surv.formula,
                variance.var, maxiter = 1000, epsilon = 1e-04, 
                quadpoint = 10, ID, RE = NULL,
                model = "interslope", print.para = FALSE) {
  
  #mdata <- as.data.frame(table(ydata[, ID]))
  #colnames(mdata) <- c(ID, "ni")
  
  # getinit <- Getinit(cdata = cdata, ydata = ydata, long.formula = long.formula, 
  #                    surv.formula = surv.formula, variance.var = variance.var,
  #                    model = model, ID = ID, RE = RE)
  
  getinit <- GetinitFake(cdata, ydata, long.formula, surv.formula, variance.var,
                         model, ID, RE)
  
  cdata <- getinit$cdata
  
  ## initialize parameters
  
  getriskset <- Getriskset(cdata = cdata, surv.formula = surv.formula)
  
  ## number of distinct survival time
  H01 <- as.matrix(getriskset$tablerisk1)
  H02 <- as.matrix(getriskset$tablerisk2)
  
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
  Sigb <- Sig[1:2, 1:2]
  
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
    
    result <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, 
                   H02, Sig, iter, convergence)
    
    names(result) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1",
                       "vee2", "H01", "H02", "Sig", "iter", "convergence")
    
    return(result)
  }
    
}