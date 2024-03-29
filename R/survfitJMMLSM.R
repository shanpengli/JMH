##' @title Prediction in Joint Models
##' @name survfitJMMLSM
##' @aliases survfitJMMLSM
##' @description This function computes the conditional probability of 
##' surviving later times than the last observed time for which a longitudinal 
##' measurement was available.
##' @param object an object inheriting from class \code{JMMLSM}.
##' @param seed a random seed number to proceed non-parametric bootstrap. Default is 100.
##' @param ynewdata a data frame that contains the longitudinal and covariate information for the subjects 
##' for which prediction of survival probabilities is required.
##' @param cnewdata a data frame that contains the survival and covariate information for the subjects 
##' for which prediction of survival probabilities is required.
##' @param u a numeric vector of times for which prediction survival probabilities are to be computed.
##' @param Last.time a numeric vector or character string. This specifies the known time at which each of 
##' the subjects in cnewdata was known to be alive. If NULL, then this is automatically taken as the 
##' survival time of each subject. If a numeric vector, then it is assumed to be greater than or equals to the 
##' last available longitudinal time point for each subject. If a character string, then it should be 
##' a variable in cnewdata.
##' @param obs.time a character string of specifying a longitudinal time variable in ynewdata.
##' @param method a character string specifying the type of probability approximation; if \code{Laplace}, then a first order estimator is computed.
##' If \code{GH}, then the standard Gauss-Hermite quadrature is used instead.
##' @param quadpoint number of quadrature points used for estimating conditional probabilities 
##' when \code{method = "GH"}. Default is NULL. If \code{method = "GH"}, then 15 is used.
##' @param ... further arguments passed to or from other methods. 
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}}
##' @export
##' 
survfitJMMLSM <- function(object, seed = 100, ynewdata = NULL, cnewdata = NULL, 
                               u = NULL, Last.time = NULL, obs.time = NULL, method = c("Laplace", "GH"), quadpoint = NULL, ...) {
  if (!inherits(object, "JMMLSM"))
    stop("Use only with 'JMMLSM' objects.\n")
  if (is.null(ynewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(cnewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(u)) 
    stop("Please specify the future time for dynamic prediction.")   
  if (!method %in% c("Laplace", "GH"))
    stop("Please specify a method for probability approximation: Laplace or GH.")
  if (!is.vector(u)) 
    stop("u must be vector typed.")
  if (is.null(quadpoint)) {
    quadpoint <- object$quadpoint
  }
  if (is.null(obs.time)) {
    stop("Please specify a vector that represents the time variable from ydatanew.")
  } else {
    if (!obs.time %in% colnames(ynewdata)) {
      stop(paste0(obs.time, " is not found in ynewdata."))
    }
  }
  
  bvar <- all.vars(object$random)
  ID <- bvar[length(bvar)]
  if (!(ID %in% colnames(ynewdata)))
    stop(paste("The ID variable", ID, "is not found in ynewdata."))
  if (!(ID %in% colnames(cnewdata)))
    stop(paste("The ID variable", ID, "is not found in cnewdata."))
  
  ynewdata <- ynewdata[, colnames(object$ydata)]
  cnewdata <- cnewdata[, colnames(object$cdata)]
  
  ydata2 <- rbind(object$ydata, ynewdata)
  cdata2 <- rbind(object$cdata, cnewdata)
  
  variance.formula <- as.formula(paste("", object$LongitudinalSubmodelvariance[3], sep = "~"))
  getdum <- getdummy(long.formula = object$LongitudinalSubmodelmean,
                     surv.formula = object$SurvivalSubmodel,
                     variance.formula = variance.formula,
                     random = object$random, ydata = ydata2, cdata = cdata2)
  
  
  ydata.mean <- getdum$ydata.mean
  ydata.variance <- getdum$ydata.variance
  cdata2 <- getdum$cdata
  
  Yvar <- colnames(ydata.mean)[-1]
  Cvar <- colnames(cdata2)[-1]
  bvar <- all.vars(object$random)
  
  ny <- nrow(ynewdata)
  nc <- nrow(cnewdata)
  Ny <- nrow(ydata2)
  Nc <- nrow(cdata2)
  
  Sig <- object$Sig
  p1a <- ncol(Sig) - 1
  quadmethod <- object$method
  getGH <- GetGHmatrix(quadpoint = quadpoint, p1a = p1a)
  
  xsmatrix <- getGH$xsmatrix
  wsmatrix <- getGH$wsmatrix
  
  nsig <- p1a + 1
  
  ynewdata.mean <- ydata.mean[c((Ny-ny+1):Ny), ]
  ynewdata.variance <- ydata.variance[c((Ny-ny+1):Ny), ]
  cnewdata <- cdata2[c((Nc-nc+1):Nc), ]
  
  if (length(bvar) > 1) bvar1 <- bvar[1:(length(bvar) - 1)]
  yID <- unique(ynewdata.mean[, ID])
  N.ID <- length(yID)
  cID <- cnewdata[, ID]
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cnewdata.")
  }
  
  if (!is.null(Last.time)) {
    if (is.character(Last.time)) {
      if (Last.time %in% colnames(cnewdata)) {
        Last.time <- cnewdata[, Last.time]
      } else {
        stop(paste(Last.time, "is not found in cnewdata."))
      }
    }
    if (is.numeric(Last.time) && (length(Last.time) != nrow(cnewdata)))
      stop("The last.time vector does not match cnewdata.")
  } else {
    Last.time <- cnewdata[, Cvar[1]]
  }
  
  Pred <- list()
  CompetingRisk <- object$CompetingRisk
  if (object$CompetingRisk) {
    
    beta <- object$beta
    tau <- object$tau
    gamma1 <- object$gamma1
    gamma2 <- object$gamma2
    alpha1 <- object$alpha1
    alpha2 <- object$alpha2
    nu1 <- object$vee1
    nu2 <- object$vee2
    H01 <- object$H01
    H02 <- object$H02
    Sig <- object$Sig
    
    Predraw1 <- matrix(0, nrow = nrow(cnewdata), ncol = length(u))
    Predraw2 <- matrix(0, nrow = nrow(cnewdata), ncol = length(u))
    y.obs <- list()
    lengthu <- length(u)
    for (j in 1:N.ID) {
      subNDy.mean <- ynewdata.mean[ynewdata.mean[, ID] == yID[j], ]
      subNDy.variance <- ynewdata.variance[ynewdata.variance[, ID] == yID[j], ]
      subNDc <- cnewdata[cnewdata[, ID] == yID[j], ]
      y.obs[[j]] <- data.frame(ynewdata[ynewdata[, ID] == yID[j], c(obs.time, Yvar[1])])
      
      s <-  as.numeric(Last.time[j])
      CH01 <- CH(H01, s)
      CH02 <- CH(H02, s)
      Y <- subNDy.mean[, Yvar[1]]
      X <- subNDy.mean[, Yvar[2:length(Yvar)]]
      X <- as.matrix(X)
      W <- subNDy.variance[, -1]
      W <- as.matrix(W)
      if (nsig == 2) {
        Z <- matrix(1, ncol = 1, nrow = length(Y))
      } else {
        Z <- data.frame(1, subNDy.mean[, bvar1])
        Z <- as.matrix(Z)
      }
      X2 <- as.matrix(subNDc[1, Cvar[3:length(Cvar)]])

      if (method == "GH") {
        for (jj in 1:lengthu) {
          ## calculate the CIF
          
          if (quadmethod == "standard") {
            CIF <- getECIF(beta, tau, gamma1, gamma2, alpha1, alpha2, nu1,
                           nu2, Sig, Z, X, W, Y, as.vector(X2), H01, H02,
                           xsmatrix, wsmatrix, CH01, CH02, s, u[jj])
          } else {
            
            data <- list(Y, X, Z, W, X2, CH01, CH02, beta, tau, gamma1, gamma2, alpha1, alpha2, nu1, nu2, Sig)
            names(data) <- c("Y", "X", "Z", "W", "X2", "CH01", "CH02", "beta", "tau",
                             "gamma1", "gamma2", "alpha1", "alpha2", "nu1", "nu2", "Sig")
            opt <- optim(rep(0, nsig), logLikCR, data = data, method = "BFGS", hessian = TRUE)
            Posmean <- opt$par
            PosCov <- solve(opt$hessian)
            
            CIF <- getECIFad(beta, tau, gamma1, gamma2, alpha1, alpha2, nu1,
                             nu2, Sig, Z, X, W, Y, as.vector(X2), H01, H02,
                             xsmatrix, wsmatrix, CH01, CH02, s, u[jj], Posmean, PosCov)
            
          }
          
          P1us <- CIF$CIF1
          P2us <- CIF$CIF2
          
          Predraw1[j, jj] <- P1us
          Predraw2[j, jj] <- P2us
        }
      } else {
        data <- list(Y, X, Z, W, X2, CH01, CH02, beta, tau, gamma1, gamma2, alpha1, alpha2, nu1, nu2, Sig)
        names(data) <- c("Y", "X", "Z", "W", "X2", "CH01", "CH02", "beta", "tau",
                         "gamma1", "gamma2", "alpha1", "alpha2", "nu1", "nu2", "Sig")
        opt <- optim(rep(0, nsig), logLikCR, data = data, method = "BFGS", hessian = TRUE)
        meanbw <- opt$par
        for (jj in 1:lengthu) {
          ## calculate the CIF
          CIF1 <- CIF1.CR(data, H01, H02, s, u[jj], meanbw)
          P1us <- Pk.us(CIF1, data, meanbw)
          Predraw1[j, jj] <- P1us
          CIF2 <- CIF2.CR(data, H01, H02, s, u[jj], meanbw)
          P2us <- Pk.us(CIF2, data, meanbw)
          Predraw2[j, jj] <- P2us
        }
        quadpoint = NULL
      }
    }
    for (jj in 1:N.ID) {
      Pred[[jj]] <- data.frame(u, Predraw1[jj, ], Predraw2[jj, ])
      colnames(Pred[[jj]]) <- c("times", "CIF1", "CIF2")
    }
    
  } else {
    
    Predraw <- matrix(0, nrow = nrow(cnewdata), ncol = length(u))
    beta <- object$beta
    tau <- object$tau
    gamma <- object$gamma1
    alpha <- object$alpha1
    nu <- object$vee1
    H01 <- object$H01
    Sig <- object$Sig
    
    y.obs <- list()
    lengthu <- length(u)
    
    for (j in 1:N.ID) {
      subNDy.mean <- ynewdata.mean[ynewdata.mean[, ID] == yID[j], ]
      subNDy.variance <- ynewdata.variance[ynewdata.variance[, ID] == yID[j], ]
      subNDc <- cnewdata[cnewdata[, ID] == yID[j], ]
      y.obs[[j]] <- data.frame(ynewdata[ynewdata[, ID] == yID[j], c(obs.time, Yvar[1])])
      
      CH0 <- CH(H01, Last.time[j])
      CH0u <- vector()
      
      for (jj in 1:lengthu) {
        CH0u[jj] <- CH(H01, u[jj])
      }
      Y <- subNDy.mean[, Yvar[1]]
      X <- subNDy.mean[, Yvar[2:length(Yvar)]]
      X <- as.matrix(X)
      W <- subNDy.variance[, -1]
      W <- as.matrix(W)
      if (nsig == 2) {
        Z <- matrix(1, ncol = 1, nrow = length(Y))
      } else {
        Z <- data.frame(1, subNDy.mean[, bvar1])
        Z <- as.matrix(Z)
      }
      X2 <- as.matrix(subNDc[1, Cvar[3:length(Cvar)]])
      
      if (method == "Laplace") {
        ## find out E(theta_i)
        data <- list(Y, X, Z, W, X2, CH0, beta, tau, gamma, alpha, nu, Sig)
        names(data) <- c("Y", "X", "Z", "W", "X2", "CH0", "beta", "tau", "gamma", "alpha", "nu", "Sig")
        opt <- optim(rep(0, nsig), logLik, data = data, method = "BFGS", hessian = TRUE)
        meanbw <- opt$par
        for (jj in 1:lengthu) {
          Pi <- P.us(data, CH0u[jj], meanbw)
          Predraw[j, jj] <- 1 - Pi
        }
        quadpoint <- NULL
      } else {
        for (jj in 1:lengthu) {
          
          if (quadmethod == "standard") {
            Predraw[j, jj] <- getES(beta, tau, gamma, alpha, nu, Sig, Z, X, W, Y, 
                                    as.vector(X2), xsmatrix, wsmatrix, CH0, CH0u[jj])
          } else {
            
            data <- list(Y, X, Z, W, X2, CH0, beta, tau, gamma, alpha, nu, Sig)
            names(data) <- c("Y", "X", "Z", "W", "X2", "CH0", "beta", "tau",
                             "gamma1", "alpha1", "nu1", "Sig")
            opt <- optim(rep(0, nsig), logLik, data = data, method = "BFGS", hessian = TRUE)
            Posmean <- opt$par
            PosCov <- solve(opt$hessian)
            
            Predraw[j, jj] <- getESad(beta, tau, gamma, alpha, nu, Sig, Z, X, W, Y, 
                                    as.vector(X2), xsmatrix, wsmatrix, CH0, CH0u[jj], Posmean, PosCov)
            
          }
        }
      }
      
    }
    for (jj in 1:N.ID) {
      Pred[[jj]] <- data.frame(u, Predraw[jj, ])
      colnames(Pred[[jj]]) <- c("times", "PredSurv")
    }
    
  }
  names(y.obs) <- names(Pred) <- yID
  Last.time <- data.frame(cID, Last.time)
  colnames(Last.time)[1] <- ID
  sum <- list(Pred = Pred, Last.time = Last.time, y.obs = y.obs, method = method, quadpoint = quadpoint,
              CompetingRisk = CompetingRisk, quadmethod = quadmethod)
  class(sum) <- "survfitJMMLSM"
  sum
  
}