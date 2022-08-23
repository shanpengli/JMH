##' @title Prediction in Joint Models
##' @name survfit5JMMLSM
##' @aliases survfit5JMMLSM
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
##' @param B the number of bootstrap samples to be generated. Default is 100.
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
survfit5JMMLSM <- function(object, seed = 100, ynewdata = NULL, cnewdata = NULL, 
                           u = NULL, B = 100, method = c("Laplace", "GH"), quadpoint = NULL, ...) {
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

  bvar <- all.vars(object$random)
  if (!(bvar[length(bvar)] %in% colnames(ynewdata)))
    stop(paste("The ID variable", bvar[length(bvar)], "is not found in ynewdata."))
  if (!(bvar[length(bvar)] %in% colnames(cnewdata)))
    stop(paste("The ID variable", bvar[length(bvar)], "is not found in cnewdata."))
  
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
  if (p1a == 2) Sigb <- Sig[1:2, 1:2]
  if (p1a == 1) Sigb <- as.matrix(Sig[1, 1])
  
  getGH <- GetGHmatrix(quadpoint = quadpoint, Sigb = Sigb)
  
  xsmatrix <- getGH$xsmatrix
  wsmatrix <- getGH$wsmatrix
  
  nsig <- p1a + 1
  
  ynewdata.mean <- ydata.mean[c((Ny-ny+1):Ny), ]
  ynewdata.variance <- ydata.variance[c((Ny-ny+1):Ny), ]
  cnewdata <- cdata2[c((Nc-nc+1):Nc), ]
  
  if (length(bvar) > 1) bvar1 <- bvar[1:(length(bvar) - 1)]
  
  yID <- unique(ynewdata.mean[, bvar[length(bvar)]])
  N.ID <- length(yID)
  cID <- cnewdata[, bvar[length(bvar)]]
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cnewdata.")
  }
  
  ID <- bvar[length(bvar)]
  variance.formula <- as.formula(paste("", object$LongitudinalSubmodelvariance[3], sep = "~"))
  
  set.seed(seed)
  Predraw1 <- array(rep(NA, B*nrow(cnewdata)*length(u)), dim=c(B, nrow(cnewdata), length(u)))
  Predraw2 <- array(rep(NA, B*nrow(cnewdata)*length(u)), dim=c(B, nrow(cnewdata), length(u)))
  CompetingRisk <- object$CompetingRisk
  for (b in 1:B) {
    
    btsp <- sample(c(1:nrow(object$cdata)), replace = TRUE)
    cdata.btsp <- object$cdata[btsp, ]
    ydata <- object$ydata
    colnames(cdata.btsp)[-which(colnames(cdata.btsp) == ID)] <- 
      paste0(colnames(cdata.btsp)[-which(colnames(cdata.btsp) == ID)], ".survival")
    colnames(ydata)[-which(colnames(ydata) == ID)] <- 
      paste0(colnames(ydata)[-which(colnames(ydata) == ID)], ".longitudinal")
    
    ydata.btsp <- dplyr::left_join(cdata.btsp, ydata, by = ID)
    Nonyname <- colnames(cdata.btsp)[-which(colnames(cdata.btsp) == ID)]
    ydata.btsp <- ydata.btsp[, !colnames(ydata.btsp) %in% Nonyname]
    
    colnames(cdata.btsp) <- gsub(pattern = ".survival", replacement = "", 
                                 x = colnames(cdata.btsp))
    colnames(ydata.btsp) <- gsub(pattern = ".longitudinal", replacement = "", 
                                 x = colnames(ydata.btsp))
    
    mdata <- as.data.frame(table(object$ydata[, ID]))
    colnames(mdata) <- c(ID, "Freq00")
    mdataboots <- dplyr::left_join(cdata.btsp, mdata, by = ID)
    mdataboots <- mdataboots[, c(ID, "Freq00")]
    
    IID <- vector()
    for (i in 1:nrow(mdataboots)) {
      num <- rep(i, mdataboots[i, 2])
      IID <- c(IID, num)
    }
    ydata.btsp[, ID] <- IID
    cdata.btsp[, ID] <- c(1:nrow(mdataboots))
    
    fitboot <- try(JMMLSM(cdata = cdata.btsp, ydata = ydata.btsp,
                          long.formula = object$LongitudinalSubmodelmean,
                          surv.formula = object$SurvivalSubmodel,
                          variance.formula = variance.formula,
                          quadpoint = object$quadpoint, random = object$random), silent = TRUE)
    
    
    if ('try-error' %in% class(fitboot)) {
      next
    } else if (fitboot$iter == maxiter) {
      next
    } else {
      if (fitboot$CompetingRisk == CompetingRisk) {
        beta <- fitboot$beta
        tau <- fitboot$tau
        gamma1 <- fitboot$gamma1
        gamma2 <- fitboot$gamma2
        alpha1 <- fitboot$alpha1
        alpha2 <- fitboot$alpha2
        nu1 <- fitboot$vee1
        nu2 <- fitboot$vee2
        H01 <- fitboot$H01
        H02 <- fitboot$H02
        Sig <- fitboot$Sig
        
        y.obs <- list()
        lengthu <- length(u)
        Last.time <- cnewdata[, c(bvar[length(bvar)], Cvar[1])]
        for (j in 1:N.ID) {
          subNDy.mean <- ynewdata.mean[ynewdata.mean[, bvar[length(bvar)]] == yID[j], ]
          subNDy.variance <- ynewdata.variance[ynewdata.variance[, bvar[length(bvar)]] == yID[j], ]
          subNDc <- cnewdata[cnewdata[, bvar[length(bvar)]] == yID[j], ]
          y.obs[[j]] <- data.frame(subNDy.mean[, c(bvar[1], Yvar[1])])
          
          s <-  as.numeric(subNDc[1, Cvar[1]])
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
          data <- list(Y, X, Z, W, X2, CH01, CH02, beta, tau, gamma1, gamma2, alpha1, alpha2, nu1, nu2, Sig)
          names(data) <- c("Y", "X", "Z", "W", "X2", "CH01", "CH02", "beta", "tau",
                           "gamma1", "gamma2", "alpha1", "alpha2", "nu1", "nu2", "Sig")
          
          opt <- optim(rep(0, nsig), logLikCR, data = data, method = "BFGS", hessian = TRUE)
          meanbw <- opt$par
          
          if (method == "GH") {
            for (jj in 1:lengthu) {
              ## calculate the CIF
              CIF <- getECIF(beta, tau, gamma1, gamma2, alpha1, alpha2, nu1,
                             nu2, Sig, Z, X, W, Y, as.vector(X2), H01, H02,
                             xsmatrix, wsmatrix, CH01, CH02, s, u[jj])
              P1us <- CIF$CIF1
              P2us <- CIF$CIF2
              
              Predraw1[b, j, jj] <- P1us
              Predraw2[b, j, jj] <- P2us
            }
          } else {
            for (jj in 1:lengthu) {
              ## calculate the CIF
              CIF1 <- CIF1.CR(data, H01, H02, s, u[jj], meanbw)
              P1us <- Pk.us(CIF1, data, meanbw)
              Predraw1[b, j, jj] <- P1us
              CIF2 <- CIF2.CR(data, H01, H02, s, u[jj], meanbw)
              P2us <- Pk.us(CIF2, data, meanbw)
              Predraw2[b, j, jj] <- P2us
            }
          }
          
          
        }
        names(y.obs) <- yID
        
      } else {
        next
      }
    } 
    

  }
  
  if (CompetingRisk) {
    Pred <- list()
    for (jj in 1:N.ID) {
      allPi1 <- Predraw1[, jj, ]
      allPi2 <- Predraw2[, jj, ]
      allPi1 <- allPi1[complete.cases(allPi1), ]
      allPi2 <- allPi2[complete.cases(allPi2), ]
      
      subCP1 <- as.data.frame(matrix(0, nrow = lengthu, ncol = 5))
      colnames(subCP1) <- c("times", "Mean", "Median", "95%HDLower", "95%HDUpper")
      for (b in 1:lengthu) {
        subCP1[b, 1] <- u[b]
        subCP1[b, 2] <- mean(allPi1[, b])
        subCP1[b, 3] <- median(allPi1[, b])
        subCP1[b, 4] <- Hmisc::hdquantile(allPi1[, b], probs = 0.025)
        subCP1[b, 5] <- Hmisc::hdquantile(allPi1[, b], probs = 0.975)
      }
      
      subCP2 <- as.data.frame(matrix(0, nrow = lengthu, ncol = 5))
      colnames(subCP2) <- c("times", "Mean", "Median", "95%HDLower", "95%HDUpper")
      for (b in 1:lengthu) {
        subCP2[b, 1] <- u[b]
        subCP2[b, 2] <- mean(allPi2[, b])
        subCP2[b, 3] <- median(allPi2[, b])
        subCP2[b, 4] <- Hmisc::hdquantile(allPi2[, b], probs = 0.025)
        subCP2[b, 5] <- Hmisc::hdquantile(allPi2[, b], probs = 0.975)
      }
      
      subCP <- list(subCP1, subCP2)
      names(subCP) <- c("Cumulative incidence probabilities for type 1 failure",
                        "Cumulative incidence probabilities for type 2 failure")
      
      Pred[[jj]] <- subCP 
    }
    
  }
  
  names(Pred) <- yID
  sum <- list(Pred = Pred, Last.time = Last.time, y.obs = y.obs, method = method, quadpoint = quadpoint,
              B = B, CompetingRisk = CompetingRisk)
  
  class(sum) <- "survfitJMMLSMboot"
  sum
  
}