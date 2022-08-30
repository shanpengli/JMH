##' @title Prediction in Joint Models
##' @name survfit2JMMLSM
##' @aliases survfit2JMMLSM
##' @description This function computes the conditional probability of 
##' surviving later times than the last observed time for which a longitudinal 
##' measurement was available.
##' @param object an object inheriting from class \code{JMMLSM}.
##' @param seed a random seed number to proceed Monte Carlo simulation. Default is 100.
##' @param ynewdata a data frame that contains the longitudinal and covariate information for the subjects 
##' for which prediction of survival probabilities is required.
##' @param cnewdata a data frame that contains the survival and covariate information for the subjects 
##' for which prediction of survival probabilities is required.
##' @param u a numeric vector of times for which prediction survival probabilities are to be computed.
##' @param M the number of Monte Carlo samples to be generated. Default is 200.
##' @param burn.in The proportion of M to discard as burn-in. Default is 0.2.
##' @param simulate logical; if \code{TRUE}, a Monte Carlo approach is used to estimate conditional probabilities. 
##' Otherwise, Gauss-Hermite quadrature rule is used for numerical integration to estimate instead. 
##' Default is \code{TRUE}.
##' @param quadpoint number of quadrature points used for estimating conditional probabilities when \code{simulate = FALSE}. Default is 20.
##' @param ... further arguments passed to or from other methods. 
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}}
##' @export
##' 
survfit2JMMLSM <- function(object, seed = 100, ynewdata = NULL, cnewdata = NULL, 
                           u = NULL, M = 200, burn.in = 0.2, simulate = TRUE, quadpoint = NULL, ...) {
  if (!inherits(object, "JMMLSM"))
    stop("Use only with 'JMMLSM' objects.\n")
  if (is.null(ynewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(cnewdata))
    stop("New longitudinal data for dynamic prediction is needed.")
  if (is.null(u)) 
    stop("Please specify the future time for dynamic prediction.")    
  if (!is.vector(u)) 
    stop("u must be vector typed.")
  if (!object$CompetingRisk) {
    H01 <- object$H01
    if (max(u) > H01[nrow(H01), 1])
      stop(paste("The current joint model cannot predict the conditional 
         survival probabilities later than the last observed time of the object. 
               The last observed time is", max(H01[, 1])))
  } else {
    H01 <- object$H01
    H02 <- object$H02
    # if (max(u) > H01[nrow(H01), 1] | max(u) > H02[nrow(H02), 1])
    #   stop(paste("The current joint model cannot predict the conditional 
    #      survival probabilities later than the last observed time of the object. 
    #            The last observed time for risk 1 and 2 is", max(H01[, 1]), "and", max(H02[, 1])))
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
  
  ynewdata.mean <- ydata.mean[c((Ny-ny+1):Ny), ]
  ynewdata.variance <- ydata.variance[c((Ny-ny+1):Ny), ]
  cnewdata <- cdata2[c((Nc-nc+1):Nc), ]
  
  ynewdata.mean2 <- ydata.mean[c(1:(Ny-ny)), ]
  ynewdata.variance2 <- ydata.variance[c(1:(Ny-ny)), ]
  cnewdata2 <- cdata2[c(1:(Nc-nc)), ]
  
  GetVar <- GetVarSurvfit(cdata = cnewdata2, ydata.mean = ynewdata.mean2,
                          ydata.variance = ynewdata.variance2, random = object$random)
  
  ## dynamic prediction 
  ## Monte Carlo simulation
  ID <- unique(ynewdata.mean[, bvar[length(bvar)]])
  N.ID <- length(ID)
  cID <- cnewdata[, bvar[length(bvar)]]
  if (prod(ID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cnewdata.")
  }
  
  if (simulate) {
    
    Sig <- object$Sig
    p1a <- ncol(Sig) - 1
    if (p1a == 2) Sigb <- Sig[1:2, 1:2]
    if (p1a == 1) Sigb <- as.matrix(Sig[1, 1])
    
    getGH <- GetGHmatrix(quadpoint = quadpoint, Sigb = Sigb)
    
    xsmatrix <- getGH$xsmatrix
    wsmatrix <- getGH$wsmatrix
    
    
    set.seed(seed)
    nbeta <- length(object$beta)
    ntau <- length(object$tau)
    ngamma <- length(object$gamma1)
    nalpha <- length(object$alpha1)
    nnu <- 1
    nsig <- nrow(object$Sig)
    lengthu <- length(u)
    if (!object$CompetingRisk) {
      Psi <- c(object$beta, object$tau, object$gamma1, object$alpha1, object$vee1)
      for (l in 1:nsig) Psi <- c(Psi, object$Sig[l, l])
      if (nsig == 2) Psi <- c(Psi, object$Sig[1, 2])
      if (nsig == 3) {
        Psi <- c(Psi, object$Sig[1, 2])
        Psi <- c(Psi, object$Sig[2, 3])
        Psi <- c(Psi, object$Sig[1, 3])
      }
      covPsi <- vcov(object)
      Psi.MC <- mvrnorm(n = M, Psi, covPsi, tol = 1e-6, empirical = FALSE)
      Pred <- list()
      y.obs <- list()
      if (length(bvar) > 1) bvar1 <- bvar[1:(length(bvar) - 1)]
      for (j in 1:N.ID) {
        subNDy.mean <- ynewdata.mean[ynewdata.mean[, bvar[length(bvar)]] == ID[j], ]
        subNDy.variance <- ynewdata.variance[ynewdata.variance[, bvar[length(bvar)]] == ID[j], ]
        subNDc <- cnewdata[cnewdata[, bvar[length(bvar)]] == ID[j], ]
        y.obs[[j]] <- data.frame(subNDy.mean[, c(1, 2)])
        allPi <- matrix(0, ncol = length(u), nrow = M)
        CH0 <- CH(H01, subNDc[, Cvar[1]])
        CH0u <- vector()
        for (jj in 1:lengthu) {
          CH0u[jj] <- CH(H01, u[jj])
        }
        Y <- subNDy.mean[, 2]
        X <- subNDy.mean[, -c(1:2)]
        X <- as.matrix(X)
        W <- subNDy.variance[, -1]
        W <- as.matrix(W)
        if (nsig == 2) {
          Z <- matrix(1, ncol = 1, nrow = length(Y))
        } else {
          Z <- data.frame(1, subNDy.mean[, bvar1])
          Z <- as.matrix(Z)
        }
        X2 <- as.matrix(subNDc[, -c(1:3)])
        for (i in 1:M) {
          ##1. draw Psi
          psil <- Psi.MC[i, ]
          betal <- psil[1:nbeta]
          taul <- psil[(nbeta+1):(nbeta+ntau)]
          gammal <- psil[(nbeta+ntau+1):(nbeta+ntau+ngamma)]
          alphal <- psil[(nbeta+ntau+ngamma+1):(nbeta+ntau+ngamma+nalpha)]
          nul <- psil[nbeta+ntau+ngamma+nalpha+1]
          Sigl <- matrix(0, ncol = nsig, nrow = nsig)
          for (l in 1:nsig) Sigl[l, l] <- psil[nbeta+ntau+ngamma+nalpha+1+l]
          if (nsig == 2) Sigl[1, 2] <- Sigl[2, 1] <- psil[nbeta+ntau+ngamma+nalpha+1+nsig+1]
          if (nsig == 3) {
            Sigl[1, 2] <- Sigl[2, 1] <- psil[nbeta+ntau+ngamma+nalpha+1+nsig+1]
            Sigl[2, 3] <- Sigl[3, 2] <- psil[nbeta+ntau+ngamma+nalpha+1+nsig+2]
            Sigl[1, 3] <- Sigl[3, 1] <- psil[nbeta+ntau+ngamma+nalpha+1+nsig+3]
          }
          data <- list(Y, X, Z, W, X2, CH0, betal, taul, gammal, alphal, nul, Sigl)
          names(data) <- c("Y", "X", "Z", "W", "X2", "CH0", "beta", "tau", "gamma", "alpha", "nu", "Sig")
          opt <- optim(rep(0, nsig), logLik, data = data, method = "BFGS", hessian = TRUE)
          meanb <- opt$par
          varb <- solve(opt$hessian)
          b.old <- meanb
          ## simulate new random effects estimates using the Metropolis Hastings algorithm
          propose.bl <- as.vector(mvtnorm::rmvt(1, delta = meanb, sigma = varb, df = 4))
          dmvt.old <- mvtnorm::dmvt(b.old, meanb, varb, df = 4, TRUE)
          dmvt.propose <- mvtnorm::dmvt(propose.bl, meanb, varb, df = 4, TRUE)
          logpost.old <- -logLik(data, b.old)
          logpost.propose <- -logLik(data, propose.bl)
          ratio <- min(exp(logpost.propose + dmvt.old - logpost.old - dmvt.propose), 1)
          if (runif(1) <= ratio) {
            bl = propose.bl
          } else {
            bl = b.old
          }
          for (jj in 1:lengthu) {
            Pi <- P.us(data, CH0u[jj], bl)
            allPi[i, jj] <- Pi
          }
        }
        allPi <- as.data.frame(allPi)
        colnames(allPi) <- u
        
        subCP <- as.data.frame(matrix(0, nrow = length(u), ncol = 5))
        colnames(subCP) <- c("times", "Mean", "Median", "95%Lower", "95%Upper")
        for (b in 1:length(u)) {
          subCP[b, 1] <- u[b]
          subCP[b, 2] <- 1 - mean(allPi[, b])
          subCP[b, 3] <- median(1 - allPi[, b])
          subCP[b, 4] <- quantile(1 - allPi[, b], probs = 0.025)
          subCP[b, 5] <- quantile(1 - allPi[, b], probs = 0.975)
        }
        Pred[[j]] <- subCP 
      }
      names(Pred) <- ID
      sum <- list()
      sum$Pred <- Pred
      class(sum) <- "survfitJMMLSM"
      sum$Last.time <- cnewdata[, c(bvar[length(bvar)], Cvar[1])]
      sum$M <- M
      names(y.obs) <- ID
      sum$y.obs <- y.obs
      sum$CompetingRisk <- FALSE
      sum$simulate <- simulate
      sum
      
    } else {
      
      
      Psi <- c(object$beta, object$tau, object$gamma1, object$gamma2, 
               object$alpha1, object$alpha2, object$vee1, object$vee2)
      for (l in 1:nsig) Psi <- c(Psi, object$Sig[l, l])
      if (nsig == 2) Psi <- c(Psi, object$Sig[1, 2])
      if (nsig == 3) {
        Psi <- c(Psi, object$Sig[1, 2])
        Psi <- c(Psi, object$Sig[2, 3])
        Psi <- c(Psi, object$Sig[1, 3])
      }
      
      covPsi <- vcov(object)
      # covPsi <- object$vcov
      
      Psi.MC <- mvrnorm(n = M, Psi, covPsi, tol = 1e-6, empirical = FALSE)
      Pred <- list()
      y.obs <- list()
      if (length(bvar) > 1) bvar1 <- bvar[1:(length(bvar) - 1)]
      
      ## Modified Monte Carlo method
      Psi.init <- Psi
      H01.init <- object$H01
      H02.init <- object$H02
      
      allPi1 <- list()
      allPi2 <- list()
      
      for (i in 1:N.ID) {
        allPi1[[i]] <- matrix(0, ncol = length(u), nrow = M)
        allPi1[[i]] <- as.data.frame(allPi1[[i]])
        colnames(allPi1[[i]]) <- u
        allPi2[[i]] <- matrix(0, ncol = length(u), nrow = M)
        allPi2[[i]] <- as.data.frame(allPi2[[i]])
        colnames(allPi2[[i]]) <- u
      }
      
      names(allPi1) <- ID
      names(allPi2) <- ID
      
      pb = txtProgressBar(min = 1, max = M, initial = 1, style = 3) 
      
      posterior <- list()
      
      theta.modes <- matrix(0, nrow = N.ID, ncol = nsig)
      theta.var <- list()
      ## extract Empirical Bayes estimates 
      for (j in 1:N.ID) {
        subNDy.mean <- ynewdata.mean[ynewdata.mean[, bvar[length(bvar)]] == ID[j], ]
        subNDy.variance <- ynewdata.variance[ynewdata.variance[, bvar[length(bvar)]] == ID[j], ]
        subNDc <- cnewdata[cnewdata[, bvar[length(bvar)]] == ID[j], ]
        y.obs[[j]] <- data.frame(subNDy.mean[, c(bvar[1], Yvar[1])])
        
        s <-  as.numeric(subNDc[1, Cvar[1]])
        Y <- subNDy.mean[, 2]
        X <- subNDy.mean[, -c(1:2)]
        X <- as.matrix(X)
        W <- subNDy.variance[, -1]
        W <- as.matrix(W)
        if (nsig == 2) {
          Z <- matrix(1, ncol = 1, nrow = length(Y))
        } else {
          Z <- data.frame(1, subNDy.mean[, bvar1])
          Z <- as.matrix(Z)
        }
        X2 <- as.matrix(subNDc[, -c(1:3)])
        
        CH01 <- CH(H01.init, s)
        CH02 <- CH(H02.init, s)
        
        data <- list(Y, X, Z, W, X2, CH01, CH02, object$beta, object$tau, object$gamma1, 
                     object$gamma2, object$alpha1, object$alpha2, object$vee1, object$vee2, 
                     object$Sig)
        names(data) <- c("Y", "X", "Z", "W", "X2", "CH01", "CH02", "beta", "tau",
                         "gamma1", "gamma2", "alpha1", "alpha2", "nu1", "nu2", "Sig")
        opt <- optim(rep(0, nsig), logLikCR, data = data, method = "BFGS", hessian = TRUE)
        theta.modes[j, ] <- opt$par
        theta.var[[j]] <- solve(opt$hessian)
      }
      theta.old <- theta.new <- theta.modes
      
      for (i in 1:M) {

        ### 0. Set the initial estimator
        H01l <- H01.init
        H02l <- H02.init
        
        ###1. draw new parameters
        psil <- Psi.MC[i, ]
        betal <- psil[1:nbeta]
        taul <- psil[(nbeta+1):(nbeta+ntau)]
        gammal1 <- psil[(nbeta+ntau+1):(nbeta+ntau+ngamma)]
        gammal2 <- psil[(nbeta+ntau+ngamma+1):(nbeta+ntau+2*ngamma)]
        alphal1 <- psil[(nbeta+ntau+2*ngamma+1):(nbeta+ntau+2*ngamma+nalpha)]
        alphal2 <- psil[(nbeta+ntau+2*ngamma+nalpha+1):(nbeta+ntau+2*ngamma+2*nalpha)]
        nul1 <- psil[nbeta+ntau+2*ngamma+2*nalpha+1]
        nul2 <- psil[nbeta+ntau+2*ngamma+2*nalpha+2]
        Sigl <- matrix(0, ncol = nsig, nrow = nsig)
        for (l in 1:nsig) Sigl[l, l] <- psil[nbeta+ntau+2*ngamma+2*nalpha+2+l]
        if (nsig == 2) Sigl[1, 2] <- Sigl[2, 1] <- psil[nbeta+ntau+2*ngamma+2*nalpha+2+nsig+1]
        if (nsig == 3) {
          Sigl[1, 2] <- Sigl[2, 1] <- psil[nbeta+ntau+2*ngamma+2*nalpha+2+nsig+1]
          Sigl[2, 3] <- Sigl[3, 2] <- psil[nbeta+ntau+2*ngamma+2*nalpha+2+nsig+2]
          Sigl[1, 3] <- Sigl[3, 1] <- psil[nbeta+ntau+2*ngamma+2*nalpha+2+nsig+3]
        }
        
        ###1. Compute the expected value
        n <- nrow(object$cdata)
        p1a <- nsig-1
        CUH01 <- rep(0, n)
        CUH02 <- rep(0, n)
        HAZ01 <- rep(0, n)
        HAZ02 <- rep(0, n)
        CumuH01 <- cumsum(H01l[, 3])
        CumuH02 <- cumsum(H02l[, 3])
        
        Z <- GetVar$Z
        X1 <- GetVar$X1
        W <- GetVar$W
        Y <- GetVar$Y
        X2 <- GetVar$X2
        survtime <- GetVar$survtime
        cmprsk <- GetVar$cmprsk
        mdata <- GetVar$mdata
        mdataS <- GetVar$mdataS
        
        getHazard(CumuH01, CumuH02, survtime, cmprsk, H01l, H02l, CUH01, CUH02, HAZ01, HAZ02)
        
        EWik = getEWik(betal, taul, gammal1, gammal2, alphal1, alphal2, nul1, nul2, 
                       Sigl, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix,
                       wsmatrix, CUH01, CUH02, HAZ01, HAZ02)
        EWik <- EWik$FUNEC
        
        ###2. Generate a perturbation
        vil <- 4*rbeta(n, shape1 = 1/2, shape2 = 3/2)
        
        ###3. Compute Breslow's estimator
        subcdata <- cbind(survtime, cmprsk)
        subcdata <- as.matrix(subcdata)
        GetNewH0 <- GetBreslowInit(subcdata, vil)
        H01l <- GetNewH0$H01
        H02l <- GetNewH0$H02
        
        GetNewH0 <- updateH0(gammal1, gammal2, X2, survtime, cmprsk, EWik, vil, H01l, H02l)
        H01l <- GetNewH0$H01
        H02l <- GetNewH0$H02
        
        postbl <- matrix(NA, nrow = N.ID, ncol = nsig)
        for (j in 1:N.ID) {
          subNDy.mean <- ynewdata.mean[ynewdata.mean[, bvar[length(bvar)]] == ID[j], ]
          subNDy.variance <- ynewdata.variance[ynewdata.variance[, bvar[length(bvar)]] == ID[j], ]
          subNDc <- cnewdata[cnewdata[, bvar[length(bvar)]] == ID[j], ]
          y.obs[[j]] <- data.frame(subNDy.mean[, c(bvar[1], Yvar[1])])
          
          s <-  as.numeric(subNDc[1, Cvar[1]])
          Y <- subNDy.mean[, 2]
          X <- subNDy.mean[, -c(1:2)]
          X <- as.matrix(X)
          W <- subNDy.variance[, -1]
          W <- as.matrix(W)
          if (nsig == 2) {
            Z <- matrix(1, ncol = 1, nrow = length(Y))
          } else {
            Z <- data.frame(1, subNDy.mean[, bvar1])
            Z <- as.matrix(Z)
          }
          X2 <- as.matrix(subNDc[, -c(1:3)])
          
          CH01 <- CH(H01l, s)
          CH02 <- CH(H02l, s)
          
          data <- list(Y, X, Z, W, X2, CH01, CH02, betal, taul, gammal1, gammal2, alphal1, alphal2, nul1, nul2, Sigl)
          names(data) <- c("Y", "X", "Z", "W", "X2", "CH01", "CH02", "beta", "tau",
                           "gamma1", "gamma2", "alpha1", "alpha2", "nu1", "nu2", "Sig")

          ## simulate new random effects estimates using the Metropolis Hastings algorithm
          propose.theta <- as.vector(mvtnorm::rmvt(1, delta = theta.old[j, ], sigma = theta.var[[j]], df = 4))
          dmvt.old <- mvtnorm::dmvt(theta.old[j, ], propose.theta, theta.var[[j]], df = 4, TRUE)
          dmvt.propose <- mvtnorm::dmvt(propose.theta, theta.old[j, ], theta.var[[j]], df = 4, TRUE)
          logpost.old <- -logLikCR(data, theta.old[j, ])
          logpost.propose <- -logLikCR(data, propose.theta)
          ratio <- min(exp(logpost.propose + dmvt.old - logpost.old - dmvt.propose), 1)
          if (runif(1) <= ratio) {
            theta.new[j, ] = propose.theta
          } 
          postbl[j, ] <- theta.new[j, ]
          for (jj in 1:lengthu) {
            ## calculate the CIF
            CIF1 <- GetCIF1CR(data$gamma1, data$gamma2, data$alpha1, data$alpha2, 
                              data$nu1, data$nu2, as.vector(data$X2),
                              H01l, H02l, s, u[jj], theta.new[j, ], nrow(data$Sig))
            P1us <- Pk.us(CIF1, data, theta.new[j, ])
            if (P1us > 1) P1us <- 1
            allPi1[[j]][i, jj] <- P1us
            CIF2 <- GetCIF2CR(data$gamma1, data$gamma2, data$alpha1, data$alpha2, 
                              data$nu1, data$nu2, as.vector(data$X2),
                              H01l, H02l, s, u[jj], theta.new[j, ], nrow(data$Sig))
            P2us <- Pk.us(CIF2, data, theta.new[j, ])
            if (P2us > 1) P2us <- 1
            allPi2[[j]][i, jj] <- P2us
          }
        }
        theta.old <- theta.new
        ## pass the current sample parameters to the next iteration
        H01.init <- H01l
        H02.init <- H02l
        
        posterior[[i]] <- postbl
        
        setTxtProgressBar(pb,i) 
      }
      names(posterior) <- c(1:M)
      close(pb)
      
      for (j in 1:N.ID) {
        
        allPi1[[j]] <- allPi1[[j]][complete.cases(allPi1[[j]]), ]
        allPi2[[j]] <- allPi2[[j]][complete.cases(allPi2[[j]]), ]
        
        nEnd1 <- nrow(allPi1[[j]])
        nEnd2 <- nrow(allPi2[[j]])
        nStart1 <- floor(nEnd1*burn.in)+1
        nStart2 <- floor(nEnd2*burn.in)+1
        allPi1[[j]] <- allPi1[[j]][nStart1:nEnd1, ]
        allPi2[[j]] <- allPi2[[j]][nStart2:nEnd2, ]
        
        subCP1 <- as.data.frame(matrix(0, nrow = length(u), ncol = 5))
        colnames(subCP1) <- c("times", "Mean", "Median", "95%HDLower", "95%HDUpper")
        for (b in 1:length(u)) {
          subCP1[b, 1] <- u[b]
          subCP1[b, 2] <- mean(allPi1[[j]][, b])
          subCP1[b, 3] <- median(allPi1[[j]][, b])
          subCP1[b, 4] <- quantile(allPi1[[j]][, b], probs = 0.025)
          subCP1[b, 5] <- quantile(allPi1[[j]][, b], probs = 0.975)
        }
        
        subCP2 <- as.data.frame(matrix(0, nrow = length(u), ncol = 5))
        colnames(subCP2) <- c("times", "Mean", "Median", "95%HDLower", "95%HDUpper")
        for (b in 1:length(u)) {
          subCP2[b, 1] <- u[b]
          subCP2[b, 2] <- mean(allPi2[[j]][, b])
          subCP2[b, 3] <- median(allPi2[[j]][, b])
          subCP2[b, 4] <- quantile(allPi2[[j]][, b], probs = 0.025)
          subCP2[b, 5] <- quantile(allPi2[[j]][, b], probs = 0.975)
        }
        
        subCP <- list(subCP1, subCP2)
        names(subCP) <- c("Cumulative incidence probabilities for type 1 failure",
                          "Cumulative incidence probabilities for type 2 failure")
        
        Pred[[j]] <- subCP 
      }
      
      names(Pred) <- ID
      sum <- list()
      sum$Pred <- Pred
      class(sum) <- "survfitJMMLSM"
      sum$Last.time <- cnewdata[, c(bvar[length(bvar)], Cvar[1])]
      sum$M <- M
      names(y.obs) <- ID
      sum$y.obs <- y.obs
      sum$CompetingRisk <- TRUE
      sum$simulate <- simulate
      sum$posterior <- posterior
      sum
      
    }
  } else {
    
    nsig <- ncol(object$Sig)
    
    if (length(bvar) > 1) bvar1 <- bvar[1:(length(bvar) - 1)]
    if (!object$CompetingRisk) {
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
      Last.time <- cnewdata[, c(bvar[length(bvar)], Cvar[1])]
      CompetingRisk <- object$CompetingRisk
      
      for (j in 1:N.ID) {
        subNDy.mean <- ynewdata.mean[ynewdata.mean[, bvar[length(bvar)]] == ID[j], ]
        subNDy.variance <- ynewdata.variance[ynewdata.variance[, bvar[length(bvar)]] == ID[j], ]
        subNDc <- cnewdata[cnewdata[, bvar[length(bvar)]] == ID[j], ]
        y.obs[[j]] <- data.frame(subNDy.mean[, c(bvar[1], Yvar[1])])
        CH0 <- CH(H01, subNDc[, Cvar[1]])
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
        
        ## find out E(theta_i)
        data <- list(Y, X, Z, W, X2, CH0, beta, tau, gamma, alpha, nu, Sig)
        names(data) <- c("Y", "X", "Z", "W", "X2", "CH0", "beta", "tau", "gamma", "alpha", "nu", "Sig")
        opt <- optim(rep(0, nsig), logLik, data = data, method = "BFGS", hessian = TRUE)
        meanbw <- opt$par
        
        for (jj in 1:lengthu) {
          Pi <- P.us(data, CH0u[jj], meanbw)
          Predraw[j, jj] <- 1 - Pi
        }
        
      }
      names(y.obs) <- ID
      Pred <- list()
      for (jj in 1:N.ID) {
        Pred[[jj]] <- data.frame(u, Predraw[jj, ])
        colnames(Pred[[jj]]) <- c("times", "PredSurv")
      }
      names(Pred) <- ID
      sum <- list(Pred = Pred, Last.time = Last.time, y.obs = y.obs, simulate = simulate, 
                  quadpoint = quadpoint, CompetingRisk = CompetingRisk)
      class(sum) <- "survfitJMMLSM"
      sum
    } else {
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
      Last.time <- cnewdata[, c(bvar[length(bvar)], Cvar[1])]
      CompetingRisk <- object$CompetingRisk
      for (j in 1:N.ID) {
        subNDy.mean <- ynewdata.mean[ynewdata.mean[, bvar[length(bvar)]] == ID[j], ]
        subNDy.variance <- ynewdata.variance[ynewdata.variance[, bvar[length(bvar)]] == ID[j], ]
        subNDc <- cnewdata[cnewdata[, bvar[length(bvar)]] == ID[j], ]
        y.obs[[j]] <- data.frame(subNDy.mean[, c(bvar[1], Yvar[1])])
        
        allPi <- matrix(0, ncol = length(u), nrow = M)
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
        
        for (jj in 1:lengthu) {
          ## calculate the CIF
          CIF1 <- CIF1.CR(data, H01, H02, s, u[jj], meanbw)
          P1us <- Pk.us(CIF1, data, meanbw)
          Predraw1[j, jj] <- P1us
          CIF2 <- CIF2.CR(data, H01, H02, s, u[jj], meanbw)
          P2us <- Pk.us(CIF2, data, meanbw)
          Predraw2[j, jj] <- P2us
        }
      }
      names(y.obs) <- ID
      
      
      Pred <- list()
      for (jj in 1:N.ID) {
        Pred[[jj]] <- data.frame(u, Predraw1[jj, ], Predraw2[jj, ])
        colnames(Pred[[jj]]) <- c("times", "CIF1", "CIF2")
      }
      names(Pred) <- ID
      sum <- list(Pred = Pred, Last.time = Last.time, y.obs = y.obs, simulate = simulate, quadpoint = quadpoint)
      class(sum) <- "survfitJMMLSM"
      sum
      
    }
  }
  
  
}