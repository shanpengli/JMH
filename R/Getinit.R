Getinit <- function(cdata, ydata, long.formula, surv.formula, variance.formula,
                    model, ID, RE, random, survinitial) {
  
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  variance <- all.vars(variance.formula)
  cnames <- colnames(cdata)
  ynames <- colnames(ydata)
  
  if (is.null(RE) & model == "interslope") {
    stop("Random effects covariates must be specified.")
  }
  yID <- unique(ydata[, ID])
  cID <- cdata[, ID]
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cdata.")
  }
  
  ydim = dim(ydata)
  cdim = dim(cdata)
  
  orderdata <- sortdata(cdata, ydata, ID, surv.formula, long.formula)
  
  ydata <- orderdata$ydata
  cdata <- orderdata$cdata
  mdata <- orderdata$mdata
  
  ##random effect covariates
  if (model == "interslope") {
    if (prod(RE %in% ynames) == 0) {
      Fakename <- which(RE %in% ynames == FALSE)
      stop(paste0("The variable ", RE[Fakename], " not found in the longitudinal dataset.\n"))
    } else if (prod(RE %in% long) == 0) {
      Fakename <- which(RE %in% long == FALSE)
      stop(paste0("The variable ", RE[Fakename], " not found in the long.formula argument. 
                  Please include this variable in the random argument.\n"))
    } else {
      p1a <- 1 + length(RE)
      Z <- ydata[, RE]
      Z <- cbind(1, Z)
      Z <- as.matrix(Z)
    }
  } else if (model == "intercept") {
    if (!is.null(RE)) {
      stop("You are fitting a mixed effects location scale model with random intercept only
           but random effects covariates are specified at the same time. Please respecify your model!")
    }
    p1a <- 1
    Z <- rep(1, ydim[1])
    Z <- as.data.frame(Z)
    Z <- as.matrix(Z)
    
  } else {
    stop("model should be one of the following options: interslope or intercept.")
  }

  getdum <- getdummy(long.formula = long.formula, surv.formula = surv.formula,
                     variance.formula = variance.formula,
                     random = random, ydata = ydata, cdata = cdata)
  
  ydata.mean <- getdum$ydata.mean
  ydata.variance <- getdum$ydata.variance
  cdata <- getdum$cdata

  ## extract covariates
  X <- ydata.mean[, -c(1:2)]
  X <- as.matrix(X)
  W <- ydata.variance[, -1]
  W <- as.matrix(W)
  
  Y <- as.vector(ydata.mean[, 2])
  X2 <- as.matrix(cdata[, -c(1:3)])
  survtime <- as.vector(cdata[, survival[1]])
  cmprsk <- as.vector(cdata[, survival[2]])

  ## get initial estimates of fixed effects in mixed effect location scale model
  longfit <- nlme::lme(fixed = long.formula, random = random, data = ydata, method = "REML", 
                       control = nlme::lmeControl(opt = "opt"))
  beta <- longfit$coefficients$fixed
  D <- as.matrix(nlme::getVarCov(longfit))
  resid <- Y - as.vector(fitted(longfit))
  ## get initial estimates of fixed effects in WS variance
  logResidsquare <- as.vector(log(resid^2))
  Tau <- OLS(W, logResidsquare)
  tau <- Tau$betahat
  names(tau) <- colnames(W)
  
  if (sum(unique(cmprsk)) <= 3) {
    if (prod(c(0, 1, 2) %in% unique(cmprsk))) {
        ## get initial estimates of fixed effects in competing risks model
        survfmla.fixed <- surv.formula[3]
        survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
        survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
        fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
        
        survfmla.out2 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==2)")
        survfmla <- as.formula(paste(survfmla.out2, survfmla.fixed, sep = "~"))
        fitSURV2 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
        if (survinitial) {
          gamma1 <- fitSURV1$coefficients
          gamma2 <- fitSURV2$coefficients
      } else {
        gamma1 = as.vector(rep(0, ncol(X2)))
        names(gamma1) <- names(fitSURV1$coefficients)
        gamma2 = as.vector(rep(0, ncol(X2)))
        names(gamma2) <- names(fitSURV1$coefficients)
      }
      
      
      vee1 = 0
      vee2 = 0
      
      if (model == "intercept") {
        alpha1 = as.vector(0)
        alpha2 = as.vector(0)
        
        Sig <- matrix(0, nrow = 2, ncol = 2)
        Sig[1, 1] <- D
        Sig[2, 2] <- 1
        Sig[1, 2] <- 0.5*D
        Sig[2, 1] <- Sig[1, 2]
        
      } else {
        alpha1 = c(0, 0)
        alpha2 = c(0, 0)
        Sig <- matrix(0, nrow = 3, ncol = 3)
        Sig[1:2, 1:2] <- D
        Sig[3, 3] <- 1
        Sig[1:2, 3] <- c(Sig[1, 1], Sig[2, 2])*0.5
        Sig[3, 1:2] <- Sig[1:2, 3] 
      }
      
      a <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, Sig, 
                Z, X, W, Y, X2, survtime, cmprsk, cdata, mdata)
      
      names(a) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1", "vee2", "Sig",
                    "Z", "X1", "W", "Y", "X2", "survtime", "cmprsk",
                    "cdata", "mdata")
      
      return(a)
      
    } 
    
    if (prod(c(0, 1) %in% unique(cmprsk))) {

      ## get initial estimates of fixed effects in survival model
      survfmla.fixed <- surv.formula[3]
      survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
      survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
      fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
      
      if (survinitial) {
        gamma1 <- fitSURV1$coefficients
      } else {
        gamma1 = as.vector(rep(0, ncol(X2)))
        names(gamma1) <- names(fitSURV1$coefficients)
      }
      
      vee1 = 0
      
      if (model == "intercept") {
        alpha1 = as.vector(0)
        Sig <- matrix(0, nrow = 2, ncol = 2)
        Sig[1, 1] <- D
        Sig[2, 2] <- 1
        Sig[1, 2] <- 0.5*D
        Sig[2, 1] <- Sig[1, 2]
      } else {
        alpha1 = c(0, 0)
        Sig <- matrix(0, nrow = 3, ncol = 3)
        Sig[1:2, 1:2] <- D
        Sig[3, 3] <- 1
        Sig[1:2, 3] <- c(Sig[1, 1], Sig[2, 2])*0.5
        Sig[3, 1:2] <- Sig[1:2, 3] 
      }
      
      a <- list(beta, tau, gamma1, alpha1, vee1, Sig, 
                Z, X, W, Y, X2, survtime, cmprsk, cdata, mdata)
      
      names(a) <- c("beta", "tau", "gamma1", "alpha1", "vee1", "Sig",
                    "Z", "X1", "W", "Y", "X2", "survtime", "cmprsk",
                    "cdata", "mdata")
      
      return(a)
      
    }
  } else {
    stop(paste0("The ", survival[2], " variable is specified incorrectly! Program stops."))
  }
  
}