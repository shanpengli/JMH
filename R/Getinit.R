
Getinit <- function(cdata, ydata, long.formula, surv.formula, variance.var,
                    model, ID, RE) {
  
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  cnames <- colnames(cdata)
  ynames <- colnames(ydata)
  
  ##variable check
  if (prod(long %in% ynames) == 0) {
    Fakename <- which(long %in% ynames == FALSE)
    stop(paste0("The variable ", long[Fakename], " not found"))
  }
  if (prod(survival %in% cnames) == 0) {
    Fakename <- which(survival %in% cnames == FALSE)
    stop(paste0("The variable ", survival[Fakename], " not found"))
  }
  if (!(ID %in% ynames)) {
    stop(paste0("ID column ", ID, " not found in long_data"))
  }
  if (!(ID %in% cnames)) {
    stop(paste0("ID column ", ID, " not found in surv_data"))
  }
  if (is.null(RE) & model == "interslope") {
    stop("Random effects covariates must be specified.")
  }
  if (prod(variance.var %in% ynames) == 0) {
    Fakename <- which(variance.var %in% ynames == FALSE)
    stop(paste0("The WS variables ", long[Fakename], " not found"))
  }
  yID <- unique(ydata[, ID])
  cID <- cdata[, ID]
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cdata.")
  }
  
  ydim = dim(ydata)
  cdim = dim(cdata)
  
  ## extract covariates
  X <- ydata[, long[-1]]
  X <- as.matrix(cbind(1, X))
  
  ##random effect covariates
  if (model == "interslope") {
    if (prod(RE %in% ynames) == 0) {
      Fakename <- which(RE %in% ynames == FALSE)
      stop(paste0("The variable ", RE[Fakename], " not found"))
    } else {
      p1a <- 1 + length(RE)
      Z <- ydata[, RE]
      Z <- cbind(1, Z)
    }
  } else if (model == "intercept") {
    if (!is.null(RE)) {
      stop("You are fitting a mixed effects location scale model with random intercept only
           but random effects covariates are specified at the same time. Please respecify your model!")
    }
    p1a <- 1
    Z <- rep(1, ydim[1])
    Z <- as.data.frame(Z)
  } else {
    stop("model should be one of the following options: interslope or intercept.")
  }
  
  W <- ydata[, variance.var]
  W <- as.matrix(cbind(1, W))
  
  Y <- as.vector(ydata[, long[1]])
  
  ## get initial estimates of fixed effects in mixed effect location scale model
  BetaResid <- OLS(X, Y)
  ## get initial estimates of fixed effects in WS variance
  
  logResidsquare <- as.vector(log(BetaResid$resid^2))
  Tau <- OLS(W, logResidsquare)
  
  ## get initial estimates of fixed effects in competing risks model
  surv_xnam <- paste(survival[3:length(survival)], sep = "")
  survfmla.fixed <- paste(surv_xnam, collapse= "+")
  survfmla.out1 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==1)")
  survfmla <- as.formula(paste(survfmla.out1, survfmla.fixed, sep = "~"))
  fitSURV1 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
  gamma1 <- as.vector(fitSURV1$coefficients)
  
  survfmla.out2 <- paste0("survival::Surv(", survival[1], ", ", survival[2], "==2)")
  survfmla <- as.formula(paste(survfmla.out2, survfmla.fixed, sep = "~"))
  fitSURV2 <- survival::coxph(formula = survfmla, data = cdata, x = TRUE)
  gamma2 <- as.vector(fitSURV2$coefficients)
  
  X2 <- as.matrix(cdata[, survival[3:length(survival)]])
  
  survtime <- as.vector(cdata[, survival[1]])
  
  cmprsk <- as.vector(cdata[, survival[2]])
  
  a <- list(BetaResid, Tau, gamma1, gamma2, Z, X, W, Y, X2, survtime, cmprsk)
  
  names(a) <- c("Beta", "Tau", "gamma1", "gamma2", "Z", "X1", "W", "Y", "X2", "survtime", "cmprsk")
  
  return(a)
}