GetinitFake <- function(cdata, ydata, long.formula, surv.formula, variance.formula,
                    model, ID, RE)
{
  
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  variance.formula <- all.vars(variance.formula)
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
  if (prod(variance.formula %in% ynames) == 0) {
    Fakename <- which(variance.formula %in% ynames == FALSE)
    stop(paste0("The WS variables ", long[Fakename], " not found"))
  }
  yID <- unique(as.character(ydata[, ID]))
  cID <- as.character(cdata[, ID])
  if (prod(yID == cID) == 0) {
    stop("The order of subjects in ydata doesn't match with cdata.")
  }
  
  ydim = dim(ydata)
  cdim = dim(cdata)
  
  orderdata <- sortdata(cdata, ydata, ID, surv.formula, long.formula)
  
  ydata <- orderdata$ydata
  cdata <- orderdata$cdata
  mdata <- orderdata$mdata
  
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
    Z <- as.matrix(Z)
  } else {
    stop("model should be one of the following options: interslope or intercept.")
  }
  
  W <- ydata[, variance.formula]
  W <- as.matrix(cbind(1, W))
  Y <- as.vector(ydata[, long[1]])
  X2 <- as.matrix(cdata[, survival[3:length(survival)]])
  survtime <- as.vector(cdata[, survival[1]])
  cmprsk <- as.vector(cdata[, survival[2]])
  
  if (sum(unique(cmprsk)) <= 3) {
    if (prod(c(0, 1, 2) %in% unique(cmprsk))) {
      beta = c(5, 1.5, 2, 1, 2)
      tau = c(0.5, 0.5, -0.2, 0.2, 0.05)
      gamma1 = c(1, 0.5, 0.5)
      gamma2 = c(-0.5, 0.5, 0.25)
      vee1 = 2
      vee2 = -2

      if (model == "intercept") {
        alpha1 = as.vector(-0.5)
        alpha2 = as.vector(0.5)
        Sig <- matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2, ncol = 2)
      } else {
        alpha1 = c(1, 0.7)
        alpha2 = c(-1, -0.5)
        Sig <- matrix(c(0.5, 0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.5), nrow = 3, ncol = 3)
      }
      
      a <- list(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, Sig, 
                Z, X, W, Y, X2, survtime, cmprsk, ydata, cdata, mdata)
      
      names(a) <- c("beta", "tau", "gamma1", "gamma2", "alpha1", "alpha2", "vee1", "vee2", "Sig",
                    "Z", "X1", "W", "Y", "X2", "survtime", "cmprsk", "ydata",
                    "cdata", "mdata")
      
      return(a)
      
    } 
    
    if (prod(c(0, 1) %in% unique(cmprsk))) {
      beta = c(5, 1.5, 2, 1, 2)
      tau = c(0.5, 0.5, -0.2, 0.2, 0.05)
      gamma1 = c(1, 0.5, 0.5)
      vee1 = 0.5
      
      if (model == "intercept") {
        alpha1 = as.vector(1)
        Sig <- matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2, ncol = 2)
      } else {
        alpha1 = c(1, 0.7)
        Sig <- matrix(c(0.5, 0.25, 0.25, 0.25, 0.5, 0.25, 0.25, 0.25, 0.5), nrow = 3, ncol = 3)
      }
      
      a <- list(beta, tau, gamma1, alpha1, vee1, Sig, 
                Z, X, W, Y, X2, survtime, cmprsk, ydata, cdata, mdata)
      
      names(a) <- c("beta", "tau", "gamma1", "alpha1", "vee1", "Sig",
                    "Z", "X1", "W", "Y", "X2", "survtime", "cmprsk", "ydata", 
                    "cdata", "mdata")
      
      return(a)
      
    }
  } else {
    stop(paste0("The ", survival[2], " variable is specified incorrectly! Program stops."))
  }
  
}