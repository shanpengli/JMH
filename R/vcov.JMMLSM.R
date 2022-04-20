##' @title Variance-covariance matrix of the estimated parameters for joint models
##' @name vcov
##' @aliases vcov.JMMLSM
##' @description Extract variance-covariance matrix for joint models.
##' @param object an object inheriting from class \code{JMMLSM}.
##' @param ... further arguments passed to or from other methods.
##' @return a matrix of variance covariance of all parameter estimates.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}}
##' @export
##' 
vcov.JMMLSM <- function(object, ...) {
  if (!inherits(object, "JMMLSM"))
    stop("Use only with 'JMMLSM' objects.\n")
  
  variance.formula <- as.formula(paste("", object$LongitudinalSubmodelvariance[3], sep = "~"))
  
  getdum <- getdummy(long.formula = object$LongitudinalSubmodelmean,
                     surv.formula = object$SurvivalSubmodel, 
                     variance.formula = variance.formula,
                     random = object$random, ydata = object$ydata, cdata = object$cdata)
  
  random <- all.vars(object$random)
  vcov <- as.data.frame(object$vcov)
  long <- names(object$beta)
  tau <- names(object$tau)
  survival <- names(object$gamma1)
  survival1 <- paste0("T.", survival)
  long <- paste0("Ymean.", long)
  tau <- paste0("Yvariance.", tau)
  p1a <- length(object$alpha1) + 1
  alpha1 <- rep("T.asso:", p1a-1)
  nu1 <- "T.asso.Residual_1:"
  sig <- "Sig"
  sig <- rep(sig, p1a*(p1a+1)/2)
  for (i in 1:p1a) sig[i] <- paste0(sig[i], i, i)
  if (p1a == 2) {
    sig[p1a+1] <- paste0(sig[p1a+1], "12")
    alpha1 <- paste0(alpha1, "(Intercept)_1")
  } else {
    sig[p1a+1] <- paste0(sig[p1a+1], "12")
    sig[p1a+2] <- paste0(sig[p1a+2], "23")
    sig[p1a+3] <- paste0(sig[p1a+3], "13")
    alpha1[1] <- paste0(alpha1[1], "(Intercept)_1")
    alpha1[2] <- paste0(alpha1[2], random[1], "_1")
  }
  
  if (object$CompetingRisk) {
    survival <- names(object$gamma2)
    survival2 <- paste0("T.", survival)
    alpha2 <- rep("T.asso:", p1a-1)
    nu2 <- "T.asso.Residual_2:"
    if (p1a == 2) {
      alpha2 <- paste0(alpha2, "(Intercept)_2")
    } else {
      alpha2[1] <- paste0(alpha2[1], "(Intercept)_2")
      alpha2[2] <- paste0(alpha2[2], random[1], "_2")
    }
    
    colnames(vcov) <- c(long, tau, survival1, survival2, alpha1, alpha2, nu1, nu2, sig)
    rownames(vcov) <- c(long, tau, survival1, survival2, alpha1, alpha2, nu1, nu2, sig)
  } else {
    colnames(vcov) <- c(long, tau, survival1, alpha1, nu1, sig)
    rownames(vcov) <- c(long, tau, survival1, alpha1, nu1, sig)
  }
  
  vcov
}