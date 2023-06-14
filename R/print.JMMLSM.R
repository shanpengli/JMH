##' Print contents of JMMLSM object.
##' @title Print JMMLSM
##' @param x Object of class 'JMMLSM'.
##' @param digits number of digits of decimal to be printed.
##' @param ... Further arguments passed to or from other methods.
##' @return a summary of data, joint model, log likelihood, and parameter estimates.
##' @seealso \code{\link{JMMLSM}}
##' @author Shanpeng Li
##' @export
##' 
print.JMMLSM <- function(x, digits = 4, ...) {
  if (!inherits(x, "JMMLSM"))
    stop("Use only with 'JMMLSM' objects.\n")
  
  cat("\nCall:\n", sprintf(format(paste(deparse(x$call, width.cutoff = 500), collapse = ""))), "\n\n")
  
  if (x$CompetingRisk) {
    cat("Data Summary:\n")
    cat("Number of observations:", nrow(x$ydata), "\n")
    cat("Number of groups:", nrow(x$cdata), "\n\n")
    cat("Proportion of competing risks: \n")
    for (i in 1:2) {
      cat("Risk", i, ":", round(x$PropEventType[i+1, 2]/nrow(x$cdata)*100, 2), "%\n")
    }
    cat("\nNumerical intergration:\n")
    cat("Method: Standard Guass-Hermite quadrature\n")
    cat("Number of quadrature points: ", x$quadpoint, "\n")
    cat("\nModel Type: joint modeling of longitudinal continuous and competing risks data with the presence of intra-individual variability", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: Mixed effects location scale model\n")
    cat("Event process: cause-specific Cox proportional hazard model with non-parametric baseline hazard\n\n")
    cat("Loglikelihood: ", x$loglike, "\n\n")
    cat("Fixed effects in mean of longitudinal submodel: ",
        sprintf(format(paste(deparse(x$LongitudinalSubmodelmean, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    
    dat <- data.frame(x$beta, x$sebeta, x$beta/x$sebeta, 2 * pnorm(-abs(x$beta/x$sebeta)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)

    cat("\nFixed effects in variance of longitudinal submodel: ",
        sprintf(format(paste(deparse(x$LongitudinalSubmodelvariance, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    
    dat <- data.frame(x$tau, x$setau, x$tau/x$setau, 2 * pnorm(-abs(x$tau/x$setau)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)

    cat("\nSurvival sub-model fixed effects: ",
        sprintf(format(paste(deparse(x$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    dat <- data.frame(x$gamma1, x$segamma1, x$gamma1/x$segamma1, 2 * pnorm(-abs(x$gamma1/x$segamma1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    dat <- data.frame(x$gamma2, x$segamma2, x$gamma2/x$segamma2, 2 * pnorm(-abs(x$gamma2/x$segamma2)))
    colnames(dat) <- NULL
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    
    cat("\nAssociation parameters:                 \n")
    dat <- data.frame(x$alpha1, x$sealpha1, x$alpha1/x$sealpha1, 2 * pnorm(-abs(x$alpha1/x$sealpha1)))
    subdat <- data.frame(x$alpha2, x$sealpha2, x$alpha2/x$sealpha2, 2 * pnorm(-abs(x$alpha2/x$sealpha2)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    colnames(subdat) <- c("Estimate", "SE", "Z value", "p-val")
    dat <- rbind(dat, subdat)
    random <- all.vars(x$random)
    if (length(x$alpha1) == 1) rownames(dat) <- c("(Intercept)_1", "(Intercept)_2")
    if (length(x$alpha1) == 2) rownames(dat) <- c("(Intercept)_1", paste0(random[1], "_1"), "(Intercept)_2", paste0(random[1], "_2"))

    datnu <- matrix(0, nrow = 2, ncol = 4)
    datnu[1, ] <- c(x$vee1, x$sevee1, x$vee1/x$sevee1, 2 * pnorm(-abs(x$vee1/x$sevee1)))
    datnu[2, ] <- c(x$vee2, x$sevee2, x$vee2/x$sevee2, 2 * pnorm(-abs(x$vee2/x$sevee2)))
    datnu <- as.data.frame(datnu)
    colnames(datnu) <- c("Estimate", "SE", "Z value", "p-val")
    rownames(datnu) <- c("Residual_1", "Residual_2")
    dat <- rbind(dat, datnu)
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    cat("  Formula:", format(as.formula(x$random)), "\n")
    
    if (nrow(x$Sig) == 2) {
      
      dat <- matrix(0, nrow = 3, ncol = 4)
      dat[1, ] <- c(x$Sig[1,1], x$seSig[1,1], x$Sig[1,1]/x$seSig[1,1], 2 * pnorm(-abs(x$Sig[1,1]/x$seSig[1,1])))
      dat[2, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
      dat[3, ] <- c(x$Sig[2,2], x$seSig[2,2], x$Sig[2,2]/x$seSig[2,2], 2 * pnorm(-abs(x$Sig[2,2]/x$seSig[2,2])))
      dat <- as.data.frame(dat)
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      rownames(dat) <- c("(Intercept)", "(Intercept):Residual", "Residual")
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
    
    } else {
      
      dat <- matrix(0, nrow = 6, ncol = 4)
      dat[1, ] <- c(x$Sig[1,1], x$seSig[1,1], x$Sig[1,1]/x$seSig[1,1], 2 * pnorm(-abs(x$Sig[1,1]/x$seSig[1,1])))
      dat[2, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
      dat[3, ] <- c(x$Sig[2,2], x$seSig[2,2], x$Sig[2,2]/x$seSig[2,2], 2 * pnorm(-abs(x$Sig[2,2]/x$seSig[2,2])))
      for (i in 1:nrow(x$Sig)) dat[3+i, ] <- c(x$Sig[i,3], x$seSig[i,3], x$Sig[i,3]/x$seSig[i,3], 
                                             2 * pnorm(-abs(x$Sig[i,3]/x$seSig[i,3])))
      dat <- as.data.frame(dat)
      slope <- all.vars(x$random)[1]
      interslope <- paste0("(Intercept):", slope)
      residslope <- paste0(slope, ":Residual")
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      rownames(dat) <- c("(Intercept)", interslope, slope, "(Intercept):Residual", residslope, "Residual")
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
      
      
    }

    
  } else {
    cat("Data Summary:\n")
    cat("Number of observations:", nrow(x$ydata), "\n")
    cat("Number of groups:", nrow(x$cdata), "\n\n")
    cat("Proportion of events:", round(x$PropEventType[2, 2]/nrow(x$cdata)*100, 2), "%\n")
    cat("\nNumerical intergration:\n")
    cat("Method: Standard Guass-Hermite quadrature\n")
    cat("Number of quadrature points: ", x$quadpoint, "\n")
    cat("\nModel Type: joint modeling of longitudinal continuous and survival data with the presence of intra-individual variability", "\n\n")
    cat("Model summary:\n")
    cat("Longitudinal process: Mixed effects location scale model\n")
    cat("Event process: Cox proportional hazard model with non-parametric baseline hazard\n\n")
    cat("Loglikelihood: ", x$loglike, "\n\n")
    cat("Fixed effects in mean of longitudinal submodel: ",
        sprintf(format(paste(deparse(x$LongitudinalSubmodelmean, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    
    dat <- data.frame(x$beta, x$sebeta, x$beta/x$sebeta, 2 * pnorm(-abs(x$beta/x$sebeta)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    cat("\nFixed effects in variance of longitudinal submodel: \n",
        sprintf(format(paste(deparse(x$LongitudinalSubmodelvariance, width.cutoff = 500), collapse=""))), "\n")
    
    dat <- data.frame(x$tau, x$setau, x$tau/x$setau, 2 * pnorm(-abs(x$tau/x$setau)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    cat("\nSurvival sub-model fixed effects: ",
        sprintf(format(paste(deparse(x$SurvivalSubmodel, width.cutoff = 500), collapse=""))), "\n")
    cat("\n")
    dat <- data.frame(x$gamma1, x$segamma1, x$gamma1/x$segamma1, 2 * pnorm(-abs(x$gamma1/x$segamma1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat$"p-val" <- sprintf(paste("%.", digits, "f", sep = ""), dat$"p-val")
    print(dat)
    
    cat("\n Association parameters:                 \n")
    dat <- data.frame(x$alpha1, x$sealpha1, x$alpha1/x$sealpha1, 2 * pnorm(-abs(x$alpha1/x$sealpha1)))
    colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
    if (length(x$alpha1) == 1) rownames(dat) <- c("(Intercept)_1")
    if (length(x$alpha1) == 2) rownames(dat) <- c("(Intercept)_1", paste0(random[1], "_1"))
    
    datnu <- matrix(0, nrow = 1, ncol = 4)
    datnu[1, ] <- c(x$vee1, x$sevee1, x$vee1/x$sevee1, 2 * pnorm(-abs(x$vee1/x$sevee1)))
    datnu <- as.data.frame(datnu)
    colnames(datnu) <- c("Estimate", "SE", "Z value", "p-val")
    rownames(datnu) <- c("Residual_1")
    dat <- rbind(dat, datnu)
    dat[, 1:3] <- round(dat[, 1:3], digits+1)
    dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
    print(dat)
    cat("\n")
    
    cat("\nRandom effects:                 \n")
    cat("  Formula:", format(as.formula(x$random)), "\n")
    
    
    if (nrow(x$Sig) == 2) {
      
      dat <- matrix(0, nrow = 3, ncol = 4)
      dat[1, ] <- c(x$Sig[1,1], x$seSig[1,1], x$Sig[1,1]/x$seSig[1,1], 2 * pnorm(-abs(x$Sig[1,1]/x$seSig[1,1])))
      dat[2, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
      dat[3, ] <- c(x$Sig[2,2], x$seSig[2,2], x$Sig[2,2]/x$seSig[2,2], 2 * pnorm(-abs(x$Sig[2,2]/x$seSig[2,2])))
      dat <- as.data.frame(dat)
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      rownames(dat) <- c("(Intercept)", "(Intercept):Residual", "Residual")
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
      
    } else {
      
      dat <- matrix(0, nrow = 6, ncol = 4)
      dat[1, ] <- c(x$Sig[1,1], x$seSig[1,1], x$Sig[1,1]/x$seSig[1,1], 2 * pnorm(-abs(x$Sig[1,1]/x$seSig[1,1])))
      dat[2, ] <- c(x$Sig[1,2], x$seSig[1,2], x$Sig[1,2]/x$seSig[1,2], 2 * pnorm(-abs(x$Sig[1,2]/x$seSig[1,2])))
      dat[3, ] <- c(x$Sig[2,2], x$seSig[2,2], x$Sig[2,2]/x$seSig[2,2], 2 * pnorm(-abs(x$Sig[2,2]/x$seSig[2,2])))
      for (i in 1:nrow(x$Sig)) dat[3+i, ] <- c(x$Sig[i,3], x$seSig[i,3], x$Sig[i,3]/x$seSig[i,3], 
                                               2 * pnorm(-abs(x$Sig[i,3]/x$seSig[i,3])))
      dat <- as.data.frame(dat)
      slope <- all.vars(x$random)[1]
      interslope <- paste0("(Intercept):", slope)
      residslope <- paste0(slope, ":Residual")
      colnames(dat) <- c("Estimate", "SE", "Z value", "p-val")
      rownames(dat) <- c("(Intercept)", interslope, slope, "(Intercept):Residual", residslope, "Residual")
      dat[, 1:3] <- round(dat[, 1:3], digits+1)
      dat[, 4] <- sprintf(paste("%.", digits, "f", sep = ""), dat[, 4])
      print(dat)
      
      
    }
    
  }
}
