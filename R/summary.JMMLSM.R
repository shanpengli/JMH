##' @export
##' 

summary.JMMLSM <- function(object, process = c("longitudinal", "survival"), digits = 4) {
  
  if (!inherits(object, "JMMLSM"))
    stop("Use only with 'JMMLSM' objects.\n")
  
  if (process == "longitudinal") {
    ##Estimates of betas
    Estimate <- object$beta
    SE <- object$sebeta
    LowerLimit <- Estimate - 1.96 * SE
    UpperLimit <- Estimate + 1.96 * SE
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out <- data.frame(Estimate, SE, LowerLimit, UpperLimit, pval)
    out <- cbind(paste0("Mean_", rownames(out)), out)
    rownames(out) <- NULL
    colnames(out)[1] <- "Parameter"
    
    ##Estimates of tau
    Estimate <- object$tau
    SE <- object$setau
    LowerLimit <- Estimate - 1.96 * SE
    UpperLimit <- Estimate + 1.96 * SE
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out2 <- data.frame(Estimate, SE, LowerLimit, UpperLimit, pval)
    out2 <- cbind(paste0("Var_", rownames(out2)), out2)
    rownames(out2) <- NULL
    colnames(out2)[1] <- "Parameter"
    
    out3 <- rbind(out, out2)
    
    rownames(out3) <- NULL
    names(out3) <- c("Longitudinal", "coef", "SE", "95%Lower", "95%Upper", "p-values")
    
    out3[, 2:ncol(out3)] <- round(out3[, 2:ncol(out3)], digits = digits)
    out3[, ncol(out3)] <- format(out3[, ncol(out3)], scientific = FALSE)
    
    return(out3)
    
  } else if (process == "survival") {
    ##gamma
    Estimate <- object$gamma1
    SE <- object$segamma1
    LowerLimit <- Estimate - 1.96 * SE
    UpperLimit <- Estimate + 1.96 * SE
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, pval)
    out <- cbind(rownames(out), out)
    rownames(out) <- NULL
    colnames(out)[1] <- "Parameter"
    
    Estimate <- object$gamma2
    SE <- object$segamma2
    LowerLimit <- Estimate - 1.96 * SE
    UpperLimit <- Estimate + 1.96 * SE
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out2 <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, pval)
    out2 <- cbind(rownames(out2), out2)
    rownames(out2) <- NULL
    colnames(out2)[1] <- "Parameter"
    outgamma <- rbind(out, out2)
    names(outgamma) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", "p-values")
    
    ##alpha
    Estimate <- object$alpha1
    if (length(Estimate) == 2) names(Estimate) <- c("alpha1_1", "alpha1_2")
    if (length(Estimate) == 1) names(Estimate) <- c("alpha1_1")
    SE <- object$sealpha1
    LowerLimit <- Estimate - 1.96 * SE
    UpperLimit <- Estimate + 1.96 * SE
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, pval)
    out <- cbind(rownames(out), out)
    rownames(out) <- NULL
    colnames(out)[1] <- "Parameter"
    
    Estimate <- object$alpha2
    if (length(Estimate) == 2) names(Estimate) <- c("alpha2_1", "alpha2_2")
    if (length(Estimate) == 1) names(Estimate) <- c("alpha2_1")
    SE <- object$sealpha2
    LowerLimit <- Estimate - 1.96 * SE
    UpperLimit <- Estimate + 1.96 * SE
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out2 <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, pval)
    out2 <- cbind(rownames(out2), out2)
    rownames(out2) <- NULL
    colnames(out2)[1] <- "Parameter"
    outalpha <- rbind(out, out2)
    names(outalpha) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", "p-values")
    
    #vee
    Estimate <- c(object$vee1, object$vee2)
    names(Estimate) <- c("vee1", "vee2")
    SE <- c(object$sevee1, object$sevee2)
    LowerLimit <- Estimate - 1.96 * SE
    UpperLimit <- Estimate + 1.96 * SE
    zval = (Estimate/SE)
    pval = 2 * pnorm(-abs(zval))
    out <- data.frame(Estimate, exp(Estimate), SE, LowerLimit, UpperLimit, pval)
    out <- cbind(rownames(out), out)
    rownames(out) <- NULL
    colnames(out)[1] <- "Parameter"
    names(out) <- c("Survival", "coef", "exp(coef)", "SE(coef)", "95%Lower", "95%Upper", "p-values")
    
    out <- rbind(outgamma, outalpha, out)
    
    out[, 2:ncol(out)] <- round(out[, 2:ncol(out)], digits = digits)
    out[, ncol(out)] <- format(out[, ncol(out)], scientific = FALSE)
    
    return(out)
  }
  
  
}