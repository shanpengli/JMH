##' @export
##'

bootsfitRI <- function(i, N = 200, lambda1 = 0.05, lambda2 = 0.1,
                     tau = c(0.5, 0.5, -0.2, 0.2, 0.05),
                     CL = 4, CU = 8, 
                     alpha1 = -0.5, alpha2 = 0.5, 
                     vee1 = 2, vee2 = -2,
                     covbw = matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2, ncol = 2),
                     seed = 10, maxiter = 1000,
                     increment = 0.25, quadpoint = 15) {
  
  data <- simJMdataRI(N = N, lambda1 = lambda1, lambda2 = lambda2,
                     tau = tau, alpha1 = alpha1, alpha2 = alpha2,
                     vee1 = vee1, vee2 = vee2, covbw = covbw,
                     CL = CL, CU = CU, seed = seed + i,
                     increment = increment)
  ydata <- data$ydata
  cdata <- data$cdata
  
  a <- proc.time()
  fit <- JMMLSM(cdata = cdata, ydata = ydata, 
                long.formula = Y ~ Z1 + Z2 + Z3 + time,
                surv.formula = Surv(survtime, cmprsk) ~ var1 + var2 + var3,
                variance.formula = ~ Z1 + Z2 + Z3 + time, maxiter = maxiter, epsilon = 1e-04, 
                quadpoint = quadpoint, random = ~ 1|ID)
  b <- proc.time()
  time <- (b - a)[3]
  
  coef <- vector()
  coefSE <- vector()
  count <- 1
  if (is.null(fit$beta)) {
    coef <- rep(NA, 25)
    coefSE <- rep(NA, 23)
    coef <- list(coef, coefSE)
    names(coef) <- c("coef", "coefSE")
    return(coef)
  } else if (fit$iter == maxiter) {
    coef <- rep(NA, 25)
    coefSE <- rep(NA, 23)
    coef <- list(coef, coefSE)
    names(coef) <- c("coef", "coefSE")
    return(coef)
  } else {
    for (j in 1:length(fit$beta)) {
      coef[count] <- fit$beta[j]
      coefSE[count] <- fit$sebeta[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$tau)) {
      coef[count] <- fit$tau[j]
      coefSE[count] <- fit$setau[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$gamma1)) {
      coef[count] <- fit$gamma1[j]
      coefSE[count] <- fit$segamma1[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$gamma2)) {
      coef[count] <- fit$gamma2[j]
      coefSE[count] <- fit$segamma2[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$alpha1)) {
      coef[count] <- fit$alpha1[j]
      coefSE[count] <- fit$sealpha1[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$alpha2)) {
      coef[count] <- fit$alpha2[j]
      coefSE[count] <- fit$sealpha2[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$vee1)) {
      coef[count] <- fit$vee1[j]
      coefSE[count] <- fit$sevee1[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$vee2)) {
      coef[count] <- fit$vee2[j]
      coefSE[count] <- fit$sevee2[j]
      count <- count + 1
    }
    
    for (j in 1:nrow(fit$Sig)) {
      coef[count] <- fit$Sig[j, j]
      coefSE[count] <- fit$seSig[j, j]
      count <- count + 1
    }
    
    coef[count] <- fit$Sig[1, 2]
    coefSE[count] <- fit$seSig[1, 2]
    count <- count + 1
    coef[count] <- time
    count <- count + 1
    coef[count] <- fit$iter
    
    coef <- list(coef, coefSE)
    names(coef) <- c("coef", "coefSE")
    
    return(coef)
  }
  
  
  
  
}