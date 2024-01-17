##' @export

bootsfitMESA <- function(i, N = 200,
                         beta = c(80, -3, -3, 0.5, 5),
                         tau = c(2, 1, 1, 0.05, 0.3),
                         gamma1 = c(0.1, -0.3),
                         gamma2 = c(0.1, -0.3),
                         alpha1 = c(0, -0.03, 0),
                         alpha2 = c(-0.02, 0, 0.02),
                         vee1 = 0.7,
                         vee2 = 0.8,
                         covbw = matrix(c(200, -70, -70, 10, 
                                        -70, 150, 70, -2, 
                                        -70, 70, 150, -2, 
                                        10, -2, -2, 0.7), 
                                      nrow = 4, 
                                      ncol = 4),
                         lambda1 = 1e-5,
                         lambda2 = 3e-5, 
                         lambdaC = 1e-2, 
                         incremin = 1, 
                         incremax = 2.5,
                         Cmax = 18,
                         seed = 100, maxiter = 1000,
                         quadpoint = 6) {
  
  data <- simJMMESA(seed = seed + i, N = N, beta = beta,
                    tau = tau,
                    gamma1 = gamma1,
                    gamma2 = gamma2,
                    alpha1 = alpha1,
                    alpha2 = alpha2,
                    vee1 = vee1,
                    vee2 = vee2,
                    lambda1 = lambda1,
                    lambda2 = lambda2,
                    lambdaC = lambdaC,
                    incremin = incremin,
                    incremax = incremax,
                    Cmax = Cmax,
                    covbw = covbw)
  ydata <- data$ydata
  cdata <- data$cdata
  rate <- data$rate
  
  a <- proc.time()
  fit <- try(JMMLSM(cdata = cdata, ydata = ydata, 
                           long.formula = Y ~ timens1 + timens2 + Z1 + Z2,
                           surv.formula = Surv(survtime, cmprsk) ~ var1 + var2,
                           variance.formula = ~ timens1 + timens2 + Z1 + Z2, 
                           quadpoint = quadpoint, random = ~ timens1 + timens2|ID,
                           initial.para = list(beta = beta,
                                               tau = tau,
                                               gamma1 = gamma1,
                                               gamma2 = gamma2,
                                               alpha1 = alpha1,
                                               alpha2 = alpha2,
                                               vee1 = vee1,
                                               vee2 = vee2,
                                               Sig = covbw),
                    epsilon = 1e-3), silent = TRUE)
  b <- proc.time()
  time <- (b - a)[3]
  
  coef <- vector()
  coefSE <- vector()
  count <- 1
  if ('try-error' %in% class(fit)) {
    coef <- rep(NA, 34)
    coefSE <- rep(NA, 32)
    coef <- list(coef, coefSE, rate)
    names(coef) <- c("coef", "coefSE", "rate")
    return(coef)
  } else if (is.null(fit$beta)) {
    coef <- rep(NA, 34)
    coefSE <- rep(NA, 32)
    coef <- list(coef, coefSE, rate)
    names(coef) <- c("coef", "coefSE", "rate")
    return(coef)
  } else if (fit$iter == maxiter) {
    coef <- rep(NA, 34)
    coefSE <- rep(NA, 32)
    coef <- list(coef, coefSE, rate)
    names(coef) <- c("coef", "coefSE", "rate")
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
    
    for (i in 1:(nrow(fit$Sig)-1)) {
      for (j in 0:(nrow(fit$Sig)-i-1)) {
        coef[count] <- fit$Sig[j+1,i+j+1]
        coefSE[count] <- fit$seSig[j+1,i+j+1]
        count <- count + 1
      }
    }

    coef[count] <- time
    count <- count + 1
    coef[count] <- fit$iter
    
    coef <- list(coef, coefSE, rate)
    names(coef) <- c("coef", "coefSE", "rate")
    
    return(coef)
  }
  
}