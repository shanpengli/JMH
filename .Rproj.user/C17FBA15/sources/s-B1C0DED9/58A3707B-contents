##' @export
##'

bootsfit <- function(i, N = 200, lambda1 = 0.05, lambda2 = 0.1,
                     tau = c(0.5, 0.5, -0.2, 0.2, 0.05),
                     CL = 4, CU = 8, seed = 10, maxiter = 1000,
                     increment = 0.25, quadpoint = 15) {
  
  data <- simJMdata(N = N, lambda1 = lambda1, lambda2 = lambda2,
                    tau = tau,
                    CL = CL, CU = CU, seed = seed + i,
                    increment = increment)
  ydata <- data$ydata
  cdata <- data$cdata
  
  a <- proc.time()
  fit <- JMMLSM(cdata = cdata, ydata = ydata, 
                long.formula = Y ~ Z1 + Z2 + Z3 + time,
                surv.formula = Surv(survtime, cmprsk) ~ var1 + var2 + var3,
                variance.var = c("Z1", "Z2", "Z3", "time"), maxiter = maxiter, epsilon = 1e-04, 
                quadpoint = quadpoint, ID = "ID", RE = "time",
                model = "interslope")
  b <- proc.time()
  time <- (b - a)[3]
  
  coef <- vector()
  count <- 1
  if (is.null(fit$beta)) {
    coef <- rep(NA, 29)
    return(coef)
  } else if (fit$iter == maxiter) {
    coef <- rep(NA, 29)
    return(coef)
  } else {
    for (j in 1:length(fit$beta)) {
      coef[count] <- fit$beta[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$tau)) {
      coef[count] <- fit$tau[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$gamma1)) {
      coef[count] <- fit$gamma1[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$gamma2)) {
      coef[count] <- fit$gamma2[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$alpha1)) {
      coef[count] <- fit$alpha1[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$alpha2)) {
      coef[count] <- fit$alpha2[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$vee1)) {
      coef[count] <- fit$vee1[j]
      count <- count + 1
    }
    
    for (j in 1:length(fit$vee2)) {
      coef[count] <- fit$vee2[j]
      count <- count + 1
    }
    
    for (j in 1:nrow(fit$Sig)) {
      coef[count] <- fit$Sig[j, j]
      count <- count + 1
    }
    
    coef[count] <- fit$Sig[1, 2]
    count <- count + 1
    coef[count] <- fit$Sig[2, 3]
    count <- count + 1
    coef[count] <- fit$Sig[1, 3]
    count <- count + 1
    coef[count] <- time
    return(coef)
  }
  
  
  
  
}