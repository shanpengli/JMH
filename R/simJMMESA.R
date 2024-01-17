##' @export
##' 

simJMMESA <- function(seed = 100, N = 200, beta = c(80, -3, -3, 0.5, 5),
                      tau = c(2, 1, 1, 0.05, 0.3),
                      gamma1 = c(0.1, -0.3),
                      gamma2 = c(0.1, -0.3),
                      alpha1 = c(0, -0.03, 0),
                      alpha2 = c(-0.02, 0, 0.02),
                      vee1 = 0.7,
                      vee2 = 0.8,
                      lambda1 = 1e-2,
                      lambda2 = 2e-2,
                      lambdaC = 6e-2,
                      Cmax = 18,
                      incremin = 0.1, 
                      incremax = 1,
                      covbw = matrix(c(200, -70, -70, 10, -70, 150, 70, -2, -70, 
                                        70, 150, -2, 10, -2, -2, 0.7), nrow = 4, ncol = 4)
) {
  
  set.seed(seed)
  
  bwi <- MASS::mvrnorm(n = N, rep(0, 4), covbw, tol = 1e-6, empirical = FALSE)
  
  ##covariate
  Z1 <- rnorm(N, mean = 60, sd = 10)
  Z2 <- sample(c(0, 1), N, replace = TRUE, prob = c(0.5, 0.5))
  Z <- cbind(Z1, Z2)
  #hazard rate of risk1 and risk2
  
  ## non-informative cencoring time 
  risk1 <- vector()
  risk2 <- vector()
  C <- vector()
  for (i in 1:N) {
    C[i] <- min(rexp(1, rate = lambdaC), Cmax)
    temp=lambda1*exp(Z[i, ] %*% gamma1 + alpha1%*% bwi[i, 1:3] + vee1*bwi[i, 4])
    risk1[i] <- rexp(1, temp)
    temp=lambda2*exp(Z[i, ] %*% gamma2 + alpha2%*% bwi[i, 1:3] + vee2*bwi[i, 4])
    risk2[i] <- rexp(1, temp)
  }
  survtimeraw <- cbind(risk1, risk2, C)
  
  cmprsk <- vector()
  survtime <- vector()
  for (i in 1:N) {
    if (min(survtimeraw[i, ]) == survtimeraw[i, 1]) {
      cmprsk[i] <- 1
      survtime[i] <- survtimeraw[i, 1]
    } else if (min(survtimeraw[i, ]) == survtimeraw[i, 2]) {
      cmprsk[i] <- 2
      survtime[i] <- survtimeraw[i, 2]
    } else {
      cmprsk[i] <- 0
      survtime[i] <- survtimeraw[i, 3]
    }
  }
  survtimeraw <- as.data.frame(cbind(survtimeraw, survtime, cmprsk))
  table <- as.data.frame(table(survtimeraw$cmprsk)/N*100)
  writeLines(paste0("The censoring rate is: ", table[1, 2], "%"))
  writeLines(paste0("The risk 1 rate is: ", table[2, 2], "%"))
  writeLines(paste0("The risk 2 rate is: ", table[3, 2], "%"))
  ID <- c(1:N)
  cdata <- cbind(ID, survtimeraw$survtime, survtimeraw$cmprsk, Z)
  colnames(cdata) <- c("ID", "survtime", "cmprsk", "var1", "var2")
  
  ##fixed effects in longitudinal mean portion
  YdataRaw <- NULL
  rawni <- sample(c(0:6), N, replace = TRUE, prob = c(0.05, 0.05, 0.10, 0.10, 0.10, 0.15, 0.45))
  rawincrement <- runif(N, min = incremin, max = incremax)
  for (i in 1:N) {

    ni <- min(floor(cdata[i, 2]/rawincrement[i]), rawni[i])
    suby <- matrix(0, nrow = ni+1, ncol = 5)
    suby[, 1] <- i
    sd <- sqrt(exp(tau[1] + tau[4]*Z[i, 1] + tau[5]*Z[i, 2] + bwi[i, 4]))
    suby[1, 2] <- beta[1] + beta[4]*Z[i, 1] + beta[5]*Z[i, 2] + bwi[i, 1] + rnorm(1, mean = 0, sd = sd)
    suby[1, 3] <- suby[1, 4] <- suby[1, 5] <- 0
    if (ni==0) {
      colnames(suby) <- c("ID", "Y", "time", "timens1", "timens2")
    } else {
      time <- seq(rawincrement[i], ni*rawincrement[i], by = rawincrement[i])
      timebs <- as.matrix(splines::ns(time, df = 2))
      for (j in 1:ni) {
        sd <- sqrt(exp(tau[1] + tau[4]*Z[i, 1] + tau[5]*Z[i, 2] + tau[2]*timebs[j, 1] + tau[3]*timebs[j, 2] + bwi[i, 4]))
        suby[j+1, 2] <- beta[1] + beta[4]*Z[i, 1] + beta[5]*Z[i, 2] + beta[2]*timebs[j, 1] + beta[3]*timebs[j, 2] + 
          bwi[i, 1] + bwi[i, 2]*timebs[j, 1] + bwi[i, 3]*timebs[j, 2] + rnorm(1, mean = 0, sd = sd) 
        suby[j+1, 3] <- j*rawincrement[i]
        suby[j+1, 4] <- timebs[j, 1]
        suby[j+1, 5] <- timebs[j, 2]
      }
    }
    YdataRaw <- rbind(YdataRaw, suby)
    
  }
  Z <- cbind(ID, Z)
  Z <- as.data.frame(Z)
  YdataRaw <- as.data.frame(YdataRaw)
  ydata <- dplyr::left_join(YdataRaw, Z, by = "ID")
  cdata <- as.data.frame(cdata)
  a <- list(cdata, ydata, table)
  names(a) <- c("cdata", "ydata", "rate")
  return(a)
}




