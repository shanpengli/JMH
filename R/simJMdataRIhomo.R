##' @export
##' 

simJMdataRIhomo <- function(seed = 100, N = 200, increment = 0.7, beta = c(5, 1.5, 2, 1, 2),
                        sigma2 = exp(0.5),
                        gamma1 = c(1, 0.5, 0.5),
                        gamma2 = c(-0.5, 0.5, 0.25),
                        alpha1 = 1,
                        alpha2 = -1,
                        lambda1 = 0.05,
                        lambda2 = 0.025,
                        CL = 5,
                        CU = 10,
                        sigb = 0.5) {
  
  set.seed(seed)
  
  bi <- rnorm(N, mean = 0, sd = sqrt(sigb))
  
  ##covariate
  Z1 <- sample(c(0, 1), N, replace = TRUE, prob = c(0.5, 0.5))
  Z2 <- runif(N, min = -1, max = 1)
  Z3 <- rnorm(N, mean = 1, sd = 2)
  Z <- cbind(Z1, Z2, Z3)
  
  
  #hazard rate of risk1 and risk2
  
  ## non-informative cencoring time 
  C <- runif(N, min = CL, max = CU)
  risk1 <- vector()
  risk2 <- vector()
  for (i in 1:N) {
    temp=lambda1*exp(Z[i, ] %*% gamma1 + alpha1*bi[i])
    risk1[i] <- rexp(1, temp)
    temp=lambda2*exp(Z[i, ] %*% gamma2 + alpha2*bi[i])
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
  colnames(cdata) <- c("ID", "survtime", "cmprsk", "var1", "var2", "var3")
  
  ##fixed effects in longitudinal mean portion
  YdataRaw <- NULL
  for (i in 1:N) {
    ni <- floor(cdata[i, 2]/increment)
    suby <- matrix(0, nrow = ni+1, ncol = 3)
    suby[, 1] <- i
    suby[1, 2] <- beta[1] + beta[2]*Z[i, 1] + beta[3]*Z[i, 2] + beta[4]*Z[i, 3] + bi[i] + rnorm(1, mean = 0, sd = sqrt(sigma2))
    suby[1, 3] <- 0
    if (ni==0) {
      colnames(suby) <- c("ID", "Y", "time")
    } else {
      for (j in 1:ni) {
        suby[j+1, 2] <- beta[1] + beta[2]*Z[i, 1] + beta[3]*Z[i, 2] + beta[4]*Z[i, 3] + beta[5]*j + 
          bi[i] + rnorm(1, mean = 0, sd = sqrt(sigma2)) 
        suby[j+1, 3] <- j
      }
    }
    YdataRaw <- rbind(YdataRaw, suby)
  }
  Z <- cbind(ID, Z)
  Z <- as.data.frame(Z)
  YdataRaw <- as.data.frame(YdataRaw)
  colnames(YdataRaw) <- c("ID", "Y", "time")
  ydata <- dplyr::left_join(YdataRaw, Z, by = "ID")
  cdata <- as.data.frame(cdata)
  a <- list(cdata, ydata, table)
  names(a) <- c("cdata", "ydata", "rate")
  return(a)
}