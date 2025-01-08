bootPredhomo <- function(i = 1, seed = 10, N = 3000, lambda1 = 0.05, lambda2 = 0.15,
                        beta = c(5, 1.5, 2, 1),
                        sigma2 = exp(0.5),
                        gamma1 = c(1, 0.5, 0.5),
                        gamma2 = c(-0.5, 0.5, 0.25),
                        alpha1 = c(0.05, 0.01, 0.02),
                        alpha2 = c(-0.05, 0.02, 0.1),
                        CL = 4, CU = 8, 
                        sigb = matrix(c(10, 3, 2, 3, 4, 2, 2, 2, 4), nrow = 3),
                        increment = 0.25,
                        landmark.time = 3, horizon.time = 4:7, n.cv = 4, 
                        metric = c("MAEQ", "Brier", "Cindex", "pseudoR2")) {
  
  
  data <- simJMdatahomo(N = N, lambda1 = lambda1, lambda2 = lambda2,
                        beta = beta,
                        sigma2 = sigma2,
                        gamma1 = gamma1,
                        gamma2 = gamma2,
                        alpha1 = alpha1,
                        alpha2 = alpha2,
                        CL = CL, CU = CU, 
                        sigb = sigb,
                        seed = seed, increment = increment)
  
  ydata <- data$ydata
  cdata <- data$cdata
  
  fit <- JMH::JMMLSM(cdata = cdata, ydata = ydata, 
                     long.formula = Y ~ Z1 + time,
                     surv.formula = Surv(survtime, cmprsk) ~ var1 + var2 + var3,
                     variance.formula = ~ Z1 + time, 
                     quadpoint = 6, random = ~ time|ID, epsilon = 1e-3, print.para = TRUE)
  
  fit.FastJM <- FastJM::jmcs(cdata = cdata, ydata = ydata, 
                             long.formula = Y ~ Z1 + timens1 + timens2,
                             surv.formula = Surv(survtime, cmprsk) ~ var1 + var2 + var3,
                             quadpoint = 6, random = ~ timens1 + timens2|ID, tol = 1e-3, print.para = TRUE)
  
  if (metric == "MAEQ") {
    Pred <- JMH::MAEQJMMLSM(seed = seed + i, object = fit, landmark.time = landmark.time, 
                            horizon.time = horizon.time,
                            obs.time = "time", method = "GH", n.cv = n.cv, initial.para = TRUE)
    
    Pred.FastJM <- FastJM::MAEQjmcs(seed = seed + i, object = fit.FastJM, landmark.time = landmark.time, 
                                    horizon.time = horizon.time,
                                    obs.time = "time", method = "GH", n.cv = n.cv, initial.para = TRUE)
  } else if (metric == "Brier") {
    Pred <- JMH::PEJMMLSM(seed = seed + i, object = fit, landmark.time = landmark.time, 
                            horizon.time = horizon.time,
                            obs.time = "time", method = "GH", n.cv = n.cv, initial.para = TRUE)
    
    Pred.FastJM <- FastJM::PEjmcs(seed = seed + i, object = fit.FastJM, landmark.time = landmark.time, 
                                    horizon.time = horizon.time,
                                    obs.time = "time", method = "GH", n.cv = n.cv, initial.para = TRUE)
  } else if (metric == "Cindex") {
    Pred <- JMH::AUCJMMLSM(seed = seed + i, object = fit, landmark.time = landmark.time, 
                          horizon.time = horizon.time, obs.time = "time", method = "GH", 
                          n.cv = n.cv, initial.para = TRUE, metric = "Cindex")
    
    Pred.FastJM <- FastJM::AUCjmcs(seed = seed + i, object = fit.FastJM, landmark.time = landmark.time, 
                                  horizon.time = horizon.time, obs.time = "time", 
                                  method = "GH", n.cv = n.cv, initial.para = TRUE, metric = "Cindex")
  } else if (metric == "pseudoR2") {
    
    Pred <- JMH::PAJMMLSM(seed = seed + i, object = fit, landmark.time = landmark.time, 
                           horizon.time = horizon.time, obs.time = "time", method = "GH", 
                           n.cv = n.cv, initial.para = TRUE)
    
    Pred.FastJM <- FastJM::PAjmcs(seed = seed + i, object = fit.FastJM, landmark.time = landmark.time, 
                                   horizon.time = horizon.time, obs.time = "time", 
                                   method = "GH", n.cv = n.cv, initial.para = TRUE)
    
    
  } else {
    Pred <- Pred.FastJM <- NULL
  }
  
  
  
  return(list(Pred = Pred, Pred.FastJM = Pred.FastJM))
  
}


