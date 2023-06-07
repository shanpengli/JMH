simfitRISF <- function(sim = 100, N = 200, lambda1 = 0.05,
                      tau = c(0.5, 0.5, -0.2, 0.2, 0.05),
                      CL = 4, CU = 8, seed = 10, maxiter = 1000,
                      increment = 0.25, 
                      alpha1 = 1,
                      vee1 = 0.5,
                      quadpoint = 15, 
                      covbw = matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2, ncol = 2),
                      ncores = 10) {
  
  ParaMatrixRaw <- parallel::mclapply(1:sim, bootsfitRISF,
                                      N = N, lambda1 = lambda1, 
                                      tau = tau, CL = CL, CU = CU, 
                                      alpha1 = alpha1, 
                                      vee1 = vee1, 
                                      covbw = covbw, seed = seed, maxiter = maxiter,
                                      increment = increment, quadpoint = quadpoint,
                                      mc.cores = ncores)
  
  paramatrix <- as.data.frame(matrix(0, nrow = sim, ncol = 20))
  paramatrixSE <- as.data.frame(matrix(0, nrow = sim, ncol = 18))
  
  for (i in 1:sim) {
    paramatrix[i, ] <- ParaMatrixRaw[[i]]$coef
    paramatrixSE[i, ] <- ParaMatrixRaw[[i]]$coefSE
  }
  
  count <- 1
  for (i in 1:5) {
    colnames(paramatrix)[count] <- paste0("beta_", i-1)
    count <- count + 1
  }
  
  for (i in 1:5) {
    colnames(paramatrix)[count] <- paste0("tau_", i-1)
    count <- count + 1
  }
  
  for (i in 1:3) {
    colnames(paramatrix)[count] <- paste0("gamma1_", i)
    count <- count + 1
  }
  
  colnames(paramatrix)[count] <- paste0("alpha1_", 1)
  count <- count + 1
  
  colnames(paramatrix)[count] <- paste0("vee1")
  count <- count + 1
  
  for (i in 1:2) {
    colnames(paramatrix)[count] <- paste0("Sig_", i, i)
    count <- count + 1
  }
  
  colnames(paramatrix)[count] <- paste0("Sig_12")
  count <- count + 1
  colnames(paramatrix)[count] <- paste0("Time")
  count <- count + 1
  colnames(paramatrix)[count] <- paste0("Iter")
  
  name <- colnames(paramatrix)[-(19:20)]
  colnames(paramatrixSE) <- paste0("se", name)
  
  result <- list(paramatrix, paramatrixSE)
  names(result) <- c("paramatrix", "paramatrixSE")
  
  return(result)
  
}