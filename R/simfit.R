simfit <- function(sim = 100, N = 200, lambda1 = 0.05, lambda2 = 0.1,
                   tau = c(0.5, 0.5, -0.2, 0.2, 0.05),
                   CL = 4, CU = 8, seed = 10, maxiter = 1000,
                   increment = 0.25, quadpoint = 15, ncores = 10) {
  
  ParaMatrixRaw <- parallel::mclapply(1:sim, bootsfit,
                                                N = N, lambda1 = lambda1, lambda2 = lambda2,
                                                tau = tau,
                                                CL = CL, CU = CU, seed = seed, maxiter = maxiter,
                                                increment = increment, quadpoint = quadpoint,
                                                mc.cores = ncores)
  
  paramatrix <- as.data.frame(matrix(0, nrow = sim, ncol = 30))
  paramatrixSE <- as.data.frame(matrix(0, nrow = sim, ncol = 28))
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
  
  for (i in 1:3) {
    colnames(paramatrix)[count] <- paste0("gamma2_", i)
    count <- count + 1
  }
  
  for (i in 1:2) {
    colnames(paramatrix)[count] <- paste0("alpha1_", i)
    count <- count + 1
  }
  
  for (i in 1:2) {
    colnames(paramatrix)[count] <- paste0("alpha2_", i)
    count <- count + 1
  }
  
  colnames(paramatrix)[count] <- paste0("vee1")
  count <- count + 1
  
  colnames(paramatrix)[count] <- paste0("vee2")
  count <- count + 1
  
  for (i in 1:3) {
    colnames(paramatrix)[count] <- paste0("Sig_", i, i)
    count <- count + 1
  }
  
  colnames(paramatrix)[count] <- paste0("Sig_12")
  count <- count + 1
  colnames(paramatrix)[count] <- paste0("Sig_23")
  count <- count + 1
  colnames(paramatrix)[count] <- paste0("Sig_13")
  count <- count + 1
  colnames(paramatrix)[count] <- paste0("Time")
  count <- count + 1
  colnames(paramatrix)[count] <- paste0("Iter")
  
  name <- colnames(paramatrix)[-(29:30)]
  colnames(paramatrixSE) <- paste0("se", name)
  
  result <- list(paramatrix, paramatrixSE)
  names(result) <- c("paramatrix", "paramatrixSE")
  
  return(result)
  
}