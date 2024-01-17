##' @export

simfitMESA <- function(sim = 100, N = 200, beta = c(80, -3, -3, 0.5, 5),
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
                       quadpoint = 6,
                     ncores = 10) {
  
  ParaMatrixRaw <- parallel::mclapply(1:sim, bootsfitMESA,
                                      N = N, lambda1 = lambda1, lambda2 = lambda2,
                                      beta = beta,
                                      tau = tau,
                                      gamma1 = gamma1,
                                      gamma2 = gamma2,
                                      alpha1 = alpha1,
                                      alpha2 = alpha2,
                                      vee1 = vee1,
                                      vee2 = vee2,
                                      lambdaC = lambdaC,
                                      incremin = incremin,
                                      incremax = incremax,
                                      Cmax = Cmax,
                                      covbw = covbw, seed = seed, maxiter = maxiter,
                                      quadpoint = quadpoint,
                                      mc.cores = ncores)
  
  paramatrix <- as.data.frame(matrix(0, nrow = sim, ncol = 34))
  paramatrixSE <- as.data.frame(matrix(0, nrow = sim, ncol = 32))
  
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
  
  for (i in 1:2) {
    colnames(paramatrix)[count] <- paste0("gamma1_", i)
    count <- count + 1
  }
  
  for (i in 1:2) {
    colnames(paramatrix)[count] <- paste0("gamma2_", i)
    count <- count + 1
  }
  
  for (i in 1:3) {
    colnames(paramatrix)[count] <- paste0("alpha1_", i)
    count <- count + 1
  }
  
  for (i in 1:3) {
    colnames(paramatrix)[count] <- paste0("alpha2_", i)
    count <- count + 1
  }
  
  colnames(paramatrix)[count] <- paste0("vee1")
  count <- count + 1
  
  colnames(paramatrix)[count] <- paste0("vee2")
  count <- count + 1
  
  for (i in 1:4) {
    colnames(paramatrix)[count] <- paste0("Sig_", i, i)
    count <- count + 1
  }
  
  for (i in 1:(4-1)) {
    for (j in 0:(4-i-1)) {
      colnames(paramatrix)[count] <- paste0("Sig_", j+1, i+j+1)
      count <- count + 1
    }
  }
  
  colnames(paramatrix)[count] <- paste0("Time")
  count <- count + 1
  colnames(paramatrix)[count] <- paste0("Iter")
  
  name <- colnames(paramatrix)[-(33:34)]
  colnames(paramatrixSE) <- paste0("se", name)
  
  ### failure rate
  table <- matrix(NA, nrow = sim, ncol = 3)
  for (i in 1:sim) {
    table[i, 1] <- ParaMatrixRaw[[i]]$rate[1, 2]
    table[i, 2] <- ParaMatrixRaw[[i]]$rate[2, 2]
    table[i, 3] <- ParaMatrixRaw[[i]]$rate[3, 2]
  }
  table <- as.data.frame(table)
  colnames(table) <- c("censoring", "risk1", "risk2")
  
  result <- list(paramatrix, paramatrixSE, table)
  names(result) <- c("paramatrix", "paramatrixSE", "table")
  
  return(result)
  
}