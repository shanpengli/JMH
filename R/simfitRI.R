##' @export

simfitRI <- function(sim = 100, N = 200, lambda1 = 0.05, lambda2 = 0.1,
                     beta = c(5, 1.5, 2, 1, 2),
                     tau = c(0.5, 0.5, -0.2, 0.2, 0.05),
                     gamma1 = c(1, 0.5, 0.5),
                     gamma2 = c(-0.5, 0.5, 0.25),
                     alpha1 = 1,
                     alpha2 = -1,
                     vee1 = 0.5,
                     vee2 = -0.5,
                     covbw = matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2, ncol = 2),
                     CL = 4, CU = 8, seed = 10, maxiter = 1000,
                     increment = 0.25, 
                     quadpoint = 15, 
                     ncores = 10) {
  
  ParaMatrixRaw <- parallel::mclapply(1:sim, bootsfitRI,
                                      N = N, lambda1 = lambda1, lambda2 = lambda2,
                                      beta = beta,
                                      tau = tau,
                                      gamma1 = gamma1,
                                      gamma2 = gamma2,
                                      alpha1 = alpha1,
                                      alpha2 = alpha2,
                                      vee1 = vee1,
                                      vee2 = vee2,
                                      CL = CL, CU = CU, 
                                      covbw = covbw, seed = seed, maxiter = maxiter,
                                      increment = increment, quadpoint = quadpoint,
                                      mc.cores = ncores)
  
  paramatrix <- as.data.frame(matrix(0, nrow = sim, ncol = 25))
  paramatrixSE <- as.data.frame(matrix(0, nrow = sim, ncol = 23))
  
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

    colnames(paramatrix)[count] <- paste0("alpha1_", 1)
    count <- count + 1
  
    colnames(paramatrix)[count] <- paste0("alpha2_", 1)
    count <- count + 1
  
  
  colnames(paramatrix)[count] <- paste0("vee1")
  count <- count + 1
  
  colnames(paramatrix)[count] <- paste0("vee2")
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
  
  name <- colnames(paramatrix)[-(24:25)]
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