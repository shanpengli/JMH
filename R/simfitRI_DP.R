##' @export
##'

simfitRI_DP <- function(sim = 10, Nt = 200, Nv = 2, seed = 10, increment = 0.7, beta = c(5, 1.5, 2, 1, 2),
                        tau = c(0.5, 0.5, -0.2, 0.2, 0.05),
                        gamma1 = c(1, 0.5, 0.5),
                        gamma2 = c(-0.5, 0.5, 0.25),
                        alpha1 = 1,
                        alpha2 = -1,
                        vee1 = 0.5,
                        vee2 = -0.5,
                        lambda1 = 0.05,
                        lambda2 = 0.1,
                        CL = 4,
                        CU = 8,
                        covbw = matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2, ncol = 2),
                        M = 100, maxiter = 1000, quadpoint = 10, 
                        u = seq(4.5, 6.5, by = 1), ncores = 10) {
  
  validated_set <- simJMdataRI_DP(seed = seed, N = Nv, increment = increment,
                                  beta = beta, tau = tau, gamma1 = gamma1,
                                  gamma2 = gamma2, alpha1 = alpha1,
                                  alpha2 = alpha2, vee1 = vee1, vee2 = vee2,
                                  lambda1 = lambda1, lambda2 = lambda2, 
                                  CL = CL, CU = CU, covbw = covbw)
  
  vali.ydata <- validated_set$ydata
  vali.cdata <- validated_set$cdata
  vali.lam1 <- validated_set$lam1
  vali.lam2 <- validated_set$lam2
  
  vali.cdata <- vali.cdata[vali.cdata$survtime>=u[1], ]
  vali.ydata <- vali.ydata[vali.ydata$ID %in% vali.cdata$ID, ]
  vali.ydata <- vali.ydata[vali.ydata$time <= u[1]/increment, ]
  
  vali.lam1 <- vali.lam1[vali.lam1$ID  %in% vali.cdata$ID, 2]
  vali.lam2 <- vali.lam2[vali.lam2$ID  %in% vali.cdata$ID, 2]
  
  vali.ydata.x <- as.matrix(cbind(1, vali.ydata[, c("Z1", "Z2", "Z3", "time")]))
  vali.ydata.xbeta <- vali.ydata.x %*% beta
  vali.ydata.wtau <- vali.ydata.x %*% tau
  vali.ydata.Y <- vali.ydata$Y
  
  getGH <- GetGHmatrix(quadpoint = quadpoint, Sigb = as.matrix(covbw[1, 1]))
  
  xsmatrix <- getGH$xsmatrix
  wsmatrix <- getGH$wsmatrix
  
  ni <- nrow(vali.ydata)/nrow(vali.cdata)
  n <- nrow(vali.cdata)
  mdata <- rep(ni, n)
  mdataS <- rep(0, n) 
  mdataS[1] <- 1
  mdataCum <- cumsum(mdata)
  mdata2 <- mdata - 1
  mdataS[2:n] <- mdataCum[2:n] - mdata2[2:n]
  vali.cdata$survtime <- rep(u[1], n)
  # return(list(vali.ydata = vali.ydata, vali.cdata = vali.cdata,
  #             vali.lam1 = vali.lam1, vali.lam2 = vali.lam2,
  #             vali.ydata.xbeta = vali.ydata.xbeta, vali.ydata.wtau = vali.ydata.wtau,
  #             vali.ydata.Y = vali.ydata.Y, xsmatrix = xsmatrix, wsmatrix = wsmatrix,
  #             ni = ni, n = n, mdata = mdata, mdataS = mdataS))
  
  ## calculate the true generated Pik
  TruePik <- getTruePik(vali.ydata.Y, vali.ydata.xbeta, vali.ydata.wtau, vali.lam1, vali.lam2,
                        alpha1, alpha2, vee1, vee2, covbw, ni, mdataS, xsmatrix, wsmatrix, u)

  TrueP1ik <- as.data.frame(cbind(vali.cdata$ID, TruePik$P1us))
  TrueP2ik <- as.data.frame(cbind(vali.cdata$ID, TruePik$P2us))
  colnames(TrueP1ik) <- c("ID", u[-1])
  colnames(TrueP2ik) <- c("ID", u[-1])

  # return(list(TrueP1ik = TrueP1ik, TrueP2ik = TrueP2ik))
  
  ParaMatrixRaw <- parallel::mclapply(1:sim, bootsfitRI_DP,
                                      seed = seed, N = Nt, increment = increment,
                                      beta = beta, tau = tau, gamma1 = gamma1,
                                      gamma2 = gamma2, alpha1 = alpha1, alpha2 = alpha2,
                                      vee1 = vee1, vee2 = vee2,
                                      lambda1 = lambda1, lambda2 = lambda2,
                                      CL = CL, CU = CU, covbw = covbw,
                                      quadpoint = quadpoint, maxiter = maxiter,
                                      u = u[-1], ynewdata = vali.ydata, cnewdata = vali.cdata,
                                      M = M, mc.cores = ncores)
  
  result <- list(TrueP1ik = TrueP1ik, TrueP2ik = TrueP2ik, ParaMatrixRaw = ParaMatrixRaw)

  return(result)
  
  
  
}