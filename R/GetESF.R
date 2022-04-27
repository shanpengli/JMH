GetESF <- function(beta, tau, gamma1, alpha1, vee1, H01, 
                 Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix) {
  
  
  n <- nrow(X2)
  p1a <- ncol(Z)
  
  CUH01 <- rep(0, n)
  HAZ01 <- rep(0, n)
  
  CumuH01 <- cumsum(H01[, 3])
  
  getHazardSF(CumuH01, survtime, cmprsk, H01, CUH01, HAZ01)
  
  
  status = getECSF(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, 
                   survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01)
  
  return(status)
  
}