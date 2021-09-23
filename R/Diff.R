Diff <- function(beta, prebeta, tau, pretau, gamma1, pregamma1, gamma2, pregamma2,
                 alpha1, prealpha1, alpha2, prealpha2, vee1, prevee1, vee2, prevee2,
                 Sig, preSig, H01, preH01, H02, preH02, epsilon) {
  
  betaAbsdiff <- max(abs(beta - prebeta))
  tauAbsdiff <- max(abs(tau - pretau))
  gamma1Absdiff <- max(abs(gamma1 - pregamma1))
  gamma2Absdiff <- max(abs(gamma2 - pregamma2))
  alpha1Absdiff <- max(abs(alpha1 - prealpha1))
  alpha2Absdiff <- max(abs(alpha2 - prealpha2))
  vee1Absdiff <- abs(vee1 - prevee1)
  vee2Absdiff <- abs(vee2 - prevee2)
  SigAbsdiff <- max(abs(Sig - preSig))
  H01Absdiff <- max(abs(H01[, 3] - preH01[, 3]))
  H02Absdiff <- max(abs(H02[, 3] - preH02[, 3]))
  
  if ((betaAbsdiff > epsilon) || (tauAbsdiff > epsilon) || (gamma1Absdiff > epsilon)
      || (gamma2Absdiff > epsilon) || (alpha1Absdiff > epsilon) || (alpha2Absdiff > epsilon)
      || (vee1Absdiff > epsilon) || (vee2Absdiff > epsilon) || (SigAbsdiff  > 10*epsilon)
      || (H01Absdiff > 10*epsilon) || (H02Absdiff > 10*epsilon)) {
    return(1)
  } else {
    return(0)
  }
  
}