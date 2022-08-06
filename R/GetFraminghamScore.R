##' @export
##' 
##' 
GetFraminghamScore <- function(object, cdata, ydata, survtime, increment = 0.25, status,
                               ytime, ID, LandmarkTime, u = NULL, M = 100, quadpoint = 15,
                               quant.width = 0.25) {

  
  cdata <- cdata[cdata[, survtime] > LandmarkTime, ]
  cID <- cdata[, ID]
  ydata <- ydata[ydata[, ID] %in% cID, ]
  cdata2 <- cdata
  cdata2[, survtime] <- LandmarkTime
  ydata2 <- ydata[ydata[, ytime]*increment <= LandmarkTime, ]
  yID <- unique(ydata2[, ID])
  cdata2 <- cdata2[cdata2[, ID] %in% yID, ]
  if (is.null(u)) {
    maxtime1 <- object$H01[nrow(object$H01), 1]
    maxtime2 <- object$H02[nrow(object$H02), 1]
    u <- c(ceiling(max(maxtime1, maxtime2)) - 1, ceiling(max(maxtime1, maxtime2)))
  } else {
    next
  }
  survfit <- survfit2JMMLSM(object = object, seed = 100, ynewdata = ydata2, cnewdata = cdata2, 
                            u = u, M = M, simulate = TRUE, quadpoint = quadpoint)
  cID <- cdata2[, ID]
  AllCIF1 <- list()
  AllCIF2 <- list()
  for (j in 1:length(u)) {
    CIF1 <- vector()
    CIF2 <- vector()
    for (i in 1:length(cID)) {
      CIF1[i] <- survfit$Pred[[i]]$`Cumulative incidence probabilities for type 1 failure`[j, 2]
      CIF2[i] <- survfit$Pred[[i]]$`Cumulative incidence probabilities for type 2 failure`[j, 2]
    }
    AllsubjectCIF <- data.frame(cID, CIF1, CIF2)
    ## group subjects based on CIF
    quant1 <- quantile(AllsubjectCIF$CIF1, probs = seq(0, 1, by = quant.width))
    EmpiricalCIF1 <- rep(NA, 1/quant.width)
    PredictedCIF1 <- rep(NA, 1/quant.width)
    for (i in 1:(1/quant.width)) {
      subquant <- AllsubjectCIF[AllsubjectCIF$CIF1 > quant1[i] &
                                  AllsubjectCIF$CIF1 <= quant1[i+1], 1:2]
      quantsubdata <- cdata[cdata[, ID] %in% subquant$cID, 2:3]
      
      quantsubCIF <- GetEmpiricalCIF(data = quantsubdata, 
                                     time = quantsubdata[, survtime],
                                     status = status)
      
      quantsubRisk1 <- quantsubCIF$H1
      ii <- 1
      while (ii <= nrow(quantsubRisk1)) {
        if (quantsubRisk1[ii, 1] > u[j]) {
          if (ii >= 2) {
            EmpiricalCIF1[i] <- quantsubRisk1[ii-1, 4]
          } else {
            EmpiricalCIF1[i] <- 0
          }
          break
        } else {
          ii <- ii + 1
        }
      }
      if (is.na(EmpiricalCIF1[i])) EmpiricalCIF1[i] <- quantsubRisk1[nrow(quantsubRisk1), 4]
      PredictedCIF1[i] <- mean(subquant$CIF1)
    }
    AllCIF1[[j]] <- data.frame(EmpiricalCIF1, PredictedCIF1)
    
    quant2 <- quantile(AllsubjectCIF$CIF2, probs = seq(0, 1, by = quant.width))
    EmpiricalCIF2 <- rep(NA, 1/quant.width)
    PredictedCIF2 <- rep(NA, 1/quant.width)
    for (i in 1:(1/quant.width)) {
      subquant <- AllsubjectCIF[AllsubjectCIF$CIF2 > quant2[i] &
                                  AllsubjectCIF$CIF2 <= quant2[i+1], c(1, 3)]
      quantsubdata <- cdata[cdata[, ID] %in% subquant$cID, 2:3]
      
      quantsubCIF <- GetEmpiricalCIF(data = quantsubdata, 
                                     time = quantsubdata[, survtime],
                                     status = status)
      
      quantsubRisk2 <- quantsubCIF$H2
      ii <- 1
      while (ii <= nrow(quantsubRisk2)) {
        if (quantsubRisk2[ii, 1] > u[j]) {
          if (ii >= 2) {
            EmpiricalCIF2[i] <- quantsubRisk2[ii-1, 4]
          } else {
            EmpiricalCIF2[i] <- 0
          }
          break
        } else {
          ii <- ii + 1
        }
      }
      if (is.na(EmpiricalCIF2[i])) EmpiricalCIF2[i] <- quantsubRisk2[nrow(quantsubRisk2), 4]
      PredictedCIF2[i] <- mean(subquant$CIF2)
    }
    AllCIF2[[j]] <- data.frame(EmpiricalCIF2, PredictedCIF2)
  }
  names(AllCIF1) <- names(AllCIF2) <- u
  
  return(list(AllCIF1 = AllCIF1, AllCIF2 = AllCIF2))
}