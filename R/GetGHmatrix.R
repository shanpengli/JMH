GetGHmatrix <- function(quadpoint = quadpoint, Sigb = Sigb) {
  
  gq_vals <- statmod::gauss.quad(n = quadpoint, kind = "hermite")
  xs <- gq_vals$nodes
  ws <- gq_vals$weights
  
  if (nrow(Sigb) == 2) {
    xsmatrix <- matrix(0, nrow = 3, ncol = quadpoint^3)
    wsmatrix <- xsmatrix
    xsmatrix[3, ] <- rep(xs, quadpoint^2)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(xs[i], quadpoint)
      Total <- c(Total, sub)
    }
    xsmatrix[2, ] <- rep(Total, quadpoint)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(xs[i], quadpoint^2)
      Total <- c(Total, sub)
    }
    xsmatrix[1, ] <- Total
    xsmatrix <- t(xsmatrix)
    
    wsmatrix[3, ] <- rep(ws, quadpoint^2)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(ws[i], quadpoint)
      Total <- c(Total, sub)
    }
    wsmatrix[2, ] <- rep(Total, quadpoint)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(ws[i], quadpoint^2)
      Total <- c(Total, sub)
    }
    wsmatrix[1, ] <- Total
    wsmatrix <- t(wsmatrix)
  } else {
    xsmatrix <- matrix(0, nrow = 2, ncol = quadpoint^2)
    wsmatrix <- xsmatrix
    xsmatrix[2, ] <- rep(xs, quadpoint)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(xs[i], quadpoint)
      Total <- c(Total, sub)
    }
    xsmatrix[1, ] <- Total
    xsmatrix <- t(xsmatrix)
    
    wsmatrix[2, ] <- rep(ws, quadpoint)
    Total <- NULL
    for (i in 1:quadpoint) {
      sub <- rep(ws[i], quadpoint)
      Total <- c(Total, sub)
    }
    wsmatrix[1, ] <- Total
    wsmatrix <- t(wsmatrix)
  }
  
  return(list(xsmatrix = xsmatrix, wsmatrix = wsmatrix))
  
}