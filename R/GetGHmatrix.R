GetGHmatrix <- function(quadpoint = quadpoint, p1a = p1a) {
  
  gq_vals <- statmod::gauss.quad(n = quadpoint, kind = "hermite")
  xs <- gq_vals$nodes
  ws <- gq_vals$weights
  
  xsmatrix <- wsmatrix <- matrix(0, nrow = (p1a+1), ncol = quadpoint^(p1a+1))
  xsmatrix[p1a+1, ] <- rep(xs, quadpoint^p1a)
  Total <- NULL
  for (i in 1:p1a) {
    for (j in 1:quadpoint) {
      sub <- rep(xs[j], quadpoint^i)
      Total <- c(Total, sub)
    }
    xsmatrix[p1a+1-i, ] <- rep(Total, quadpoint^(p1a-i))
    Total <- NULL
  }
  xsmatrix <- t(xsmatrix)
  
  wsmatrix[p1a+1, ] <- rep(ws, quadpoint^p1a)
  Total <- NULL
  for (i in 1:p1a) {
    for (j in 1:quadpoint) {
      sub <- rep(ws[j], quadpoint^i)
      Total <- c(Total, sub)
    }
    wsmatrix[p1a+1-i, ] <- rep(Total, quadpoint^(p1a-i))
    Total <- NULL
  }
  wsmatrix <- t(wsmatrix)
  
  return(list(xsmatrix = xsmatrix, wsmatrix = wsmatrix))
  
}