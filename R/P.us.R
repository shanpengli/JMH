P.us <- function(data, CH0u, bwl) {
  p1a <- nrow(data$Sig) - 1
  bl <- as.vector(bwl[1:p1a])
  wl <- bwl[p1a+1]
  (exp(-data$CH0*exp(data$X2%*%data$gamma + data$alpha%*%bl + data$nu*wl)) - 
     exp(-CH0u*exp(data$X2%*%data$gamma + data$alpha%*%bl + data$nu*wl)))/
    exp(-data$CH0*exp(data$X2%*%data$gamma + data$alpha%*%bl + data$nu*wl))
}