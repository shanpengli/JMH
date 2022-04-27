logLik <- function(data, bw) {
  p1a <- nrow(data$Sig) - 1
  b <- as.vector(bw[1:p1a])
  w <- bw[p1a+1]
  sigma <- exp(data$W%*%data$tau + w)
  sum((data$Y - data$X%*%data$beta - data$Z%*%b)^2/(2*sigma) + 0.5*log(sigma)) + data$CH0*exp(data$X2%*%data$gamma + data$alpha%*%b + data$nu*w) +
    t(bw)%*%solve(data$Sig)%*%bw/2
}