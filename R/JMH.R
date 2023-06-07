#' @useDynLib JMH, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats  fitted model.frame model.matrix optim quantile rexp rnorm runif as.formula  pnorm  pchisq complete.cases
#' @importFrom statmod  gauss.quad
#' @importFrom utils  read.table
#' @importFrom survival coxph Surv survfit
#' @importFrom parallel mclapply parLapply makeCluster stopCluster
#' @importFrom dplyr left_join
#' @importFrom MASS mvrnorm
#' @importFrom nlme lme getVarCov lmeControl
#' @importFrom caret groupKFold
#' @importFrom graphics axis mtext par segments title
NULL