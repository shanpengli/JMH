#' @useDynLib JMH, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats  as.formula  pnorm  pchisq complete.cases
#' @importFrom statmod  gauss.quad
#' @importFrom utils  read.table
#' @importFrom survival coxph
#' @importFrom parallel mclapply parLapply makeCluster stopCluster
#' @importFrom dplyr left_join
#' @importFrom MASS mvrnorm
#' @importFrom lme4 lmer fixef VarCorr 
NULL