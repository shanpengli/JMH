# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

GetrisksetC <- function(cdata) {
    .Call(`_JMH_GetrisksetC`, cdata)
}

GetrisksetCSF <- function(cdata) {
    .Call(`_JMH_GetrisksetCSF`, cdata)
}

OLS <- function(X, Y) {
    .Call(`_JMH_OLS`, X, Y)
}

getdeterminant <- function(H) {
    .Call(`_JMH_getdeterminant`, H)
}

MultVV <- function(x, y) {
    .Call(`_JMH_MultVV`, x, y)
}

MultVVoutprod <- function(x) {
    .Call(`_JMH_MultVVoutprod`, x)
}

MultVV2outprod <- function(x, y) {
    .Call(`_JMH_MultVV2outprod`, x, y)
}

MultVVinprod <- function(x) {
    .Call(`_JMH_MultVVinprod`, x)
}

CH <- function(H, t) {
    .Call(`_JMH_CH`, H, t)
}

HAZ <- function(H, t) {
    .Call(`_JMH_HAZ`, H, t)
}

MultMM <- function(x, y) {
    .Call(`_JMH_MultMM`, x, y)
}

GetCIF1CR <- function(gamma1, gamma2, alpha1, alpha2, nu1, nu2, X2, H01, H02, s, u, bwi, q) {
    .Call(`_JMH_GetCIF1CR`, gamma1, gamma2, alpha1, alpha2, nu1, nu2, X2, H01, H02, s, u, bwi, q)
}

GetCIF2CR <- function(gamma1, gamma2, alpha1, alpha2, nu1, nu2, X2, H01, H02, s, u, bwi, q) {
    .Call(`_JMH_GetCIF2CR`, gamma1, gamma2, alpha1, alpha2, nu1, nu2, X2, H01, H02, s, u, bwi, q)
}

getCov <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW) {
    .Call(`_JMH_getCov`, beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW)
}

getCovSF <- function(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW) {
    .Call(`_JMH_getCovSF`, beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW)
}

getEC <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, CUH02, HAZ01, HAZ02) {
    .Call(`_JMH_getEC`, beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, CUH02, HAZ01, HAZ02)
}

getECIF <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, nu1, nu2, Sig, Z, X1, W, Y, X2, H01, H02, xsmatrix, wsmatrix, CH01, CH02, s, u) {
    .Call(`_JMH_getECIF`, beta, tau, gamma1, gamma2, alpha1, alpha2, nu1, nu2, Sig, Z, X1, W, Y, X2, H01, H02, xsmatrix, wsmatrix, CH01, CH02, s, u)
}

getECIFad <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, nu1, nu2, Sig, Z, X1, W, Y, X2, H01, H02, xsmatrix, wsmatrix, CH01, CH02, s, u, Posmean, Poscov) {
    .Call(`_JMH_getECIFad`, beta, tau, gamma1, gamma2, alpha1, alpha2, nu1, nu2, Sig, Z, X1, W, Y, X2, H01, H02, xsmatrix, wsmatrix, CH01, CH02, s, u, Posmean, Poscov)
}

getECSF <- function(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01) {
    .Call(`_JMH_getECSF`, beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01)
}

getECSFad <- function(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01, Posbwi, Poscov) {
    .Call(`_JMH_getECSFad`, beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01, Posbwi, Poscov)
}

getECad <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, CUH02, HAZ01, HAZ02, Posbwi, Poscov) {
    .Call(`_JMH_getECad`, beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, CUH02, HAZ01, HAZ02, Posbwi, Poscov)
}

getES <- function(beta, tau, gamma1, alpha1, nu1, Sig, Z, X1, W, Y, X2, xsmatrix, wsmatrix, CH0s, CH0u) {
    .Call(`_JMH_getES`, beta, tau, gamma1, alpha1, nu1, Sig, Z, X1, W, Y, X2, xsmatrix, wsmatrix, CH0s, CH0u)
}

getESad <- function(beta, tau, gamma1, alpha1, nu1, Sig, Z, X1, W, Y, X2, xsmatrix, wsmatrix, CH0s, CH0u, Posmean, Poscov) {
    .Call(`_JMH_getESad`, beta, tau, gamma1, alpha1, nu1, Sig, Z, X1, W, Y, X2, xsmatrix, wsmatrix, CH0s, CH0u, Posmean, Poscov)
}

getHazard <- function(CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02) {
    .Call(`_JMH_getHazard`, CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02)
}

getHazardSF <- function(CumuH01, survtime, cmprsk, H01, CUH01, HAZ01) {
    .Call(`_JMH_getHazardSF`, CumuH01, survtime, cmprsk, H01, CUH01, HAZ01)
}

getMC <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW) {
    .Call(`_JMH_getMC`, beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW)
}

getMCSF <- function(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW) {
    .Call(`_JMH_getMCSF`, beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW)
}

getloglikeC <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, CUH02, HAZ01, HAZ02) {
    .Call(`_JMH_getloglikeC`, beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, CUH02, HAZ01, HAZ02)
}

getloglikeCSF <- function(beta, tau, gamma1, alpha1, vee1, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01) {
    .Call(`_JMH_getloglikeCSF`, beta, tau, gamma1, alpha1, vee1, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01)
}

getloglikeCSFad <- function(beta, tau, gamma1, alpha1, vee1, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01, Posbwi, Poscov) {
    .Call(`_JMH_getloglikeCSFad`, beta, tau, gamma1, alpha1, vee1, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01, Posbwi, Poscov)
}

getloglikeCad <- function(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, CUH02, HAZ01, HAZ02, Posbwi, Poscov) {
    .Call(`_JMH_getloglikeCad`, beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, CUH02, HAZ01, HAZ02, Posbwi, Poscov)
}

