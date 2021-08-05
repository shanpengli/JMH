// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// GetrisksetC
Rcpp::List GetrisksetC(const Eigen::MatrixXd& cdata);
RcppExport SEXP _JMH_GetrisksetC(SEXP cdataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type cdata(cdataSEXP);
    rcpp_result_gen = Rcpp::wrap(GetrisksetC(cdata));
    return rcpp_result_gen;
END_RCPP
}
// GetrisksetCSF
Rcpp::List GetrisksetCSF(const Eigen::MatrixXd& cdata);
RcppExport SEXP _JMH_GetrisksetCSF(SEXP cdataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type cdata(cdataSEXP);
    rcpp_result_gen = Rcpp::wrap(GetrisksetCSF(cdata));
    return rcpp_result_gen;
END_RCPP
}
// OLS
Rcpp::List OLS(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y);
RcppExport SEXP _JMH_OLS(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(OLS(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// MultVV
double MultVV(const Eigen::VectorXd& x, const Eigen::VectorXd& y);
RcppExport SEXP _JMH_MultVV(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(MultVV(x, y));
    return rcpp_result_gen;
END_RCPP
}
// MultVVoutprod
Eigen::MatrixXd MultVVoutprod(const Eigen::VectorXd& x);
RcppExport SEXP _JMH_MultVVoutprod(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(MultVVoutprod(x));
    return rcpp_result_gen;
END_RCPP
}
// MultVVinprod
double MultVVinprod(const Eigen::VectorXd& x);
RcppExport SEXP _JMH_MultVVinprod(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(MultVVinprod(x));
    return rcpp_result_gen;
END_RCPP
}
// CH
double CH(const Eigen::MatrixXd& H, double t);
RcppExport SEXP _JMH_CH(SEXP HSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type H(HSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(CH(H, t));
    return rcpp_result_gen;
END_RCPP
}
// HAZ
double HAZ(const Eigen::MatrixXd& H, double t);
RcppExport SEXP _JMH_HAZ(SEXP HSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type H(HSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(HAZ(H, t));
    return rcpp_result_gen;
END_RCPP
}
// MultMM
Eigen::MatrixXd MultMM(const Eigen::MatrixXd& x, const Eigen::MatrixXd& y);
RcppExport SEXP _JMH_MultMM(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(MultMM(x, y));
    return rcpp_result_gen;
END_RCPP
}
// haha
Rcpp::List haha(const Eigen::VectorXd& x, const Eigen::MatrixXd& y);
RcppExport SEXP _JMH_haha(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(haha(x, y));
    return rcpp_result_gen;
END_RCPP
}
// getCov
Rcpp::List getCov(const Eigen::VectorXd& beta, const Eigen::VectorXd& tau, const Eigen::VectorXd& gamma1, const Eigen::VectorXd& gamma2, const Eigen::VectorXd& alpha1, const Eigen::VectorXd& alpha2, const double vee1, const double vee2, const Eigen::MatrixXd& H01, const Eigen::MatrixXd& H02, const Eigen::MatrixXd& Sig, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& X1, const Eigen::MatrixXd& W, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X2, const Eigen::VectorXd& survtime, const Eigen::VectorXd& cmprsk, const Eigen::VectorXd& mdata, const Eigen::VectorXd& mdataS, const Eigen::VectorXd& FUNENW, const Eigen::MatrixXd& FUNBENW, const Eigen::MatrixXd& FUNBS, const Eigen::MatrixXd& FUNBW, const Eigen::VectorXd& FUNWS, const Eigen::MatrixXd& FUNBSENW, const Eigen::MatrixXd& FUNEC, const Eigen::MatrixXd& FUNBEC, const Eigen::MatrixXd& FUNBSEC, const Eigen::MatrixXd& FUNWEC, const Eigen::MatrixXd& FUNWSEC, const Eigen::MatrixXd& FUNB, const Eigen::VectorXd& FUNW);
RcppExport SEXP _JMH_getCov(SEXP betaSEXP, SEXP tauSEXP, SEXP gamma1SEXP, SEXP gamma2SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP vee1SEXP, SEXP vee2SEXP, SEXP H01SEXP, SEXP H02SEXP, SEXP SigSEXP, SEXP ZSEXP, SEXP X1SEXP, SEXP WSEXP, SEXP YSEXP, SEXP X2SEXP, SEXP survtimeSEXP, SEXP cmprskSEXP, SEXP mdataSEXP, SEXP mdataSSEXP, SEXP FUNENWSEXP, SEXP FUNBENWSEXP, SEXP FUNBSSEXP, SEXP FUNBWSEXP, SEXP FUNWSSEXP, SEXP FUNBSENWSEXP, SEXP FUNECSEXP, SEXP FUNBECSEXP, SEXP FUNBSECSEXP, SEXP FUNWECSEXP, SEXP FUNWSECSEXP, SEXP FUNBSEXP, SEXP FUNWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type vee1(vee1SEXP);
    Rcpp::traits::input_parameter< const double >::type vee2(vee2SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type H01(H01SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type H02(H02SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Sig(SigSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type survtime(survtimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type cmprsk(cmprskSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdata(mdataSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdataS(mdataSSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNENW(FUNENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBENW(FUNBENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBS(FUNBSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBW(FUNBWSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNWS(FUNWSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBSENW(FUNBSENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNEC(FUNECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBEC(FUNBECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBSEC(FUNBSECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNWEC(FUNWECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNWSEC(FUNWSECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNB(FUNBSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNW(FUNWSEXP);
    rcpp_result_gen = Rcpp::wrap(getCov(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW));
    return rcpp_result_gen;
END_RCPP
}
// getCovSF
Rcpp::List getCovSF(const Eigen::VectorXd& beta, const Eigen::VectorXd& tau, const Eigen::VectorXd& gamma1, const Eigen::VectorXd& alpha1, const double vee1, const Eigen::MatrixXd& H01, const Eigen::MatrixXd& Sig, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& X1, const Eigen::MatrixXd& W, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X2, const Eigen::VectorXd& survtime, const Eigen::VectorXd& cmprsk, const Eigen::VectorXd& mdata, const Eigen::VectorXd& mdataS, const Eigen::VectorXd& FUNENW, const Eigen::MatrixXd& FUNBENW, const Eigen::MatrixXd& FUNBS, const Eigen::MatrixXd& FUNBW, const Eigen::VectorXd& FUNWS, const Eigen::MatrixXd& FUNBSENW, const Eigen::MatrixXd& FUNEC, const Eigen::MatrixXd& FUNBEC, const Eigen::MatrixXd& FUNBSEC, const Eigen::MatrixXd& FUNWEC, const Eigen::MatrixXd& FUNWSEC, const Eigen::MatrixXd& FUNB, const Eigen::VectorXd& FUNW);
RcppExport SEXP _JMH_getCovSF(SEXP betaSEXP, SEXP tauSEXP, SEXP gamma1SEXP, SEXP alpha1SEXP, SEXP vee1SEXP, SEXP H01SEXP, SEXP SigSEXP, SEXP ZSEXP, SEXP X1SEXP, SEXP WSEXP, SEXP YSEXP, SEXP X2SEXP, SEXP survtimeSEXP, SEXP cmprskSEXP, SEXP mdataSEXP, SEXP mdataSSEXP, SEXP FUNENWSEXP, SEXP FUNBENWSEXP, SEXP FUNBSSEXP, SEXP FUNBWSEXP, SEXP FUNWSSEXP, SEXP FUNBSENWSEXP, SEXP FUNECSEXP, SEXP FUNBECSEXP, SEXP FUNBSECSEXP, SEXP FUNWECSEXP, SEXP FUNWSECSEXP, SEXP FUNBSEXP, SEXP FUNWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type vee1(vee1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type H01(H01SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Sig(SigSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type survtime(survtimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type cmprsk(cmprskSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdata(mdataSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdataS(mdataSSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNENW(FUNENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBENW(FUNBENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBS(FUNBSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBW(FUNBWSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNWS(FUNWSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBSENW(FUNBSENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNEC(FUNECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBEC(FUNBECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBSEC(FUNBSECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNWEC(FUNWECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNWSEC(FUNWSECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNB(FUNBSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNW(FUNWSEXP);
    rcpp_result_gen = Rcpp::wrap(getCovSF(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW));
    return rcpp_result_gen;
END_RCPP
}
// getEC
Rcpp::List getEC(const Eigen::VectorXd& beta, const Eigen::VectorXd& tau, const Eigen::VectorXd& gamma1, const Eigen::VectorXd& gamma2, const Eigen::VectorXd& alpha1, const Eigen::VectorXd& alpha2, const double vee1, const double vee2, const Eigen::MatrixXd& H01, const Eigen::MatrixXd& H02, const Eigen::MatrixXd& Sig, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& X1, const Eigen::MatrixXd& W, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X2, const Eigen::VectorXd& survtime, const Eigen::VectorXd& cmprsk, const Eigen::VectorXd& mdata, const Eigen::VectorXd& mdataS, const Eigen::MatrixXd& xsmatrix, const Eigen::MatrixXd& wsmatrix, const Eigen::VectorXd& CUH01, const Eigen::VectorXd& CUH02, const Eigen::VectorXd& HAZ01, const Eigen::VectorXd& HAZ02);
RcppExport SEXP _JMH_getEC(SEXP betaSEXP, SEXP tauSEXP, SEXP gamma1SEXP, SEXP gamma2SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP vee1SEXP, SEXP vee2SEXP, SEXP H01SEXP, SEXP H02SEXP, SEXP SigSEXP, SEXP ZSEXP, SEXP X1SEXP, SEXP WSEXP, SEXP YSEXP, SEXP X2SEXP, SEXP survtimeSEXP, SEXP cmprskSEXP, SEXP mdataSEXP, SEXP mdataSSEXP, SEXP xsmatrixSEXP, SEXP wsmatrixSEXP, SEXP CUH01SEXP, SEXP CUH02SEXP, SEXP HAZ01SEXP, SEXP HAZ02SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type vee1(vee1SEXP);
    Rcpp::traits::input_parameter< const double >::type vee2(vee2SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type H01(H01SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type H02(H02SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Sig(SigSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type survtime(survtimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type cmprsk(cmprskSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdata(mdataSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdataS(mdataSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type xsmatrix(xsmatrixSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type wsmatrix(wsmatrixSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type CUH01(CUH01SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type CUH02(CUH02SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type HAZ01(HAZ01SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type HAZ02(HAZ02SEXP);
    rcpp_result_gen = Rcpp::wrap(getEC(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, CUH02, HAZ01, HAZ02));
    return rcpp_result_gen;
END_RCPP
}
// getECSF
Rcpp::List getECSF(const Eigen::VectorXd& beta, const Eigen::VectorXd& tau, const Eigen::VectorXd& gamma1, const Eigen::VectorXd& alpha1, const double vee1, const Eigen::MatrixXd& H01, const Eigen::MatrixXd& Sig, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& X1, const Eigen::MatrixXd& W, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X2, const Eigen::VectorXd& survtime, const Eigen::VectorXd& cmprsk, const Eigen::VectorXd& mdata, const Eigen::VectorXd& mdataS, const Eigen::MatrixXd& xsmatrix, const Eigen::MatrixXd& wsmatrix, const Eigen::VectorXd& CUH01, const Eigen::VectorXd& HAZ01);
RcppExport SEXP _JMH_getECSF(SEXP betaSEXP, SEXP tauSEXP, SEXP gamma1SEXP, SEXP alpha1SEXP, SEXP vee1SEXP, SEXP H01SEXP, SEXP SigSEXP, SEXP ZSEXP, SEXP X1SEXP, SEXP WSEXP, SEXP YSEXP, SEXP X2SEXP, SEXP survtimeSEXP, SEXP cmprskSEXP, SEXP mdataSEXP, SEXP mdataSSEXP, SEXP xsmatrixSEXP, SEXP wsmatrixSEXP, SEXP CUH01SEXP, SEXP HAZ01SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< const double >::type vee1(vee1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type H01(H01SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Sig(SigSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type survtime(survtimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type cmprsk(cmprskSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdata(mdataSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdataS(mdataSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type xsmatrix(xsmatrixSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type wsmatrix(wsmatrixSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type CUH01(CUH01SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type HAZ01(HAZ01SEXP);
    rcpp_result_gen = Rcpp::wrap(getECSF(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, xsmatrix, wsmatrix, CUH01, HAZ01));
    return rcpp_result_gen;
END_RCPP
}
// getHazard
int getHazard(const Eigen::Map<Eigen::VectorXd>& CumuH01, const Eigen::Map<Eigen::VectorXd>& CumuH02, const Eigen::Map<Eigen::VectorXd>& survtime, const Eigen::Map<Eigen::VectorXd>& cmprsk, const Eigen::Map<Eigen::MatrixXd>& H01, const Eigen::Map<Eigen::MatrixXd>& H02, Eigen::Map<Eigen::VectorXd>& CUH01, Eigen::Map<Eigen::VectorXd>& CUH02, Eigen::Map<Eigen::VectorXd>& HAZ01, Eigen::Map<Eigen::VectorXd>& HAZ02);
RcppExport SEXP _JMH_getHazard(SEXP CumuH01SEXP, SEXP CumuH02SEXP, SEXP survtimeSEXP, SEXP cmprskSEXP, SEXP H01SEXP, SEXP H02SEXP, SEXP CUH01SEXP, SEXP CUH02SEXP, SEXP HAZ01SEXP, SEXP HAZ02SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type CumuH01(CumuH01SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type CumuH02(CumuH02SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type survtime(survtimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type cmprsk(cmprskSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type H01(H01SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type H02(H02SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd>& >::type CUH01(CUH01SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd>& >::type CUH02(CUH02SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd>& >::type HAZ01(HAZ01SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd>& >::type HAZ02(HAZ02SEXP);
    rcpp_result_gen = Rcpp::wrap(getHazard(CumuH01, CumuH02, survtime, cmprsk, H01, H02, CUH01, CUH02, HAZ01, HAZ02));
    return rcpp_result_gen;
END_RCPP
}
// getHazardSF
int getHazardSF(const Eigen::Map<Eigen::VectorXd>& CumuH01, const Eigen::Map<Eigen::VectorXd>& survtime, const Eigen::Map<Eigen::VectorXd>& cmprsk, const Eigen::Map<Eigen::MatrixXd>& H01, Eigen::Map<Eigen::VectorXd>& CUH01, Eigen::Map<Eigen::VectorXd>& HAZ01);
RcppExport SEXP _JMH_getHazardSF(SEXP CumuH01SEXP, SEXP survtimeSEXP, SEXP cmprskSEXP, SEXP H01SEXP, SEXP CUH01SEXP, SEXP HAZ01SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type CumuH01(CumuH01SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type survtime(survtimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type cmprsk(cmprskSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type H01(H01SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd>& >::type CUH01(CUH01SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd>& >::type HAZ01(HAZ01SEXP);
    rcpp_result_gen = Rcpp::wrap(getHazardSF(CumuH01, survtime, cmprsk, H01, CUH01, HAZ01));
    return rcpp_result_gen;
END_RCPP
}
// getMC
Rcpp::List getMC(Eigen::VectorXd& beta, Eigen::VectorXd& tau, Eigen::VectorXd& gamma1, Eigen::VectorXd& gamma2, Eigen::VectorXd& alpha1, Eigen::VectorXd& alpha2, double vee1, double vee2, Eigen::MatrixXd& H01, Eigen::MatrixXd& H02, Eigen::MatrixXd& Sig, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& X1, const Eigen::MatrixXd& W, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X2, const Eigen::VectorXd& survtime, const Eigen::VectorXd& cmprsk, const Eigen::VectorXd& mdata, const Eigen::VectorXd& mdataS, const Eigen::VectorXd& FUNENW, const Eigen::MatrixXd& FUNBENW, const Eigen::MatrixXd& FUNBS, const Eigen::MatrixXd& FUNBW, const Eigen::VectorXd& FUNWS, const Eigen::MatrixXd& FUNBSENW, const Eigen::MatrixXd& FUNEC, const Eigen::MatrixXd& FUNBEC, const Eigen::MatrixXd& FUNBSEC, const Eigen::MatrixXd& FUNWEC, const Eigen::MatrixXd& FUNWSEC, const Eigen::MatrixXd& FUNB, const Eigen::VectorXd& FUNW);
RcppExport SEXP _JMH_getMC(SEXP betaSEXP, SEXP tauSEXP, SEXP gamma1SEXP, SEXP gamma2SEXP, SEXP alpha1SEXP, SEXP alpha2SEXP, SEXP vee1SEXP, SEXP vee2SEXP, SEXP H01SEXP, SEXP H02SEXP, SEXP SigSEXP, SEXP ZSEXP, SEXP X1SEXP, SEXP WSEXP, SEXP YSEXP, SEXP X2SEXP, SEXP survtimeSEXP, SEXP cmprskSEXP, SEXP mdataSEXP, SEXP mdataSSEXP, SEXP FUNENWSEXP, SEXP FUNBENWSEXP, SEXP FUNBSSEXP, SEXP FUNBWSEXP, SEXP FUNWSSEXP, SEXP FUNBSENWSEXP, SEXP FUNECSEXP, SEXP FUNBECSEXP, SEXP FUNBSECSEXP, SEXP FUNWECSEXP, SEXP FUNWSECSEXP, SEXP FUNBSEXP, SEXP FUNWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< double >::type vee1(vee1SEXP);
    Rcpp::traits::input_parameter< double >::type vee2(vee2SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type H01(H01SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type H02(H02SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Sig(SigSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type survtime(survtimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type cmprsk(cmprskSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdata(mdataSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdataS(mdataSSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNENW(FUNENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBENW(FUNBENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBS(FUNBSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBW(FUNBWSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNWS(FUNWSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBSENW(FUNBSENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNEC(FUNECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBEC(FUNBECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBSEC(FUNBSECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNWEC(FUNWECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNWSEC(FUNWSECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNB(FUNBSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNW(FUNWSEXP);
    rcpp_result_gen = Rcpp::wrap(getMC(beta, tau, gamma1, gamma2, alpha1, alpha2, vee1, vee2, H01, H02, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW));
    return rcpp_result_gen;
END_RCPP
}
// getMCSF
Rcpp::List getMCSF(Eigen::VectorXd& beta, Eigen::VectorXd& tau, Eigen::VectorXd& gamma1, Eigen::VectorXd& alpha1, double vee1, Eigen::MatrixXd& H01, Eigen::MatrixXd& Sig, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& X1, const Eigen::MatrixXd& W, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X2, const Eigen::VectorXd& survtime, const Eigen::VectorXd& cmprsk, const Eigen::VectorXd& mdata, const Eigen::VectorXd& mdataS, const Eigen::VectorXd& FUNENW, const Eigen::MatrixXd& FUNBENW, const Eigen::MatrixXd& FUNBS, const Eigen::MatrixXd& FUNBW, const Eigen::VectorXd& FUNWS, const Eigen::MatrixXd& FUNBSENW, const Eigen::MatrixXd& FUNEC, const Eigen::MatrixXd& FUNBEC, const Eigen::MatrixXd& FUNBSEC, const Eigen::MatrixXd& FUNWEC, const Eigen::MatrixXd& FUNWSEC, const Eigen::MatrixXd& FUNB, const Eigen::VectorXd& FUNW);
RcppExport SEXP _JMH_getMCSF(SEXP betaSEXP, SEXP tauSEXP, SEXP gamma1SEXP, SEXP alpha1SEXP, SEXP vee1SEXP, SEXP H01SEXP, SEXP SigSEXP, SEXP ZSEXP, SEXP X1SEXP, SEXP WSEXP, SEXP YSEXP, SEXP X2SEXP, SEXP survtimeSEXP, SEXP cmprskSEXP, SEXP mdataSEXP, SEXP mdataSSEXP, SEXP FUNENWSEXP, SEXP FUNBENWSEXP, SEXP FUNBSSEXP, SEXP FUNBWSEXP, SEXP FUNWSSEXP, SEXP FUNBSENWSEXP, SEXP FUNECSEXP, SEXP FUNBECSEXP, SEXP FUNBSECSEXP, SEXP FUNWECSEXP, SEXP FUNWSECSEXP, SEXP FUNBSEXP, SEXP FUNWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< double >::type vee1(vee1SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type H01(H01SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Sig(SigSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type survtime(survtimeSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type cmprsk(cmprskSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdata(mdataSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mdataS(mdataSSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNENW(FUNENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBENW(FUNBENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBS(FUNBSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBW(FUNBWSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNWS(FUNWSSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBSENW(FUNBSENWSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNEC(FUNECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBEC(FUNBECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNBSEC(FUNBSECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNWEC(FUNWECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNWSEC(FUNWSECSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type FUNB(FUNBSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type FUNW(FUNWSEXP);
    rcpp_result_gen = Rcpp::wrap(getMCSF(beta, tau, gamma1, alpha1, vee1, H01, Sig, Z, X1, W, Y, X2, survtime, cmprsk, mdata, mdataS, FUNENW, FUNBENW, FUNBS, FUNBW, FUNWS, FUNBSENW, FUNEC, FUNBEC, FUNBSEC, FUNWEC, FUNWSEC, FUNB, FUNW));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_JMH_GetrisksetC", (DL_FUNC) &_JMH_GetrisksetC, 1},
    {"_JMH_GetrisksetCSF", (DL_FUNC) &_JMH_GetrisksetCSF, 1},
    {"_JMH_OLS", (DL_FUNC) &_JMH_OLS, 2},
    {"_JMH_MultVV", (DL_FUNC) &_JMH_MultVV, 2},
    {"_JMH_MultVVoutprod", (DL_FUNC) &_JMH_MultVVoutprod, 1},
    {"_JMH_MultVVinprod", (DL_FUNC) &_JMH_MultVVinprod, 1},
    {"_JMH_CH", (DL_FUNC) &_JMH_CH, 2},
    {"_JMH_HAZ", (DL_FUNC) &_JMH_HAZ, 2},
    {"_JMH_MultMM", (DL_FUNC) &_JMH_MultMM, 2},
    {"_JMH_haha", (DL_FUNC) &_JMH_haha, 2},
    {"_JMH_getCov", (DL_FUNC) &_JMH_getCov, 33},
    {"_JMH_getCovSF", (DL_FUNC) &_JMH_getCovSF, 29},
    {"_JMH_getEC", (DL_FUNC) &_JMH_getEC, 26},
    {"_JMH_getECSF", (DL_FUNC) &_JMH_getECSF, 20},
    {"_JMH_getHazard", (DL_FUNC) &_JMH_getHazard, 10},
    {"_JMH_getHazardSF", (DL_FUNC) &_JMH_getHazardSF, 6},
    {"_JMH_getMC", (DL_FUNC) &_JMH_getMC, 33},
    {"_JMH_getMCSF", (DL_FUNC) &_JMH_getMCSF, 29},
    {NULL, NULL, 0}
};

RcppExport void R_init_JMH(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
