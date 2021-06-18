// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include "basics.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

// another simple example: outer product of a vector,
// returning a matrix
//
// [[Rcpp::export]]
double  MultVV(const Eigen::VectorXd & x, const Eigen::VectorXd & y) {
    double v = x.transpose() * y;
    return v;
}

// [[Rcpp::export]]
Eigen::MatrixXd MultVVoutprod(const Eigen::VectorXd & x) {
    Eigen::MatrixXd m = x * x.transpose();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double  MultVVinprod(const Eigen::VectorXd & x) {
    double v = x.transpose() * x;
    return v;
}

// [[Rcpp::export]]
double  CH(const Eigen::MatrixXd & H, double t) {
    
    int a = H.rows();
    int i;
    double ch;
    if (t < H(0,0)) ch=0;
    else {
        ch=0;
        i=0;
        do {
            ch+=H(i, 2);
            i+=1;
            
        } while (i<=a-1 && t>= H(i,0));
        }
    
    return (ch);
    
}

// [[Rcpp::export]]
double HAZ(const Eigen::MatrixXd & H, double t) {
    
    int a = H.rows();
    int i;
    double temp=0;
    for (i=0;i<a;i++) {
        if (t == H(i, 0)) temp = H(i,2);
        }
    
    return(temp);
    
    }

// [[Rcpp::export]]
Eigen::MatrixXd MultMM(const Eigen::MatrixXd & x, const Eigen::MatrixXd & y) {
    Eigen::MatrixXd v = x * y;
    return v;
}

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List haha(const Eigen::VectorXd & x, const Eigen::MatrixXd & y) {
    double aa = MultVVinprod(x);
    Eigen::MatrixXd m = MultVVoutprod(x);
    Eigen::LLT<Eigen::MatrixXd> LLT_of_K(y);
    if ( !LLT_of_K.info() ) { 
        Eigen::MatrixXd Rooti = LLT_of_K.matrixL();
        Eigen::MatrixXd Rootiinv = Rooti.inverse();
        Eigen::VectorXd v = x + x;
        return Rcpp::List::create(Rcpp::Named("outer")=aa,
                                  Rcpp::Named("inner")=v,
                                  Rcpp::Named("chol")=Rooti,
                                  Rcpp::Named("cholinv")=Rootiinv);
    } else {
        return ( -1.0 );
    }
    
    
}
