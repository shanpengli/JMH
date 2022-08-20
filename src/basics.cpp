// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include "basics.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operation on them
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

// [[Rcpp::export]]
Eigen::MatrixXd MultVV2outprod(const Eigen::VectorXd & x, const Eigen::VectorXd & y) {
    Eigen::MatrixXd m = x * y.transpose();
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

// [[Rcpp::export]]
double GetCIF1CR(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                 const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                 const double nu1, const double nu2,
                 const Eigen::VectorXd & X2,
                 const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                 const double s, const double u, const Eigen::VectorXd & bwi,
                 const int q){
  
  int a = H01.rows();
  int b = H02.rows();
  
  int i = 0;
  Eigen::VectorXd CH01 = Eigen::VectorXd::Zero(a);
  
  double hazard=0;
  for (i=0;i<a;i++) {
    hazard = hazard + H01(i, 2);
    CH01(i) = hazard;
  }
  
  Eigen::VectorXd FCH02 = Eigen::VectorXd::Zero(a);
  int count = 0;
  i = 0;
  while (count < b && i < a) {
    if (H02(count, 0) < H01(i, 0)) {
      FCH02(i) = FCH02(i) + H02(count, 2);
      count++;
    } else {
      i++;
    }
  }
  Eigen::VectorXd CH02 = Eigen::VectorXd::Zero(a);
  hazard=0;
  for (i=0;i<a;i++) {
    hazard = hazard + FCH02(i);
    CH02(i) = hazard;
  }
  
  double CIF1=0;
  int p1a=q-1;
  Eigen::VectorXd bi = Eigen::VectorXd::Zero(p1a);
  
  for (i=0;i<p1a;i++) bi(i)=bwi(i);
  double wi=bwi(p1a);
  
  for (i=0;i<a;i++) {
    if (s < H01(i, 0) && u >= H01(i, 0)) {
      if (i >= 1) {
        CIF1 = CIF1 + H01(i, 2)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi)+nu1*wi)*
          exp(-CH01(i-1)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi)+nu1*wi)-
          CH02(i-1)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi)+nu2*wi));
      } else {
        CIF1 = CIF1 + H01(i, 2)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi)+nu1*wi);
      }
    } else continue;
    
  }
  
  return CIF1;
}

// [[Rcpp::export]]
double GetCIF2CR(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                 const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                 const double nu1, const double nu2,
                 const Eigen::VectorXd & X2,
                 const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                 const double s, const double u, const Eigen::VectorXd & bwi,
                 const int q){
  
  int a = H01.rows();
  int b = H02.rows();
  
  int i = 0;
  Eigen::VectorXd CH02 = Eigen::VectorXd::Zero(b);
  
  double hazard=0;
  for (i=0;i<b;i++) {
    hazard = hazard + H02(i, 2);
    CH02(i) = hazard;
  }
  
  Eigen::VectorXd FCH01 = Eigen::VectorXd::Zero(b);
  int count = 0;
  i = 0;
  while (count < a && i < b) {
    if (H01(count, 0) < H02(i, 0)) {
      FCH01(i) = FCH01(i) + H01(count, 2);
      count++;
    } else {
      i++;
    }
  }
  Eigen::VectorXd CH01 = Eigen::VectorXd::Zero(b);
  hazard=0;
  for (i=0;i<b;i++) {
    hazard = hazard + FCH01(i);
    CH01(i) = hazard;
  }
  
  double CIF2=0;
  int p1a=q-1;
  Eigen::VectorXd bi = Eigen::VectorXd::Zero(p1a);
  
  for (i=0;i<p1a;i++) bi(i)=bwi(i);
  double wi=bwi(p1a);
  
  for (i=0;i<b;i++) {
    if (s < H02(i, 0) && u >= H02(i, 0)) {
      if (i >= 1) {
        CIF2 = CIF2 + H02(i, 2)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi)+nu2*wi)*
          exp(-CH01(i-1)*exp(MultVV(X2,gamma1)+MultVV(alpha1,bi)+nu1*wi)-
          CH02(i-1)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi)+nu2*wi));
      } else {
        CIF2 = CIF2 + H02(i, 2)*exp(MultVV(X2,gamma2)+MultVV(alpha2,bi)+nu2*wi);
      }
    } else continue;
    
  }
  
  return CIF2;
}



