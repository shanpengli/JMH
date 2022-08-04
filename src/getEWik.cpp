#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getEWik(const Eigen::VectorXd & beta, const Eigen::VectorXd & tau, 
                   const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                   const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                   const double vee1, const double vee2, const Eigen::MatrixXd & Sig, 
                   const Eigen::MatrixXd & Z, const Eigen::MatrixXd & X1, 
                   const Eigen::MatrixXd & W, const Eigen::VectorXd & Y, 
                   const Eigen::MatrixXd & X2, const Eigen::VectorXd & survtime, 
                   const Eigen::VectorXd & cmprsk, const Eigen::VectorXd & mdata, 
                   const Eigen::VectorXd & mdataS, const Eigen::MatrixXd & xsmatrix, 
                   const Eigen::MatrixXd & wsmatrix, 
                   const Eigen::VectorXd & CUH01, 
                   const Eigen::VectorXd & CUH02, 
                   const Eigen::VectorXd & HAZ01, 
                   const Eigen::VectorXd & HAZ02){ 
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,t,db;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int k=mdata.size();
  int p1a=Z.cols();
  
  double dem,cuh01,cuh02,haz01,haz02,xgamma1,xgamma2,temp,mu,sigma,zb,wi;
  Eigen::VectorXd bwi(p1a+1);
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd weightbwi(p1a+1);
  Eigen::VectorXd ri(p1a+1);
  
  /* Define functions Wik*/
  //exp(alpha*b + vee*w)
  Eigen::MatrixXd FUNEC = Eigen::MatrixXd::Zero(2,k);
  
  int point=wsmatrix.rows();
  
  for(j=0;j<k;j++)
  {
    dem=0;
    q=mdata(j);
    cuh01=CUH01(j);
    cuh02=CUH02(j);
    haz01=HAZ01(j);
    haz02=HAZ02(j);
    
    xgamma1=MultVV(X2.row(j),gamma1);
    xgamma2=MultVV(X2.row(j),gamma2);
    
    for (db=0;db<point;db++) {
      
      bwi = xsmatrix.row(db);
      weightbwi = wsmatrix.row(db);
      ri = sqrt(2)*SigSQRT*bwi;
      temp=exp(10);
      
      for (i=0;i<p1a;i++) bi(i)=ri(i);
      wi=ri(p1a);
      
      for (i=0;i<q;i++) {
        mu=MultVV(X1.row(mdataS(j)-1+i),beta);
        zb=MultVV(Z.row(mdataS(j)-1+i),bi);
        sigma=exp(MultVV(W.row(mdataS(j)-1+i),tau) + wi);
        temp*=1/sqrt(sigma)*exp(-1/(2*sigma)*pow((Y(mdataS(j)-1+i) - mu - zb), 2)); 
      }
      
      if(cmprsk(j)==1)  temp*=haz01*exp(xgamma1+MultVV(alpha1,bi)+vee1*wi);
      if(cmprsk(j)==2)  temp*=haz02*exp(xgamma2+MultVV(alpha2,bi)+vee2*wi);
      
      temp*=exp(0-cuh01*exp(xgamma1+MultVV(alpha1,bi)+vee1*wi)-cuh02*exp(xgamma2+MultVV(alpha2,bi)+vee2*wi));
      for (i=0;i<(p1a+1);i++) temp*=weightbwi(i);
      
      dem+=temp;
      
      //calculate FUNEC
      FUNEC(0,j)+=temp*exp(MultVV(alpha1,bi)+vee1*wi);
      FUNEC(1,j)+=temp*exp(MultVV(alpha2,bi)+vee2*wi);
    }
    
    if(dem==0) {
      Rprintf("E step ran into issue for the %dth subject. Program stops.\n", j);
      return ( 100.0 );
    } 
    
    FUNEC.col(j)/=dem;
    
  }
  
  return Rcpp::List::create(Rcpp::Named("FUNEC")=FUNEC);
}
