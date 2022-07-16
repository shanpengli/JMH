#include <RcppEigen.h>
#include "basics.h"

// [[Rcpp::depends(RcppEigen)]]

//
// [[Rcpp::export]]
Rcpp::List getTruePik(const Eigen::VectorXd & Y, 
                      const Eigen::VectorXd & xbeta, const Eigen::VectorXd & wtau, 
                      const Eigen::VectorXd & lam1, const Eigen::VectorXd & lam2, 
                      const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                      const double vee1, const double vee2, const Eigen::MatrixXd & Sig, 
                      const int ni, const Eigen::VectorXd & mdataS, 
                      const Eigen::MatrixXd & xsmatrix, const Eigen::MatrixXd & wsmatrix, 
                      const Eigen::VectorXd & landmarktime){ 
  
  //calculate the square root of random effect covariance matrix 
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sig, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd eigenSQ = svd.singularValues();
  int i,j,q,t,db;
  for (i=0;i<eigenSQ.size();i++) {
    eigenSQ(i) = sqrt(eigenSQ(i));
  }
  Eigen::MatrixXd SigSQRT  = svd.matrixU() * eigenSQ.asDiagonal() * svd.matrixV().transpose();
  
  int k=lam1.size();
  int p1a=1;
  
  double dem,temp,mu,sigma,zb,wi,CIFtemp,Stemp,lam1temp,lam2temp,timediff,survdiff;
  Eigen::VectorXd bi(p1a);
  Eigen::VectorXd bwi(p1a+1);
  Eigen::VectorXd weightbwi(p1a+1);
  Eigen::VectorXd ri(p1a+1);
  
  int landmarklength = landmarktime.size();
  Eigen::MatrixXd P1us = Eigen::MatrixXd::Zero(k, landmarklength-1);
  Eigen::MatrixXd P2us = Eigen::MatrixXd::Zero(k, landmarklength-1);
  
  int point=wsmatrix.rows();
  
  for(j=0;j<k;j++)
  {
    dem=0;
    for (db=0;db<point;db++) {

      bwi = xsmatrix.row(db);
      weightbwi = wsmatrix.row(db);
      ri = sqrt(2)*SigSQRT*bwi;
      temp=1;

      for (i=0;i<p1a;i++) bi(i)=ri(i);
      wi=ri(p1a);

      for (i=0;i<ni;i++) {
        mu=xbeta(mdataS(j)-1+i);
        zb=bi(0);
        sigma=exp(wtau(mdataS(j)-1+i) + wi);
        temp*=1/sqrt(sigma)*exp(-1/(2*sigma)*pow((Y(mdataS(j)-1+i) - mu - zb), 2));
      }

      temp=exp(0-lam1(j)*exp(MultVV(alpha1,bi)+vee1*wi)*landmarktime(0)
                 -lam2(j)*exp(MultVV(alpha2,bi)+vee2*wi)*landmarktime(0));

      for (i=0;i<(p1a+1);i++) temp*=weightbwi(i);

      dem+=temp;

      // calculate CIF and S
      lam1temp=lam1(j)*exp(MultVV(alpha1,bi)+vee1*wi);
      lam2temp=lam2(j)*exp(MultVV(alpha2,bi)+vee2*wi);
      for (t=0;t<(landmarklength-1);t++)
      {
        timediff=landmarktime(t+1)-landmarktime(0);
        survdiff=1-exp(0-(lam1temp+lam2temp)*timediff);
        CIFtemp=lam1temp/(lam1temp+lam2temp)*survdiff;
        P1us(j,t)+=temp*CIFtemp;
        CIFtemp=lam2temp/(lam1temp+lam2temp)*survdiff;
        P2us(j,t)+=temp*CIFtemp;
      }

    }

    if(dem==0) {
      Rprintf("E step ran into issue for the %dth subject. Program stops.\n", j);
      return ( 100.0 );
    }

    P1us.row(j)/=dem;
    P2us.row(j)/=dem;
  }
  
  return Rcpp::List::create(Rcpp::Named("P1us")=P1us,
                            Rcpp::Named("P2us")=P2us);
}


