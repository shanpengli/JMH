#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

double  MultVV(const Eigen::VectorXd & x, const Eigen::VectorXd & y);

Eigen::MatrixXd MultVVoutprod(const Eigen::VectorXd & x);

double MultVVinprod(const Eigen::VectorXd & x);

double CH(const Eigen::MatrixXd & H, double t);

double HAZ(const Eigen::MatrixXd & H, double t);

double getdeterminant(const Eigen::MatrixXd & H);

Eigen::MatrixXd MultVV2outprod(const Eigen::VectorXd & x, const Eigen::VectorXd & y);

double GetCIF1CR(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                 const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                 const double nu1, const double nu2,
                 const Eigen::VectorXd & X2,
                 const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                 const double s, const double u, const Eigen::VectorXd & bwi,
                 const int q);

double GetCIF2CR(const Eigen::VectorXd & gamma1, const Eigen::VectorXd & gamma2, 
                 const Eigen::VectorXd & alpha1, const Eigen::VectorXd & alpha2, 
                 const double nu1, const double nu2,
                 const Eigen::VectorXd & X2,
                 const Eigen::MatrixXd & H01, const Eigen::MatrixXd & H02,
                 const double s, const double u, const Eigen::VectorXd & bwi,
                 const int q);
