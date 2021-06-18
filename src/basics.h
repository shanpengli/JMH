#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

double  MultVV(const Eigen::VectorXd & x, const Eigen::VectorXd & y);

Eigen::MatrixXd MultVVoutprod(const Eigen::VectorXd & x);

double MultVVinprod(const Eigen::VectorXd & x);

double CH(const Eigen::MatrixXd & H, double t);

double HAZ(const Eigen::MatrixXd & H, double t);