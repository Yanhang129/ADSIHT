#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <algorithm>
#include <vector>
#include <iostream>
using namespace Rcpp;
using namespace std;

void constant_warning_ith_variable(int i)
{
  Rcout << "Warning: the variable " << i + 1 << " is constant. ";
  Rcout << "It may cause NAN in the result. Please drop this variable or disable the normalization.\n";
}

void Normalize(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXd &meanx, double &meany,
               Eigen::VectorXd &normx) {
  int n = X.rows();
  int p = X.cols();
  Eigen::VectorXd tmp(n);
  for (int i = 0; i < p; i++) {
    meanx(i) = weights.dot(X.col(i)) / double(n);
  }
  meany = (y.dot(weights)) / double(n);
  for (int i = 0; i < p; i++) {
    X.col(i) = X.col(i).array() - meanx(i);
  }
  y = y.array() - meany;

  for (int i = 0; i < p; i++) {
    tmp = X.col(i);
    tmp = tmp.array().square();
    normx(i) = sqrt(weights.dot(tmp));
    if (normx(i) == 0) {
      constant_warning_ith_variable(i);
    }
  }
  for (int i = 0; i < p; i++) {
    X.col(i) = sqrt(double(n)) * X.col(i) / normx(i);
  }
}
