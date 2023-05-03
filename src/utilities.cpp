#include "utilities.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace std;

double Delta(double s_0, int m, int d) {
  return (1/s_0*log(exp(1)*m)+log(exp(1)*d/s_0));
}

Eigen::VectorXd tau(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXi &gindex, Eigen::VectorXi &gsize, double lambda, double s_0, int m, int p) {
    Eigen::VectorXd temp = beta+X.transpose()*(y-X*beta)/X.rows();
     for (int i = 0; i < p; i++) {
        if (abs(temp(i)) < lambda) temp(i) = 0.0;
     }
     for (int i = 0; i < m; i++) {
        if (temp.segment(gindex(i), gsize(i)).squaredNorm() < s_0*pow(lambda, 2)) {
          temp.segment(gindex(i), gsize(i)) = Eigen::VectorXd::Zero(gsize(i));
        }
     }
     return temp;
}


int group_support_size(Eigen::VectorXd &beta, Eigen::VectorXi &gindex, Eigen::VectorXi &gsize, int m) {
  int size = 0;
  Eigen::VectorXd temp = Eigen::VectorXd::Zero(m);
  for (int i = 0; i < m; i++) {
    temp(i) = beta.segment(gindex(i), gsize(i)).squaredNorm();
    if (temp(i) != 0 ) {
        size = size + 1;
    }
  }
  return size;
}

Eigen::VectorXi support_set(Eigen::VectorXd &beta, int p, int size) {
  Eigen::VectorXi temp = Eigen::VectorXi::Zero(size);
  int flag = 0;
  for (int i = 0; i < p; i++) {
    if (beta(i) != 0) {
      temp(flag) = i;
      flag = flag+1;
    }
  }
  return temp;
}

Eigen::VectorXd least_square(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, int p) {
  int size = (beta.array() != 0).count();
  if (size >= X.rows() || size == 0) {
    return beta;
  }
  else{
  Eigen::VectorXi set = support_set(beta, p, size);
  Eigen::MatrixXd X_temp = Eigen::MatrixXd::Zero(X.rows(), size);
  for (int i = 0; i < size; i++) {
    X_temp.col(i) = X.col(set(i));
  }
  Eigen::VectorXd temp = X_temp.colPivHouseholderQr().solve(y);
  Eigen::VectorXd beta_hat = Eigen::VectorXd::Zero(p);
  for (int i = 0; i < size; i++) {
    beta_hat(set(i)) = temp(i);
  }
  return(beta_hat);
  }
}


double IC(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXi &gindex, Eigen::VectorXi &gsize, double s_0, int n, int m, int p, int d, double delta_t, double ic_coef) {
  int size1 = group_support_size(beta, gindex, gsize, m);
  int size2 = (beta.array() != 0).count();
  double size3;
  if (size1 > size2/s_0) {
    size3 = size1;
  }
  else {
    size3 = size2/s_0;
  }
  double omega = size3*log(exp(1)*m/size3)+s_0*size3*log(exp(1)*d/s_0);
  double ic = (y-X*beta).squaredNorm()/n+ic_coef*0.5*omega*pow(delta_t, 2)/n;
  if (isnan(ic)) {
    return 1e10;
  }
  else {
    return ic;
  }
}

