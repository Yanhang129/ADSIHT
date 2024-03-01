#ifndef SRC_ALGORITHM_H
#define SRC_ALGORITHM_H
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>
#include "utilities.h"
#include "Data.h"
#include "math.h"
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]
using namespace std;

class Algorithm {
public:
  Data data;
  Eigen::VectorXd beta;
  Eigen::VectorXi A_out;
  int size;
  int group_size;
  double lam;
  double rho;
  double ic;
  double s0;

  Algorithm() = default;

  Algorithm(Data &data)
  {
    this->data = data;
    this->beta = Eigen::VectorXd::Zero(data.p);
    this->lam = 0.0;
    this->ic = 0.0;
    this->size = 0;
    this->s0 = 1;
  };

  void update_rho(double rho) {
    this->rho = rho;
  }

  double get_lambda() {
    return this->lam;
  }

  double get_s0() {
    return this->s0;
  }

  double get_ic() {
    return this->ic;
  }

  Eigen::VectorXd get_beta() {
    return this->beta;
  }

  double get_support_size() {
    return this->size;
  }

  double get_group_support_size() {
    return this->group_size;
  }

  Eigen::VectorXi get_A_out() {
    return this->A_out;
  }

  void fit_eta(double s_0, double ic_coef, double coef1, double coef2, double eta, int max_iter) {
    Eigen::MatrixXd x = data.x;
    Eigen::VectorXd y = data.y;
    int n = data.n;
    int p = data.p;
    int m = data.g_num;
    Eigen::VectorXi gindex = data.get_g_index();
    Eigen::VectorXi gsize = data.get_g_size();
    int d = gsize.maxCoeff();
    Eigen::VectorXd beta0 = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd beta1 = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd gradient = x.transpose()*(y-x*beta0)/n;

    double delta = Delta(s_0, m, d);
    double lam1, lam0 = pow(rho, ceil(s_0/3)-1)*max(std::sqrt(delta*y.squaredNorm()/n/n), ((x.transpose()*y/n).cwiseAbs()).maxCoeff());
    while(1) {
      beta1 = tau_eta(x, y, beta0, gradient, eta, gindex, gsize, lam0, s_0, m, p, max_iter);
      lam1 = rho*lam0;
      if (lam1 >= coef1*sqrt((y-x*beta1).squaredNorm()/n*delta/n)) {
        beta0 = beta1;
        gradient = x.transpose()*(y-x*beta0)/n;
        lam0 = lam1;
      } else {
        break;
      }
    }
    double delta_tbar = sqrt((y-x*beta1).squaredNorm()/n);
    double ic0, ic1 = IC(x, y, beta1, gindex, gsize, s_0, n, m, p, d, delta_tbar, ic_coef);
    this->ic = ic1;
    this->beta = beta1;
    this->lam = lam1;
    gradient = x.transpose()*(y-x*beta1)/n;
    while (1) {
      beta0 = tau_eta(x, y, beta1, gradient, eta, gindex, gsize, lam0, s_0, m, p, max_iter);
      if ((beta0.array() != 0).count() >= p/5 || (beta0.array() != 0).count() >= x.rows()) break;
      lam0 = rho*lam1;
      ic0 = IC(x, y, beta0, gindex, gsize, s_0, n, m, p, d, delta_tbar, ic_coef);
      if (ic0 < this->ic) {
        this->beta = beta0;
        this->lam = lam0;
        this->ic = ic0;
        beta1 = beta0;
        gradient = x.transpose()*(y-x*beta1)/n;
      }
      lam1 = lam0;
      if (lam1 < coef2*log(exp(1)*d/s_0)*delta_tbar/n) break;
    }
    this->s0 = s_0;
    this->group_size = group_support_size(this->beta, gindex, gsize, m);
    this->size = (this->beta.array() != 0).count();
    this->A_out = support_set(this->beta, p, this->size)+Eigen::VectorXi::Ones(this->size);
  }

  void fit1(double s_0, double ic_coef, double coef1, double coef2) {
    Eigen::MatrixXd x = data.x;
    Eigen::VectorXd y = data.y;
    int n = data.n;
    int p = data.p;
    int m = data.g_num;
    Eigen::VectorXi gindex = data.get_g_index();
    Eigen::VectorXi gsize = data.get_g_size();
    int d = gsize.maxCoeff();
    Eigen::VectorXd beta0 = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd beta1 = Eigen::VectorXd::Zero(p);
    double delta = Delta(s_0, m, d);
    double lam1, lam0 = pow(rho, ceil(s_0/3)-1)*max(std::sqrt(delta*y.squaredNorm()/n/n), ((x.transpose()*y/n).cwiseAbs()).maxCoeff());
    while(1) {
      beta1 = tau(x, y, beta0, gindex, gsize, lam0, s_0, m, p);
      lam1 = rho*rho*lam0;
//    coef1 = 2
    if (lam1 >= coef1*sqrt((y-x*beta1).squaredNorm()/n*delta/n)) {
        beta0 = beta1;
        lam0 = lam1;
      } else {
        break;
      }
    }
    beta1 = least_square(x, y, beta1, p);
    double delta_tbar = sqrt((y-x*beta1).squaredNorm()/n);
    double ic0, ic1 = IC(x, y, beta1, gindex, gsize, s_0, n, m, p, d, delta_tbar, ic_coef);
    this->ic = ic1;
    this->beta = beta1;
    this->lam = lam1;
    while (1) {
      beta0 = tau(x, y, beta1, gindex, gsize, lam1, s_0, m, p);
      beta0 = least_square(x, y, beta0, p);
      if ((beta0.array() != 0).count() >= p/5 || (beta0.array() != 0).count() >= x.rows()) break;
      lam0 = rho*lam1;
      ic0 = IC(x, y, beta0, gindex, gsize, s_0, n, m, p, d, delta_tbar, ic_coef);
      // Rcout<<"-IC0:"<<ic0<<"\n";
      // Rcout<<"-IC1:"<<this->ic<<"\n";
      if (ic0 < this->ic) {
        this->beta = beta0;
        this->lam = lam0;
        this->ic = ic0;
        beta1 = beta0;
      }
      lam1 = lam0;
      //coef2 = 1
      if (lam1 < coef2*log(exp(1)*d/s_0)*delta_tbar/n) break;
    }
    this->group_size = group_support_size(this->beta, gindex, gsize, m);
    this->size = (this->beta.array() != 0).count();
    this->A_out = support_set(this->beta, p, this->size)+Eigen::VectorXi::Ones(this->size);
  }
};

class DSIHTLm : public Algorithm {
public:
  DSIHTLm(Data &data) : Algorithm(data){};
};

#endif //SRC_ALGORITHM_H
