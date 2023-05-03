#ifndef SRC_DATA_H
#define SRC_DATA_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include "normalize.h"
// [[Rcpp::plugins("cpp11")]]
using namespace std;

class Data {
public:
  Eigen::MatrixXd x;
  Eigen::VectorXd y;
  Eigen::VectorXd weight;
  Eigen::VectorXd x_mean;
  Eigen::VectorXd x_norm;
  double y_mean;
  int n;
  int p;
  int g_num;
  Eigen::VectorXi g_index;
  Eigen::VectorXi g_size;

  Data() = default;

  Data(Eigen::MatrixXd& x, Eigen::VectorXd& y, Eigen::VectorXd& weight, Eigen::VectorXi& g_index) {
    this->x = x;
    this->y = y;
    this->n = x.rows();
    this->p = x.cols();
    this->weight = weight;
    this->x_mean = Eigen::VectorXd::Zero(this->p);
    this->x_norm = Eigen::VectorXd::Zero(this->p);
    this->weight = weight;
    this->g_index = g_index;
    this->g_num = (g_index).size();
    if (g_num > 1) {
      Eigen::VectorXi temp = Eigen::VectorXi::Zero(g_num);
      temp.head(g_num-1) = g_index.tail(g_num-1);
      temp(g_num-1) = this->p;
      this->g_size =  temp-g_index;
    }
    this->normalize();
  };

  void add_weight() {
    for(int i=0;i<this->n;i++){
      this->x.row(i) = this->x.row(i)*sqrt(this->weight(i));
      this->y(i) = this->y(i)*sqrt(this->weight(i));
    }
  };

  void normalize() {
    Normalize(this->x, this->y, this->weight, this->x_mean, this->y_mean, this->x_norm);
  };

  Eigen::VectorXi get_g_index() {
    return this->g_index;
  };

  Eigen::VectorXi get_g_size() {
    return this->g_size;
  };
};
#endif //SRC_DATA_H
