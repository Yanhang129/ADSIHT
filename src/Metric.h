#ifndef SRC_METRICS_H
#define SRC_METRICS_H

#include "Data.h"
#include "Algorithm.h"
#include <vector>
#include <random>
#include <algorithm>
using namespace std;
// [[Rcpp::plugins("cpp11")]]

class Metric {
public:
  int ic_type;
  double ic_coef;
  Metric() = default;

  Metric(int ic_type, double ic_coef = 1.0) {
    this->ic_type = ic_type;
    this->ic_coef = ic_coef;
  };

  virtual double loss(Algorithm *algorithm, Data &data) = 0;

  virtual double ic(Algorithm *algorithm, Data &data) = 0;
};

class LmMetric : public Metric {
public:

  LmMetric(int ic_type, double ic_coef) : Metric(ic_type, ic_coef) {};

  double loss(Algorithm *algorithm, Data &data) {
    return (data.y - data.x * algorithm->get_beta()).squaredNorm() / (data.n);
  }

  double ic(Algorithm *algorithm, Data &data) {
    if (ic_type == 0) {
      return this->loss(algorithm, data);
    }  else if (ic_type == 1) {
      return  double(data.n) * log(this->loss(algorithm, data)) + this->ic_coef*(algorithm->get_group_support_size()*log(exp(1)*double(data.g_num)/algorithm->get_group_support_size())+
                            algorithm->get_support_size()*log(exp(1)*double(data.p)/double(data.g_num)));

    }  else return 0;
  }
};

#endif //SRC_METRICS_H
