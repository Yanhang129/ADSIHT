#ifndef SRC_PATH_H
#define SRC_PATH_H
#include <RcppEigen.h>
#include <Rcpp.h>
#include "Data.h"
#include "Algorithm.h"
#include "Metric.h"
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

List sequential_path(Data &data, Algorithm *algorithm, Metric *metric, Eigen::VectorXd sequence, double rho, double ic_coef, double coef1, double coef2);
List sequential_path_eta(Data &data, Algorithm *algorithm, Metric *metric, Eigen::VectorXd sequence, double rho, double ic_coef, double coef1, double coef2, double eta, int max_iter);

#endif //SRC_PATH_H
