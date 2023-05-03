#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]
#include <iostream>
#include "Data.h"
#include "Algorithm.h"
#include "path.h"
#include "utilities.h"
#include "Metric.h"
#include <vector>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List DSIHT_Cpp(Eigen::MatrixXd &x, Eigen::VectorXd &y, Eigen::VectorXd &weight, int ic_type,
               double ic_scale, Eigen::VectorXd &sequence, double kappa, Eigen::VectorXi &g_index, double ic_coef)
{
  Data data(x, y, weight, g_index);
  data.add_weight();
  Algorithm *algorithm = new DSIHTLm(data);
  Metric *metric = new LmMetric(ic_type, ic_scale);
  List result = sequential_path(data, algorithm, metric, sequence, kappa, ic_coef);
  algorithm -> ~Algorithm();
  metric -> ~Metric();
  return result;
}
