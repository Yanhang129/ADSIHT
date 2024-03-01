#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]
#include <iostream>
#include "Data.h"
#include "Algorithm.h"
#include "utilities.h"
#include "Metric.h"
using namespace Eigen;
using namespace std;

List sequential_path(Data &data, Algorithm *algorithm, Metric *metric, Eigen::VectorXd sequence, double rho, double ic_coef, double coef1, double coef2)
{
  int p = data.p;
  int m = data.g_num;
  int sequence_size = sequence.size();
  List A_out(sequence_size);
  Eigen::VectorXd ic_sequence = Eigen::VectorXd::Zero(sequence_size);
  Eigen::VectorXd lam_sequence = Eigen::VectorXd::Zero(sequence_size);
  Eigen::MatrixXd beta_matrix= Eigen::MatrixXd::Zero(p, sequence_size);
  Eigen::VectorXd intercept_sequence = Eigen::VectorXd::Zero(sequence_size);
  Eigen::VectorXi support_sequence = Eigen::VectorXi::Zero(sequence_size);
  algorithm->update_rho(rho);
  for (int i = 0; i < sequence_size; i++) {
    algorithm->fit1(sequence(i), ic_coef, coef1, coef2);
    beta_matrix.col(i) = algorithm->get_beta();
    lam_sequence(i) = algorithm->get_lambda();
    support_sequence(i) = algorithm->get_support_size();
    A_out[i] = algorithm->get_A_out();
    ic_sequence(i) = metric->ic(algorithm, data);
  }
  beta_matrix = beta_matrix.array().colwise()/data.x_norm.array()*sqrt(data.n);
  intercept_sequence = data.y_mean*Eigen::VectorXd::Ones(sequence_size) - beta_matrix.transpose()*data.x_mean;
  int min_ic = 0;
  ic_sequence.minCoeff(&min_ic);

  return List::create(Named("beta") = beta_matrix,
                      Named("intercept") = intercept_sequence,
                      Named("lambda") = lam_sequence,
                      Named("support_size") = support_sequence,
                      Named("A_out") = A_out,
                      Named("ic") = ic_sequence);
}

List sequential_path_eta(Data &data, Algorithm *algorithm, Metric *metric, Eigen::VectorXd sequence, double rho, double ic_coef, double coef1, double coef2, double eta, int max_iter)
{
  int p = data.p;
  int m = data.g_num;
  int sequence_size = sequence.size();
  List A_out(sequence_size);
  Eigen::VectorXd ic_sequence = Eigen::VectorXd::Zero(sequence_size);
  Eigen::VectorXd lam_sequence = Eigen::VectorXd::Zero(sequence_size);
  Eigen::MatrixXd beta_matrix= Eigen::MatrixXd::Zero(p, sequence_size);
  Eigen::VectorXd intercept_sequence = Eigen::VectorXd::Zero(sequence_size);
  Eigen::VectorXi support_sequence = Eigen::VectorXi::Zero(sequence_size);
  algorithm->update_rho(rho);
  for (int i = 0; i < sequence_size; i++) {
    algorithm->fit_eta(sequence(i), ic_coef, coef1, coef2, eta, max_iter);
    beta_matrix.col(i) = algorithm->get_beta();
    lam_sequence(i) = algorithm->get_lambda();
    support_sequence(i) = algorithm->get_support_size();
    A_out[i] = algorithm->get_A_out();
    ic_sequence(i) = metric->ic(algorithm, data);
  }
  beta_matrix = beta_matrix.array().colwise()/data.x_norm.array()*sqrt(data.n);
  intercept_sequence = data.y_mean*Eigen::VectorXd::Ones(sequence_size) - beta_matrix.transpose()*data.x_mean;
  int min_ic = 0;
  ic_sequence.minCoeff(&min_ic);

  return List::create(Named("beta") = beta_matrix,
                      Named("intercept") = intercept_sequence,
                      Named("lambda") = lam_sequence,
                      Named("support_size") = support_sequence,
                      Named("A_out") = A_out,
                      Named("ic") = ic_sequence);
}
