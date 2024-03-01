#include <algorithm>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
using namespace Rcpp;

double Delta(double s_0, int m, int d);
Eigen::VectorXd tau_eta(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXd &gadient, double eta, Eigen::VectorXi &gindex, Eigen::VectorXi &gsize, double lambda, double s_0, int m, int p, int max_iter);
Eigen::VectorXd tau(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXi &gindex, Eigen::VectorXi &gsize, double lambda, double s_0, int m, int p);
int group_support_size(Eigen::VectorXd &beta, Eigen::VectorXi &gindex, Eigen::VectorXi &gsize, int m);
Eigen::VectorXi support_set(Eigen::VectorXd &beta, int p, int size);
double IC(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXi &gindex, Eigen::VectorXi &gsize, double s_0, int n, int m, int p, int d, double delta_t, double ic_coef);
Eigen::VectorXd least_square(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, int p);
