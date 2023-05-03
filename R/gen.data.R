#' @title Generate simulated data
#'
#' @description Generate simulated data for sparse group linear model.
#'
#' @param n The number of observations.
#' @param m The number of groups of interest.
#' @param d The group size of each group. Only even group structure is allowed here.
#' @param s The number of important groups in the underlying regression model.
#' @param s0 The number of important variables in each important group.
#' @param rho A parameter used to characterize the pairwise correlation in
#' predictors. Default is \code{0.5}..
#' @param cor.type The structure of correlation.
#' \code{cor.type = 1} denotes the independence structure,
#' where the covariance matrix has \eqn{(i,j)} entry equals \eqn{I(i \neq j)}.
#' \code{cor.type = 2} denotes the exponential structure,
#' where the covariance matrix has \eqn{(i,j)} entry equals \eqn{rho^{|i-j|}}.
#' \code{cor.type = 3} denotes the constant structure,
#' where the non-diagnoal entries of covariance
#' matrix are \eqn{rho} and diagonal entries are 1.\
#' @param beta.type The structure of coefficients.
#' \code{beta.type = 1} denotes the homogenous setup,
#' where each entry has the same magnitude.
#' \code{beta.type = 2} denotes the heterogeneous structure,
#' where the coefficients are drawn from a normal distribution.
#' @param sigma1 The value controlling the strength of the gaussian noise. A large value implies strong noise. Default \code{sigma1 = 1}.
#' @param sigma2 The value controlling the strength of the coefficients. A large value implies large coefficients. Default \code{sigma2 = 1}.
#' @param seed random seed. Default: \code{seed = 1}.
#'
#' @return A \code{list} object comprising:
#' \item{x}{Design matrix of predictors.}
#' \item{y}{Response variable.}
#' \item{beta}{The coefficients used in the underlying regression model.}
#' \item{group}{The group index of each variable.}
#' \item{true.group}{The important groups in the sparse group linear model.}
#' \item{true.variable}{The important variables in the sparse group linear model.}
#'
#' @author Yanhang Zhang, Zhifan Li, Jianxin Yin.
#'
#' @export
#'
#' @examples
#'
#' # Generate simulated data
#' n <- 200
#' m <- 100
#' d <- 10
#' s <- 5
#' s0 <- 5
#' data <- gen.data(n, m, d, s, s0)
#' str(data)

gen.data <- function(n,
                     m,
                     d,
                     s,
                     s0,
                     cor.type = 1,
                     beta.type = 1,
                     rho = 0.5,
                     sigma1 = 1,
                     sigma2 = 1,
                     seed = 1) {
  set.seed(seed)
  group <- rep(1:m, each = d)
  p <- m*d

  if (cor.type == 1) {
    Sigma <- diag(p)
  } else if (cor.type == 2) {
    Sigma <- matrix(0, p, p)
    Sigma <- rho ^ (abs(row(Sigma) - col(Sigma)))
  } else if (cor.type == 3) {
    Sigma <- matrix(rho, p, p)
    diag(Sigma) <- 1
  }
  if (cor.type == 1) {
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  } else {
    x <- mvnfast::rmvn(n, rep(0, p), Sigma)
  }

  ind_group <- sort(sample(1:m, s))
  coef <- rep(0, p)
  if (beta.type == 1) {
    temp1 <- sigma2*(2*rbinom(s*s0, 1, 0.5)-1)
  } else {
    temp1 <- rnorm(s*s0, 0, sigma2)
  }
  for (i in 1:s) {
    temp2 <- rep(0, d)
    temp2[sample(1:d, s0)] <- temp1[((i-1)*s0+1):(i*s0)]
    coef[((ind_group[i]-1)*d+1):(ind_group[i]*d)] <- temp2
  }

  y <- x %*% coef + rnorm(n, 0, sigma1)
  set.seed(NULL)
  colnames(x) <- paste0("x", 1:p)
  return(list(
    x = x,
    y = y,
    beta = coef,
    group = group,
    true.group = ind_group,
    true.variable = which(coef != 0)
  ))
}
