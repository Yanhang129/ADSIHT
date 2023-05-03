rm(list = ls())
library(DSIHT)
library(grpreg)
library(sparsegl)
library(stringr)
library(mccr)
library(mvnfast)

## Evaluation criteria
f <- function(y, pre){
  tp <- sum(pre==1 & y==1)
  fp <- sum(pre==1 & y!=1)
  fn <- sum(pre!=1 & y==1)
  tn <- sum(pre!=1 & y!=1)
  result <- c(tp+fp, tp/(tp+fn), fp/(tn+fp))
  return(result)
}

## Generate coefficients
gen.coef <- function(s, k, Tn, J, i, signal){
  set.seed(i)
  ind_group <- sort(sample(1:J, Tn))
  coef <- rep(0, k*J)
  # temp1 <- 2*rbinom(s*Tn, 1, 0.5)-1
  temp1 <- runif(s*Tn, signal, 100*signal)*(2*rbinom(s*Tn, 1, 0.5)-1)
  for (i in 1:Tn) {
    temp2 <- rep(0, k)
    temp2[sample(1:k, s)] <- temp1[((i-1)*s+1):(i*s)]
    coef[((ind_group[i]-1)*k+1):(ind_group[i]*k)] <- temp2
  }
  return(coef)
}

sgl <- function(x, y, group){
  fit <- sparsegl(x, y, group, eps = 1e-3, asparse = 0.1)
  temp <- estimate_risk(fit, type = "EBIC1", approx_df = TRUE)
  ind <- which.min(temp[, 3])
  beta <- fit$beta[, ind]
  inter <- fit$b0[ind]
  return(list(beta = beta, intercept = inter))
}

grp <- function(x, y, group, penalty){
  fit <- grpreg(x, y, group = group, penalty = penalty, lambda.min = 0.01, eps = 1e-8)
  lam <- grpreg::select(fit, crit = "EBIC")$lambda
  return(list(model = fit, lam = lam))
}



##############################################################################
## Simulation 1:varying sparsity in group
##############################################################################

i <- 4
n <- 500
s <- 5
k <- 10
J <- 500
Tn <- 8
snr <- 10
p <- J*k
sigma <- 1


result <- sapply(1:100, function(i){
  print(i)
  set.seed(i)

  beta <- gen.coef(s, k, Tn, J, i, sigma*sqrt(2*log(p)/n))
  R <- matrix(rnorm(p*n, 0, 1), n, p)
  Sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] <- 0.3^(abs(i-j))
    }
  }
  x <- rmvn(n, rep(0, p), Sigma)
  pe <- x %*% beta
  y <- pe+rnorm(n, 0, sigma)
  true <- rep(0, p)
  true[which(beta != 0)] <- 1
  group <- rep(1:J, each = k)
  t <- as.numeric(system.time(fit <- sgl(x, y, group = group))[3])
  lasso <- rep(0, J*k)
  gr <- which(fit$beta != 0)
  lasso[as.numeric(gr)] <- 1
  beta1 <- fit$beta
  intercept <- fit$intercept
  lasso1 <- c(f(true, lasso), 10*sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)


  t <- as.numeric(system.time(fit <- grp(x, y, group = group, penalty = "gel"))[3])
  explasso <- rep(0, J*k)
  gr <- predict(fit$model, type = "vars", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  explasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  explasso1 <- c(f(true, explasso), 10*sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)

  t <- as.numeric(system.time(fit <- grp(x, y, group = group, penalty = "cMCP"))[3])
  MCP <- rep(0, J*k)
  gr <- predict(fit$model, type = "vars", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1 <- c(f(true, MCP), 10*sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)
#
#   t <- as.numeric(system.time(fit <- DSIHT(x, y, group = group, L = 5, rho = 0.95))[3])
#   l <- which.min(fit$ic)
#   coef <- c(fit$intercept[l], fit$beta[, l])
#   best.group <- which(fit$beta[, l] != 0)
#   dsiht <- rep(0, p)
#   dsiht[best.group] <- 1
#   dsiht1 <- c(f(true, dsiht), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
#   round(dsiht1, 3)

  t <- as.numeric(system.time(fit <- DSIHT(x, y, group = group, L = 5, rho = 0.95, method = "fast", ic.type = "ebic"))[3])
  l <- which.min(fit$ic)
  coef <- c(fit$intercept[l], fit$beta[, l])
  best.group <- which(fit$beta[, l] != 0)
  dsiht <- rep(0, p)
  dsiht[best.group] <- 1
  dsiht2 <- c(f(true, dsiht), 10*sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)), t)
  # round(dsiht2, 3)

  return(c(lasso1, explasso1, MCP1, dsiht2))
      # return(dsiht2)
})

matrix(round(apply(result, 1, mean), 2), ncol = 4)
matrix(round(apply(result, 1, sd), 2), ncol = 4)

# # s0 = 8
# matrix(round(apply(result[, -c(14, 19, 34, 42, 100)], 1, mean), 3), ncol = 4)
# matrix(round(apply(result[, -c(14, 19, 34, 42, 100)], 1, sd), 3), ncol = 4)


round(result[16:20, ], 3)
round(result[6:10, ], 3)

##############################################3
result <- sapply(1:100, function(i){
  print(i)
  set.seed(i)

  beta <- gen.coef(s, k, Tn, J, i, sigma*sqrt(2*log(p)/n))
  R <- matrix(rnorm(p*n, 0, 1), n, p)
  Sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] <- 0.3^(abs(i-j))
    }
  }
  x <- rmvn(n, rep(0, p), Sigma)
  pe <- x %*% beta
  y <- pe+rnorm(n, 0, sigma)
  true <- rep(0, p)
  true[which(beta != 0)] <- 1
  group <- rep(1:J, each = k)
  t <- as.numeric(system.time(fit <- sgl(x, y, group = group))[3])
  lasso <- rep(0, J*k)
  gr <- which(fit$beta != 0)
  lasso[as.numeric(gr)] <- 1
  beta1 <- fit$beta
  intercept <- fit$intercept
  lasso1 <- c(f(true, lasso), 10*sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)

  return(lasso1)
})

ind <- which(result[1, ] == 0)
round(apply(result[, -ind], 1, mean), 2)
round(apply(result[, -ind], 1, sd), 2)


result <- sapply(1:100, function(i){
  print(i)
  set.seed(i)

  beta <- gen.coef(s, k, Tn, J, i, sigma*sqrt(2*log(p)/n))
  R <- matrix(rnorm(p*n, 0, 1), n, p)
  Sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] <- 0.3^(abs(i-j))
    }
  }
  x <- rmvn(n, rep(0, p), Sigma)
  pe <- x %*% beta
  y <- pe+rnorm(n, 0, sigma)
  true <- rep(0, p)
  true[which(beta != 0)] <- 1
  group <- rep(1:J, each = k)
  t <- as.numeric(system.time(fit <- sgl(x, y, group = group))[3])
  lasso <- rep(0, J*k)
  gr <- which(fit$beta != 0)
  lasso[as.numeric(gr)] <- 1
  beta1 <- fit$beta
  intercept <- fit$intercept
  lasso1 <- c(f(true, lasso), 10*sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)


  t <- as.numeric(system.time(fit <- grp(x, y, group = group, penalty = "gel"))[3])
  explasso <- rep(0, J*k)
  gr <- predict(fit$model, type = "vars", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  explasso[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  explasso1 <- c(f(true, explasso), 10*sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)

  t <- as.numeric(system.time(fit <- grp(x, y, group = group, penalty = "cMCP"))[3])
  MCP <- rep(0, J*k)
  gr <- predict(fit$model, type = "vars", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1 <- c(f(true, MCP), 10*sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)
  #
  #   t <- as.numeric(system.time(fit <- DSIHT(x, y, group = group, L = 5, rho = 0.95))[3])
  #   l <- which.min(fit$ic)
  #   coef <- c(fit$intercept[l], fit$beta[, l])
  #   best.group <- which(fit$beta[, l] != 0)
  #   dsiht <- rep(0, p)
  #   dsiht[best.group] <- 1
  #   dsiht1 <- c(f(true, dsiht), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)/sum(beta^2)), t)
  #   round(dsiht1, 3)

  t <- as.numeric(system.time(fit <- DSIHT(x, y, group = group, L = 5, rho = 0.95, method = "fast", ic.type = "ebic"))[3])
  l <- which.min(fit$ic)
  coef <- c(fit$intercept[l], fit$beta[, l])
  best.group <- which(fit$beta[, l] != 0)
  dsiht <- rep(0, p)
  dsiht[best.group] <- 1
  dsiht2 <- c(f(true, dsiht), 10*sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)), t)
  # round(dsiht2, 3)

  return(c(lasso1, explasso1, MCP1, dsiht2))
  # return(dsiht2)
})

matrix(round(apply(result, 1, mean), 2), ncol = 4)
matrix(round(apply(result, 1, sd), 2), ncol = 4)

# # s0 = 8
# matrix(round(apply(result[, -c(14, 19, 34, 42, 100)], 1, mean), 3), ncol = 4)
# matrix(round(apply(result[, -c(14, 19, 34, 42, 100)], 1, sd), 3), ncol = 4)


round(result[16:20, ], 3)
round(result[6:10, ], 3)



result <- sapply(1:50, function(i){
  print(i)
  set.seed(i)

  beta <- gen.coef(s, k, Tn, J, i, sigma*sqrt(2*log(p)/n))
  R <- matrix(rnorm(p*n, 0, 1), n, p)
  Sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] <- 0.5^(abs(i-j))
    }
  }
  x <- rmvn(n, rep(0, p), Sigma)
  pe <- x %*% beta
  y <- pe+rnorm(n, 0, sigma)
  true <- rep(0, p)
  true[which(beta != 0)] <- 1
  group <- rep(1:J, each = k)
  t <- as.numeric(system.time(fit <- grp(x, y, group = group, penalty = "cMCP"))[3])
  MCP <- rep(0, J*k)
  gr <- predict(fit$model, type = "vars", lambda = fit$lam)
  coef <- predict(fit$model, type = "coef", lambda = fit$lam)
  MCP[as.numeric(gr)] <- 1
  beta1 <- coef[-1]
  intercept <- coef[1]
  MCP1 <- c(f(true, MCP), 10*sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)

  t <- as.numeric(system.time(fit <- DSIHT(x, y, group = group, L = 5, rho = 0.9, method = "fast2", ic.type = "ebic", ic.coef2 = 0.1))[3])
  l <- which.min(fit$ic)
  coef <- c(fit$intercept[l], fit$beta[, l])
  best.group <- which(fit$beta[, l] != 0)
  dsiht <- rep(0, p)
  dsiht[best.group] <- 1
  dsiht2 <- c(f(true, dsiht), 10*sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)), t)
  # round(dsiht2, 3)

  return(c(MCP1, dsiht2))
  # return(dsiht2)
})

matrix(round(apply(result, 1, mean), 2), ncol = 2)
matrix(round(apply(result, 1, sd), 2), ncol = 2)

# # s0 = 8
# matrix(round(apply(result[, -c(14, 19, 34, 42, 100)], 1, mean), 3), ncol = 4)
# matrix(round(apply(result[, -c(14, 19, 34, 42, 100)], 1, sd), 3), ncol = 4)


round(result[16:20, ], 3)
round(result[6:10, ], 3)
