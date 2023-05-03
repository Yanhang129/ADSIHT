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
  mcc <- mccr(y, pre)
  result <- c(tp+fp, tp/(tp+fn), fp/(tn+fp), mcc)
  return(result)
}

## Generate coefficients
gen.coef <- function(s, k, Tn, J, i, signal){
  set.seed(i)
  ind_group <- sort(sample(1:J, Tn))
  coef <- rep(0, k*J)
  # temp1 <- 2*signal*rbinom(s*Tn, 1, 0.5)-1
  temp1 <- rnorm(s*Tn, 0, 1)
  for (i in 1:Tn) {
    temp2 <- rep(0, k)
    temp2[sample(1:k, s)] <- temp1[((i-1)*s+1):(i*s)]
    coef[((ind_group[i]-1)*k+1):(ind_group[i]*k)] <- temp2
  }
  return(coef)
}

# sgl <- function(x, y, group){
#   fit <- sparsegl(x, y, group)
#   temp <- estimate_risk(fit, type = "EBIC", approx_df = TRUE)
#   ind <- which.min(temp[, 3])
#   beta <- fit$beta[, ind]
#   inter <- fit$b0[ind]
#   return(list(beta = beta, intercept = inter))
# }

sgl <- function(x, y, group){
  fit <- cv.sparsegl(x, y, group, nfolds = 5)
  coeff <- coef(fit, s = "lambda.min")
  beta <- coeff[-1]
  inter <- coeff[1]
  return(list(beta = beta, intercept = inter))
}

grp <- function(x, y, group, penalty){
  fit <- grpreg(x, y, group = group, penalty = penalty, dfmax = 100)
  lam <- grpreg::select(fit, crit = "EBIC")$lambda
  return(list(model = fit, lam = lam))
}

gB <- function(x, y, group){
  fit <- gBridge(x, y, group = group)
  lam <- grpreg::select(fit, crit = "EBIC")$lambda
  return(list(model = fit, lam = lam))
}

detect.group <- function(beta, J, k){
  size <- 0
  for (i in 1:J) {
    temp <- beta[((i-1)*k+1):(i*k)]
    size <- size + ifelse(any(temp != 0 ), 1, 0)
  }
  return(size)
}

i <- 20
s <- 10
k <- 20
J <- 150
Tn <- 5
p <- J*k
snr <- 5


res <- list()
n_seq <- seq(300, 1000, 100)
for (n in n_seq) {
  result <- sapply(1:100, function(i){
    print(i)
    set.seed(i)

    beta <- gen.coef(s, k, Tn, J, i, 1)
    R <- matrix(rnorm(p*n, 0, 1), n, p)
    Sigma <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma[i, j] <- 0.5^(abs(i-j))
        # Sigma[i, j] <- ifelse(i==j, 1, 0.5)
      }
    }
    x <- rmvn(n, rep(0, p), Sigma)
    pe <- x %*% beta
    sigma <- sqrt(var(pe)/snr)
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
    lasso1 <- c(f(true, lasso), detect.group(beta1, J, k), sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)
    round(lasso1, 3)

    t <- as.numeric(system.time(fit <- gB(x, y, group = group))[3])
    gbridge <- rep(0, J*k)
    gr <- predict(fit$model, type = "vars", lambda = fit$lam)
    coef <- predict(fit$model, type = "coef", lambda = fit$lam)
    gbridge[as.numeric(gr)] <- 1
    beta1 <- coef[-1]
    intercept <- coef[1]
    gbridge1 <- c(f(true, gbridge), detect.group(beta1, J, k), sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)
    round(gbridge1, 3)

    t <- as.numeric(system.time(fit <- grp(x, y, group = group, penalty = "gel"))[3])
    explasso <- rep(0, J*k)
    gr <- predict(fit$model, type = "vars", lambda = fit$lam)
    coef <- predict(fit$model, type = "coef", lambda = fit$lam)
    explasso[as.numeric(gr)] <- 1
    beta1 <- coef[-1]
    intercept <- coef[1]
    explasso1 <- c(f(true, explasso), detect.group(beta1, J, k), sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)
    round(explasso1, 3)

    t <- as.numeric(system.time(fit <- grp(x, y, group = group, penalty = "cMCP"))[3])
    MCP <- rep(0, J*k)
    gr <- predict(fit$model, type = "vars", lambda = fit$lam)
    coef <- predict(fit$model, type = "coef", lambda = fit$lam)
    MCP[as.numeric(gr)] <- 1
    beta1 <- coef[-1]
    intercept <- coef[1]
    MCP1 <- c(f(true, MCP), detect.group(beta1, J, k), sqrt(sum(sum((beta1-beta)^2)+intercept^2)), t)
    round(MCP1, 3)

    t <- as.numeric(system.time(fit <- DSIHT(x, y, group = group, L = 5, rho = 0.9, method = "fast2", ic.type = "ebic", ic.coef2 = 0.5))[3])
    l <- which.min(fit$ic)
    coef <- c(fit$intercept[l], fit$beta[, l])
    best.group <- which(fit$beta[, l] != 0)
    dsiht <- rep(0, p)
    dsiht[best.group] <- 1
    dsiht1 <- c(f(true, dsiht), detect.group(coef[-1], J, k), sqrt(sum(sum((coef[-1]-beta)^2)+coef[1]^2)), t)
    round(dsiht1, 3)

    return(c(lasso1, gbridge1, explasso1, MCP1, dsiht1))
  })
  res[[n/100-2]] <- matrix(round(apply(result, 1, mean), 2), ncol = 5)
}

matrix(round(apply(result, 1, mean), 2), ncol = 5)
matrix(round(apply(result, 1, sd), 2), ncol = 5)
round(result[1:5, ], 3)
round(result[8:14, ], 3)

save(res, file = "C:/Users/Yanhang/Desktop/code/result/res_n_normal.rda")
