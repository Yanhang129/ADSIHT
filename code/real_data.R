rm(list = ls()); gc(reset = TRUE)
library(DSIHT)
library(grpreg)
library(splines)
library(sparsegl)
library(Ball)
library(abess)
library(SIS)
library(dplyr)
library(qqplotr)
load("C://Users//Yanhang//Desktop//code//trim32.rda")

## Data pre-processing
t <-500
# z <- apply(x, 2, function(x) bcor(x, y))
z <- abs(apply(x, 2, function(x) cor(x, y)))
id <- which(z >= sort(z, decreasing = T)[t])
temp <- x[, id]
X <- ns(temp[, 1], df = 6)
for (i in 2:ncol(temp)) {
  X <- cbind(X, ns(temp[, i], df = 6))
}
group = rep(1:ncol(temp), each = 6)
sgl <- function(x, y, group){
  fit <- cv.sparsegl(x, y, group, nfolds = 5)
  coeff <- coef(fit, s = "lambda.min")
  beta <- coeff[-1]
  inter <- coeff[1]
  return(list(beta = beta, intercept = inter))
}

grp <- function(x, y, group, penalty){
  fit <- grpreg(x, y, group = group, penalty = penalty)
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

result <- sapply(101:200, function(i) {
  set.seed(i)
  print(i)
  ind <- sample(1:120, 100)
  fit1 <- sgl(X[ind, ], y[ind], group = group)
  l1 <- length(which(fit1$beta != 0))
  g1 <- detect.group(fit1$beta, 500, 6)
  res1 <- y[-ind] - X[-ind, ]%*%fit1$beta-fit1$intercept
  sum(res1^2)

  fit2 <- gB(X[ind, ], y[ind], group = group)
  coef <- predict(fit2$model, type = "coef", lambda = fit2$lam)
  l2 <- predict(fit2$model, type = "nvars", lambda = fit2$lam)
  g2 <- predict(fit2$model, type = "ngroups", lambda = fit2$lam)
  res2 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res2^2)

  fit3 <- grp(X[ind, ], y[ind], group = group, penalty = "gel")
  coef <- predict(fit3$model, type = "coef", lambda = fit3$lam)
  l3 <- predict(fit3$model, type = "nvars", lambda = fit3$lam)
  g3 <- predict(fit3$model, type = "ngroups", lambda = fit3$lam)
  res3 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res3^2)

  fit4 <- grp(X[ind, ], y[ind], group = group, penalty = "cMCP")
  coef <- predict(fit4$model, type = "coef", lambda = fit4$lam)
  l4 <- predict(fit4$model, type = "nvars", lambda = fit4$lam)
  g4 <- predict(fit4$model, type = "ngroups", lambda = fit4$lam)
  res4 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res4^2)

  fit5 <- DSIHT(X[ind, ], y[ind], group = group, L = 5, kappa = 0.95, ic.coef = 1.6, ic.scale = 0.3)
  l <- which.min(fit5$ic)
  coef <- c(fit5$intercept[l], fit5$beta[, l])
  l5 <- length(which(fit5$beta[, l] != 0))
  g5 <- detect.group(fit5$beta[, l], 500, 6)
  res5 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res5^2)

  return(c(sum(res1^2),
           sum(res2^2),
           sum(res3^2),
           sum(res4^2),
           sum(res5^2),
           g1, g2, g3, g4, g5,
           l1, l2, l3, l4, l5))})
# save(result, file = "C:/Users/Yanhang/Desktop/code/result/res_realdata.rda")

# result2 <- result
# result[1:5, ] <- sqrt(result[1:5, ])
round(matrix(apply(result, 1, mean), ncol = 3), 2)
round(matrix(apply(result, 1, sd), ncol = 3), 2)
library(ggplot2)
dat <- data.frame(rep(c("SGLasso", "GBridge", "GEL", "cMCP", "ADSIHT"), each = 100), as.numeric(as.vector(t(result[1:5, ]))))
colnames(dat) <- c("Methods", "V2")
dat[, 1] <- factor(dat[, 1], levels = c("SGLasso", "GBridge", "GEL", "cMCP", "ADSIHT"))
dat[401:500, 2] <- dat[401:500, 2] - 0.01
# dat <- dat[-which(dat$V2 > 1), ]
p <- ggplot(dat, aes(x=Methods, y=log(V2), fill = Methods)) +
  geom_boxplot()+ylab('The logrithm of prediction error')+theme_bw()+
  theme(legend.position = "bottom", legend.key.size = unit(25, "pt"), legend.box.spacing = unit(0, 'pt'))
dat %>% group_by(Methods) %>% summarise(mean(V2))
dat %>% group_by(Methods) %>% summarise(median(V2))
# ggsave(file = "C:\\Users\\Yanhang\\Desktop\\boxplot.pdf", dpi = 600, p, height = 5, width = 8)


fit1 <- sgl(X, y, group = group)
res1 <- y - X%*%fit1$beta-fit1$intercept
sort(colnames(temp)[unique(which(fit1$beta != 0) %/% 6)])

fit2 <- gB(X, y, group = group)
coef <- predict(fit2$model, type = "coef", lambda = fit2$lam)
colnames(temp)[predict(fit2$model, type = "groups", lambda = fit2$lam)]
res2 <- y - X%*%coef[-1]-coef[1]

fit3 <- grp(X, y, group = group, penalty = "gel")
coef <- predict(fit3$model, type = "coef", lambda = fit3$lam)
sort(colnames(temp)[predict(fit3$model, type = "groups", lambda = fit3$lam)])
res3 <- y - X%*%coef[-1]-coef[1]


fit4 <- grp(X, y, group = group, penalty = "cMCP")
coef <- predict(fit4$model, type = "coef", lambda = fit4$lam)
sort(colnames(temp)[predict(fit4$model, type = "groups", lambda = fit4$lam)])
res4 <- y - X%*%coef[-1]-coef[1]
summary(lm(y~X[, predict(fit4$model, type = "groups", lambda = fit4$lam)]))

fit5 <- DSIHT(X, y, group = group, L = 5, rho = 0.95, ic.coef2 = 1.6, ic.coef = 0.3)
l <- which.min(fit5$ic)
coef <- c(fit5$intercept[l], fit5$beta[, l])
res5 <- y - X%*%coef[-1]-coef[1]
summary(lm(y~X[, unlist(fit5$A_out[l])]))
sort(colnames(temp)[unique((unlist(fit5$A_out[l])) %/% 6)])

dat <- data.frame(rep(c("SGLasso", "GBridge", "GEL", "cMCP", "ADSIHT"), each = 120), c(res1, res2, res3, res4, res5))
colnames(dat) <- c("Method", "V2")
dat[, 1] <- factor(dat[, 1], levels = c("SGLasso", "GBridge", "GEL", "cMCP", "ADSIHT"))
# tibble(y = res5) %>%
#   ggplot(aes(sample = y)) +
#   geom_qq() + geom_qq_line()

p <- ggplot(dat, aes(sample = V2, colour = Method)) +
  geom_qq(alpha=0.7) + geom_qq_line()+
  facet_wrap(vars(Method)) +
  labs(x = "Theoretical Quantiles ", y = "Sample Quantiles")+theme_bw()+theme(legend.position = c(0.8,0.2))
ggsave(file = "C:\\Users\\Yanhang\\Desktop\\qq.pdf", dpi = 600, p, height = 5.5, width = 8)


M <- 100
result2 <- lapply(1:M, function(i) {
  set.seed(i)
  print(i)
  ind <- sample(1:120, 60)
  fit <- DSIHT(X[ind, ], y[ind], group = group, L = 10, rho = 0.95, ic.coef2 = 1.6, ic.coef = 0.3)
  res.list <- lapply(fit$A_out, function(x) {
    x <- x%/%6+1
    return(unique(x))
  })
  return(res.list)}
)

res.s1 <- sort(table(unlist(sapply(result2, function(x) x[1]))))
sort(table(unlist(sapply(result2, function(x) x[9]))))


res <- lapply(result2, function(x) {
  x <- x%/%6+1
  return(unique(x))
})
sort(table(unlist(res)))
# return(c(sum(res5^2), sum(res6^2),l5, l6))})
matrix(apply(result2, 1, mean), ncol = 2)
matrix(apply(result2[, 101:200], 1, mean), ncol = 2)
matrix(apply(result2[, 101:200], 1, mean), ncol = 2)


s1 <- seq(1, 1.5, 0.1)
s2 <- seq(0.1, 0.5, 0.1)
ab <- matrix(0, length(s1), length(s2))
for (a in 1:length(s1)) {
  print(a)
  for (b in 1:length(s2)) {
    result2 <- sapply(101:200, function(i) {
      set.seed(i)
      # print(i)
      ind <- sample(1:120, 100)
      fit5 <- DSIHT(X[ind, ], y[ind], group = group, L = 5, rho = 0.9, ic.coef2 = s1[a], ic.coef = s2[b])
      l <- which.min(fit5$ic)
      coef <- c(fit5$intercept[l], fit5$beta[, l])
      l5 <- length(which(fit5$beta[, l] != 0))
      res5 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
      sum(res5^2)
      return(c(sum(res5^2), l5))})
    ab[a, b] <- mean(result2[1, ])
  }
}
View(ab)
min(ab)

result2 <- sapply(101:200, function(i) {
  set.seed(i)
  print(i)
  ind <- sample(1:120, 100)
  fit5 <- DSIHT(X[ind, ], y[ind], group = group, L = 5, rho = 0.95, ic.coef2 = 1.6, ic.coef = 0.3)
  l <- which.min(fit5$ic)
  coef <- c(fit5$intercept[l], fit5$beta[, l])
  l5 <- length(which(fit5$beta[, l] != 0))
  res5 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res5^2)
  return(c(sum(res5^2), l5))})
matrix(apply(result2, 1, mean), ncol = 2)

result <- sapply(101:200, function(i) {
  set.seed(i)
  print(i)
  ind <- sample(1:120, 100)
  fit5 <- DSIHT(X[ind, ], y[ind], group = group, L = 5, rho = 0.9, ic.coef2 = 1.6, ic.coef = 0.3)
  l <- which.min(fit5$ic)
  coef <- c(fit5$intercept[l], fit5$beta[, l])
  l5 <- length(which(fit5$beta[, l] != 0))
  res5 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res5^2)

  fit6 <- abess(X[ind, ], y[ind], support.size = 1:20, ic.scale = 1)
  coef <- coef(fit6, support.size = fit6$best.size)
  l6 <- fit6$best.size
  res6 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res6^2)

  return(c(sum(res5^2), sum(res6^2),l5, l6))})
matrix(apply(result, 1, mean), ncol = 2)
####################################################33
library(hdi)
data(riboflavin)
y <- riboflavin$y
x <- riboflavin$x
## Data pre-processing
t <- 500
# z <- apply(x, 2, function(x) bcor(x, y))
z <- abs(apply(x, 2, function(x) cor(x, y)))
id <- which(z >= sort(z, decreasing = T)[t])
temp <- x[, id]
X <- ns(temp[, 1], df = 5)
for (i in 2:ncol(temp)) {
  X <- cbind(X, ns(temp[, i], df = 5))
}
group = rep(1:ncol(temp), each = 5)
sgl <- function(x, y, group){
  fit <- cv.sparsegl(x, y, group, nfolds = 5)
  coeff <- coef(fit, s = "lambda.min")
  beta <- coeff[-1]
  inter <- coeff[1]
  return(list(beta = beta, intercept = inter))
}

grp <- function(x, y, group, penalty){
  fit <- grpreg(x, y, group = group, penalty = penalty, dfmax = 56)
  lam <- grpreg::select(fit, crit = "EBIC")$lambda
  return(list(model = fit, lam = lam))
}

gB <- function(x, y, group){
  fit <- gBridge(x, y, group = group, dfmax = 56)
  lam <- grpreg::select(fit, crit = "EBIC")$lambda
  return(list(model = fit, lam = lam))
}

library(abess)
result <- sapply(1:50, function(i) {
  set.seed(i)
  print(i)
  ind <- sample(1:71, 60)
  m <- 71-length(ind)
  fit1 <- sgl(X[ind, ], y[ind], group = group)
  l1 <- length(which(fit1$beta != 0))
  res1 <- y[-ind] - X[-ind, ]%*%fit1$beta-fit1$intercept
  sum(res1^2)

  fit2 <- gB(X[ind, ], y[ind], group = group)
  coef <- predict(fit2$model, type = "coef", lambda = fit2$lam)
  l2 <- predict(fit2$model, type = "nvars", lambda = fit2$lam)
  res2 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res2^2)

  fit3 <- grp(X[ind, ], y[ind], group = group, penalty = "gel")
  coef <- predict(fit3$model, type = "coef", lambda = fit3$lam)
  l3 <- predict(fit3$model, type = "nvars", lambda = fit3$lam)
  res3 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res3^2)

  fit4 <- grp(X[ind, ], y[ind], group = group, penalty = "cMCP")
  coef <- predict(fit4$model, type = "coef", lambda = fit4$lam)
  l4 <- predict(fit4$model, type = "nvars", lambda = fit4$lam)
  res4 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res4^2)

  fit5 <- DSIHT(X[ind, ], y[ind], group = group, L = 10, rho = 0.95, ic.coef2 = 1)
  l <- which.min(fit5$ic)
  coef <- c(fit5$intercept[l], fit5$beta[, l])
  l5 <- length(which(fit5$beta[, l] != 0))
  res5 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res5^2)

  # fit6 <- abess(X[ind, ], y[ind])
  # coef <- coef(fit6, support.size = fit6$best.size)
  # l6 <- length(which(coef[-1] != 0))
  # res6 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  # sum(res6^2)

  return(c(sum(res1^2),
           sum(res2^2),
           sum(res3^2),
           sum(res4^2),
           sum(res5^2),
           l1, l2, l3, l4, l5))})
matrix(apply(result, 1, mean), ncol = 2)

result <- sapply(1:100, function(i) {
  set.seed(i)
  print(i)
  ind <- sample(1:71, 56)
  fit5 <- DSIHT(X[ind, ], y[ind], group = group, L = 5, rho = 0.95, ic.coef2 = 0.6, ic.coef = 2)
  l <- which.min(fit5$ic)
  coef <- c(fit5$intercept[l], fit5$beta[, l])
  l5 <- length(which(fit5$beta[, l] != 0))
  res5 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res5^2)

  fit6 <- abess(X[ind, ], y[ind], support.size = 1:20, ic.scale = 1)
  coef <- coef(fit6, support.size = fit6$best.size)
  l6 <- fit6$best.size
  res6 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res6^2)

  return(c(sum(res5^2), sum(res6^2),l5, l6))})
matrix(apply(result, 1, mean), ncol = 2)
matrix(apply(result[, 101:200], 1, mean), ncol = 2)
table(which(result[1, ] < result[2, ])%/%100)

#############################
load("C:/Mr.Zhang/RUC/data/market_x.rda")
load("C:/Mr.Zhang/RUC/data/market_y.rda")
## Data pre-processing
t <- 1000
z <- apply(x, 2, function(x) bcor(x, y))
id <- which(z >= sort(z, decreasing = T)[t])
temp <- x[, id]
X <- ns(temp[, 1], df = 5)
for (i in 2:ncol(temp)) {
  X <- cbind(X, ns(temp[, i], df = 5))
}
group = rep(1:ncol(temp), each = 5)
sgl <- function(x, y, group){
  fit <- cv.sparsegl(x, y, group, nfolds = 5)
  coeff <- coef(fit, s = "lambda.min")
  beta <- coeff[-1]
  inter <- coeff[1]
  return(list(beta = beta, intercept = inter))
}

grp <- function(x, y, group, penalty){
  fit <- grpreg(x, y, group = group, penalty = penalty, dfmax = 200)
  lam <- grpreg::select(fit, crit = "EBIC")$lambda
  return(list(model = fit, lam = lam))
}

gB <- function(x, y, group){
  fit <- gBridge(x, y, group = group, dfmax = 200)
  lam <- grpreg::select(fit, crit = "EBIC")$lambda
  return(list(model = fit, lam = lam))
}

##############################################################################
## Real-world dataset analysis
##############################################################################
i <- 2

result <- sapply(1:100, function(i) {
  set.seed(i)
  print(i)
  ind <- sample(1:464, round(0.6*464))
  m <- 464-length(ind)
  fit1 <- sgl(X[ind, ], y[ind], group = group)
  l1 <- length(which(fit1$beta != 0))
  res1 <- y[-ind] - X[-ind, ]%*%fit1$beta-fit1$intercept
  sum(res1^2)

  fit2 <- gB(X[ind, ], y[ind], group = group)
  coef <- predict(fit2$model, type = "coef", lambda = fit2$lam)
  l2 <- predict(fit2$model, type = "nvars", lambda = fit2$lam)
  res2 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res2^2)

  fit3 <- grp(X[ind, ], y[ind], group = group, penalty = "gel")
  coef <- predict(fit3$model, type = "coef", lambda = fit3$lam)
  l3 <- predict(fit3$model, type = "nvars", lambda = fit3$lam)
  res3 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res3^2)

  fit4 <- grp(X[ind, ], y[ind], group = group, penalty = "cMCP")
  coef <- predict(fit4$model, type = "coef", lambda = fit4$lam)
  l4 <- predict(fit4$model, type = "nvars", lambda = fit4$lam)
  res4 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res4^2)

  fit5 <- DSIHT(X[ind, ], y[ind], group = group, L = 5, rho = 0.9, ic.coef2 = 1)
  l <- which.min(fit5$ic)
  coef <- c(fit5$intercept[l], fit5$beta[, l])
  l5 <- length(which(fit5$beta[, l] != 0))
  res5 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res5^2)

  return(c(sum(res1^2),
           sum(res2^2),
           sum(res3^2),
           sum(res4^2),
           sum(res5^2),
           l1, l2, l3, l4, l5))})
matrix(apply(result, 1, mean), ncol = 2)

result <- sapply(1:100, function(i) {
  set.seed(i)
  print(i)
  ind <- sample(1:71, 56)
  m <- 71-length(ind)

  fit5 <- DSIHT(X[ind, ], y[ind], group = group, L = 10, rho = 0.95, ic.coef2 = 0.1)
  l <- which.min(fit5$ic)
  coef <- c(fit5$intercept[l], fit5$beta[, l])
  l5 <- length(which(fit5$beta[, l] != 0))
  res5 <- y[-ind] - X[-ind, ]%*%coef[-1]-coef[1]
  sum(res5^2)

  return(c(sum(res5^2), l5))})
apply(result, 1, mean)

