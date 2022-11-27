## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
x <- rnorm(100, 0, 4)
y <- x^3 + 1
plot(x,y)

## -----------------------------------------------------------------------------
dataset <- iris
head(dataset)
knitr::kable(head(dataset))

## -----------------------------------------------------------------------------
rpareto = function(num, a, b){
  stopifnot(a > 0 && b > 0)
  u = runif(num)
  return(b*(1-u)^(-1/a))
}

dpareto = function(x, a, b){
  stopifnot(a > 0 && b > 0)
  sapply(x, function(u) if (u<b) {0} else {a*b^a*u^(-a-1)})
}

N <- 10000
a <- 2
b <- 2

sample = rpareto(N, a, b)

hist(sample, probability = TRUE)
ran = range(sample)
x = seq(ran[1], ran[2], 0.01)
y = dpareto(x, a, b)
lines(x=x, y=y)

## -----------------------------------------------------------------------------
rbeta = function(num, df, dgen, rgen, c){
  ct = 1
  n = 1
  res = numeric(num)
  while(n <= num) {
    y = rgen(1)
    u = runif(1)
    if (u < df(y)/dgen(y)/c) {
      res[n] = y
      n = n + 1
    }
    ct = ct + 1
  }
  print(paste0("Acceptance rate: ",(n-1)/(ct-1)))
  return(res)
}

rgen = function(size){
    runif(size, min = u.min, max = u.max)
}

dgen = function(x){
    dunif(x, min = u.min, max = u.max)
}

df = function(x){
    dbeta(x = x, shape1 = a, shape2 = b)
}


a = 3
b = 2
N = 10000
 
xs = seq(0, 1, 0.01)
u.max = max(xs)
u.min = min(xs)
y.max = max(df(xs))
c = y.max / (u.max - u.min)

sample = rbeta(num = N, df = df, rgen = rgen, dgen = dgen, c)  

ys = df(xs)
  
nrBins = 50
# Compute the maximum value in the histogram for a better visualization.
bins = seq(min(xs), max(xs), by = dist(range(xs))/nrBins)
hist.vals = table(cut(x = sample, bins))
# adapt to a total area of 1 (probability histogram).
hist.vals = hist.vals/sum(hist.vals) * nrBins
# y limits as maximum of distribution and sample.
ylim = c(min(c(ys, hist.vals)), max(c(ys, hist.vals)))
  
hist(sample, probability = TRUE, ylim = ylim, breaks = nrBins)
lines(x = xs, y = ys)
  

## -----------------------------------------------------------------------------
N = 1000
r = 4
beta = 2

rmix = function(size, r, beta){
  rexp(n = N, rgamma(n = N, r, beta))
}

sample = rmix(size = size, r = r, beta = beta)

hist (sample, probability = TRUE, breaks = 100, ylim = c(0, 2))

dpareto = function(x, a, b){
  stopifnot(a > 0 && b > 0)
  sapply(x, function(u) {a*b^a*(b+u)^(-a-1)})
}

xs = seq(min(sample), max(sample), 0.01)
lines(x = xs, y = dpareto(xs, a = r, b = beta))

## -----------------------------------------------------------------------------
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}
}
n <- c(10^4, 2*10^4, 4*10^4, 6*10^4, 8*10^4)

time <- matrix(rep(0,500), nrow = 100)

for (j in 1:length(n)) {
  num <- n[j]
  for (i in 1:100) {
    data <- sample(1:num)
    time_start<-Sys.time()
    quick_sort(data)
    exc_time<-difftime(Sys.time(),time_start,units = 'mins')
    time[i,j] <- exc_time
  }
  # print(paste0(num,' running time：',round(exc_time,2),'mins'))
}

df <- cbind((colMeans(time)), matrix(rep(sapply(n, function(x){x*log(x)}), 1),nrow=5)) # to create a dataframe used for regression
colnames(df) <- c('an','tn')
df <- data.frame(df)

time_lm <- lm(an~tn, df)
plot(df$tn,df$an)
lines(df$tn,fitted(time_lm))



## -----------------------------------------------------------------------------
N <- 100000
X <- numeric(N)
Y <- numeric(N)
for (i in 1:N) {
  U <- runif(2)
  Y[i] <- (exp(U[1])+exp(1-U[1]))/2
  X[i] <- (exp(U[1])+exp(U[2]))/2
}
varX <- var(X)
varY <- var(Y)
print(paste0('The theoretical value:',round(exp(1)-1,4)))
print(paste0('The estimation calculated by antithetic variate approch:',round(mean(X),4)))
print(paste0('The estimation calculated by simple MC:',round(mean(Y),4)))
print(paste0('The empirical estimate of the percent reduction in variance:',round((varX-varY)/varX,4)*100,'%.'))

## -----------------------------------------------------------------------------
g <- function (x){
  x ^ 2 / sqrt(2*pi) * exp(-x^2/2)
}

f1 <- function(x){
  exp(-(x-1))
}

f2 <- function(x) {
  3 * x^(-4)
}


xs = seq(1,10,0.1)
y_g <- g(xs)
y_f1 <- f1(xs)
y_f2 <- f2(xs)
lim = max(c(y_g, y_f1, y_f2))

plot(xs, y_g, type = "l", ylim = c(0, lim))
lines(xs, y_f1, col="red", ylim = c(0, lim))
lines(xs, y_f2, col="blue", ylim = c(0, lim))
legend('topright', legend=c('g', 'f1', 'f2'), fill=c('black','red','blue'))
title(ylab = 'y')


## -----------------------------------------------------------------------------
m <- 10000
g <- function(x){
  exp(-x - log(1 + x^2)) * (x > 0) * (x < 1)
}
u <- runif(m)  # f3, inverse transform method
x <- -log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat510 <- mean(fg)
se510 <- sd(fg)

# Stratified importance sampling
k <- 5  # number of strata
N <- m / k  # replicates per stratum
theta.hat513 <- se513 <- numeric(k)

for (j in 1:k) {
  u <- runif(N, (j - 1) / k, j / k)
  x <- -log(1 - u * (1 - exp(-1)))
  fg <- g(x) / (k * exp(-x) / (1 - exp(-1)))
  theta.hat513[j] <- mean(fg)
  se513[j] <- sd(fg)
}
theta.hat513 <- sum(theta.hat513)
se513 <- sum(se513)

print(paste0('The estimation obtained in Example 5.10:', round(theta.hat510,6), '; The estimation obtained in Example 5.13: ', round(theta.hat513,6), '.'))
print(paste0('The estimation  se obtained in Example 5.10:', round(se510,6), '; The estimation se obtained in Example 5.13: ', round(se513,6), '.'))

## -----------------------------------------------------------------------------
exercise_6_4 <- function(seed,n){
set.seed(seed)
alpha <- 0.05
m <- 10000
cv.t <- sapply(1:m,FUN= function(o){
  y <- rnorm(n)
  c <- qt(0.975,n-1) 
  m <- mean(y) 
  se <- sqrt(var(y)) 
  as.numeric((m-c*se/sqrt(n)<0)&(m+c*se/sqrt(n)>0)) 
})
level <- mean(cv.t) 

return(level)
}

for (n in c(10,50,100,1000,10000)) {
  print(paste0('When n is ', n, ', the empirical estimate of the confidence level is: ', exercise_6_4(1,n)))
}


## -----------------------------------------------------------------------------
exercise_6_8 <- function(){
count5test <- function(x,y){
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx,outy)) > 5))
}
n <- c(20,200,1000)#分别对应小样本、中样本和大样本
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 1.5
m <- 10000
power1 <- power2 <- numeric(length(n))
set.seed(1234)
for(i in 1:length(n)){
  power1[i] <- mean(replicate(m,expr = {
    x <- rnorm(n[i],mu1,sigma1)
    y <- rnorm(n[i],mu2,sigma2)
    x <- x - mean(x)
    y <- y - mean(y)
    count5test(x,y)
  }))
  pvalues <- replicate(m,expr={
    x <- rnorm(n[i],mu1,sigma1)
    y <- rnorm(n[i],mu2,sigma2)
    Ftest <- var.test(x, y, ratio = 1,
                      alternative = c("two.sided"),
                      conf.level = 0.945, ...)
    Ftest$p.value})
  power2[i] <- mean(pvalues<=0.055)
}
df <- data.frame(power1,power2)
colnames(df) <- c('Count Five test','F test')
return(df)
}
exercise_6_8()


## -----------------------------------------------------------------------------
table <- matrix(c(6510, 3490, 10000, 6760, 3240, 10000, 13270, 6730, 20000), 3, 3,
         dimnames = list(c("Rejected", "Accepted", "total"), c("method A", "method B", "total")))
table

## -----------------------------------------------------------------------------
example_7_4 <- function(seed=123){
set.seed(seed)
library(boot)
hours <- aircondit$hours
B = 1000
# MLE 
mle_lambda <- function (x) {
  return(mean(x))
}

lambda_hat_mle <- mle_lambda(hours)

lambda_hat_b <- numeric(B)

for (b in 1:B) {
  hours_b <- sample(hours, length(hours), replace = TRUE)
  lambda_hat_b[b] <- mle_lambda(hours_b)
}

lambda_hat_b_mean <- mean(lambda_hat_b)
bias <- lambda_hat_b_mean - lambda_hat_mle
standard_error <- sd(lambda_hat_b)
print(paste0('The MLE of the rate: ', round(lambda_hat_mle,4)))
print(paste0('The estimate of the rate using bootstrap: ', round(lambda_hat_b_mean,4)))
print(paste0('The bias of estimate using bootstrap: ', bias))
print(paste0('The standard error of the estimate using bootstrap: ', round(standard_error,4)))
}

example_7_4()



## -----------------------------------------------------------------------------
example_7_5 <- function(seed=123){
set.seed(123)
library(boot) #for boot and boot.ci
hours <- aircondit$hours
theta_boot <- function(hours, index) {
  #function to compute the statistic
  mean(hours[index])
}
boot_obj <- boot(hours, statistic = theta_boot, R = 10000)
print(boot_obj)
print(boot.ci(boot_obj,
              type = c("basic", "norm", "perc",'bca')))
}
example_7_5()

## -----------------------------------------------------------------------------
mu = 4
sigma = 3
n = 10000

xs = rnorm(n, mean = mu, sd = sqrt(sigma))

# sample mean.
mu.hat = mean(xs)

# compute bootstrap sample of the mean.
B = 200
mu.hats.b = numeric(B)
ts = numeric(B)

for (b in 1:B) {
  i = sample(1:n, n, replace = TRUE)
  xs.b = xs[i]
  mu.hats.b[b] = mean(xs.b)
  
  for (b2 in 1:B) {
    i2 = sample(1:n, n, replace = TRUE)
    ts[b] = (mu.hats.b[b] - mu.hat) / sd(xs.b[i])
  }
}

se.hat = sd(mu.hats.b)

# visualize.
par(mfrow = c(2,1))
# hist(mu.hats.b, breaks = 100)
# hist(ts, breaks = 100)

# compute CIs.
alpha = 0.05
probs = c(alpha/2, 1-alpha/2)

names = sapply(probs, function(p) paste(p*100, '%', sep = ''))
setCINames = function (object) {
  return (setNames(object = object, names))
}

# standard normal.
qs.norm = qnorm(probs)
ci.sn = setCINames(mu.hat - rev(qs.norm)*se.hat)

# basic bootstrap.
qs.mu.hats.b = quantile(x = mu.hats.b, probs = probs)
ci.basic = setCINames(2*mu.hat - rev(qs.mu.hats.b))

# percentile.
ci.percentile = setCINames(quantile(mu.hats.b, probs = probs))

# set up data for the MC study.
mc.study = data.frame(rbind(ci.sn, ci.basic, ci.percentile))
colnames(mc.study) = names
mc.study['miss.left'] = rep.int(0, times = nrow(mc.study))
mc.study['miss.right'] = rep.int(0, times = nrow(mc.study))

# compute coverage rates for sample mean when sampling from the normal population xs.
size = n
rep = 10000
miss.l = 0
miss.r = 0

for(r in 1:rep) {
  i = sample(1:n, size, replace = TRUE)
  mu.sample = mean (xs[i])
  for(y in 1:nrow(mc.study)) {
    lower = mc.study[y,names[1]]
    upper = mc.study[y,names[2]]
    if (mu.sample < lower) {
      mc.study[y,'miss.left'] = mc.study[y,'miss.left'] + 1
    } else if (mu.sample > upper) {
      mc.study[y,'miss.right'] = mc.study[y,'miss.right'] + 1
    }
  }
}

mc.study$miss.left = mc.study$miss.left/rep
mc.study$miss.right = mc.study$miss.right/rep

mc.study

## -----------------------------------------------------------------------------
library(bootstrap)
Excercise_7_8 <- function(seed=1012){
sc <- scor
set.seed(seed)

theta <- function(x){ 
  sigma <- cov(x)
  pca.sigma <- prcomp(sigma) 
  theta <- pca.sigma$sdev[1] / sum(pca.sigma$sdev) 
  theta
}

n <- NROW(sc)
theta.j<- numeric(n)
for (i in 1:n){    
  theta.j[i] <- theta(sc[-i,])
}
theta.hat <- theta(sc)
bias <- (n-1) * (mean(theta.j) - theta.hat) #BIAS
se <- sqrt((n-1)*var(theta.j)) #SE
round(c(bias,se),3)
return(list(Bias.jack=bias,SE.jack=se))
}

Excercise_7_8()

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)

n <- length(magnetic)
N <- choose(n,2)
e1 <- e2 <- e3 <- e4 <- e5 <- numeric(N)

index <- 1
for (i in 1:(n-1)) for(j in (i+1):n) {
  leave <- c(i, j)
  y <- magnetic[-leave]
  x <- chemical[-leave]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coefficients[1] + J1$coefficients[2] * chemical[leave]
  e1[index] <- sum((magnetic[leave] - yhat1)^2)
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coefficients[1] + J2$coefficients[2] * chemical[leave] + 
    J2$coefficients[3] * chemical[leave] ^ 2
  e2[index] <- sum((magnetic[leave] - yhat2)^2)
  
  J3 <- lm(log(y) ~ x)
  yhat3 <- exp(J3$coefficients[1] + J3$coefficients[2] * chemical[leave])
  e3[index] <- sum((magnetic[leave] - yhat3)^2)
  
  J4 <- lm(log(y) ~ log(x))
  yhat4 <- exp(J4$coefficients[1] + J4$coefficients[2] * log(chemical[leave]))
  e4[index] <- sum((magnetic[leave] - yhat4)^2)
  
  x2 <- x^2
  x3 <- x^3
  J5 <- lm(y ~ x + x2 + x3)
  yhat5 <- J5$coefficients[1] + J5$coefficients[2] * chemical[leave] + 
    J5$coefficients[3] * chemical[leave] ^ 2 + 
    J5$coefficients[4] * chemical[leave] ^ 3
  e5[index] <- sum((magnetic[leave] - yhat5)^2)
  
  index <- index + 1
}
print(c(mean(e1), mean(e2), mean(e3), mean(e4), mean(e5)))

## -----------------------------------------------------------------------------
Excercise_8_2 <- function(){
  Spearman_rank_test <- function(x, y, B = 1e4){
    t0 = cor(x,y,method = "spearman")
    perm = numeric(B)
    z = c(x,y)
    for(i in 1:B){
      samp = sample(z)
      perm[i] = cor(samp[1:length(x)], samp[(length(x)+1):length(z)], method = "spearman")
    }
    p_value = mean(abs(perm)>=abs(t0))
    return(list(statistic = t0, 'p.value' = p_value))
  }
  print('When x and y are independent:')
  n = 1000
  x = rnorm(n, 0, 1)
  y = rnorm(n, 0, 1)
  print(Spearman_rank_test(x, y))
  print(cor.test(x, y, method = "spearman"))
  print('When x and y are dependent:')
  n = 1000
  x = rnorm(n, 5, 4)
  y = x + rnorm(n, 0, 1)
  print(Spearman_rank_test(x, y))
  print(cor.test(x, y, method = "spearman"))
}
Excercise_8_2()

## -----------------------------------------------------------------------------
sds = c(0.5, 1, 2, 4, 25, 100)

exercise_9.4 <- function (sd) {
 
  N <- 10000
  
  # standard laplace distribution.
  df <- function (x) {
    1/2 * exp(-abs(x))
  }
  
  rg = function (mean, sd) {
    rnorm(n = 1, mean = mean, sd = sd)
  }
  
  rw = function (N, df, rg) {
    
    x = numeric(N)
    x[1] = rg(0, sd)
    k = 0
    us = runif(N)
    
    for (i in 2:N) {
      xt = x[i-1]
      y = rg(xt, sd)
      res = df(y) / df(xt)
      if (us[i] <= res) {
        x[i] = y
        k = k + 1
      } else {
        x[i] = xt
      }
    }
    print(paste0('The acceptance rates is: ', k/N, ' (sd=', sd, ')'))
    return(x)
  }
  x <- rw(N, df, rg)
  return(x)
}
X <- matrix(0, nrow=length(sds), ncol=10000)
for (i in 1:length(sds))
    X[i, ] <- exercise_9.4(sds[i])

# par(mfrow = c(3,2))
# plot(1:ncol(X), X[1,], type="l", ylab = paste("sd = ", sds[1], sep = ''))
# for (i in 2:length(sds)) {
#   plot(1:ncol(X), X[i,], type="l", ylab = paste("sd = ", sds[i], sep = ''), col=i)
# }

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
}

sds = c(0.5, 1, 2, 4, 25, 100)

exercise_9.4_1 <- function (sds) {
 
  # standard laplace distribution.
  df <- function (x) {
    1/2 * exp(-abs(x))
  }
  
  rg = function (mean, sd) {
    rnorm(n = 1, mean = mean, sd = sd)
  }
  
  k <- 1
  acc_rate <- 0
  r_hat <- c()
  x <- matrix(0,nrow = length(sds), ncol = 1)
  for (i in 1:length(sds)) {
    x[i,1] <- rg(0, sds[i])
  }
  psi <- t(apply(x, 1, cumsum))
  for (i in 1:nrow(psi)){
    psi[i,] <- psi[i,] / (1:ncol(psi))
  }
  
  while (k == 1 | Gelman.Rubin(psi) >= 1.2) {
    k = k + 1
    temp <- matrix(0,nrow = length(sds), ncol = 1)
    for (i in 1:length(sds)) {
      xt <- x[i,k-1]
      y <- rg(xt, sds[i])
      res = df(y) / df(xt)
      if (runif(1) <= res) {
        temp[i,1] = y
      } else {
        temp[i,1] = xt
      }
    }
    x <- cbind(x,temp)
    psi <- t(apply(x, 1, cumsum))
    for (i in 1:nrow(psi)){
      psi[i,] <- psi[i,] / (1:ncol(psi))
    }
    r_hat[k-1] <- Gelman.Rubin(psi)
  }
  plot(1:length(r_hat), r_hat[1:length(r_hat)], type="l", ylab = 'R_hat')
  return(x)
}
x <- exercise_9.4_1(sds)





## -----------------------------------------------------------------------------
m <- 10000
burn <- 1000

x <- matrix(0, m, 2)

rho <- 0.9
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

mean12 <- function (x2) mu1 + rho*sigma1/sigma2*(x2 - mu2)
mean21 <- function (x1) mu2 + rho*sigma2/sigma1*(x1 - mu1)

x[1,] <- c(0, 0)

for (i in 2:m) {
  xt <- x[i-1,]
  xt[1] <- rnorm(1, mean12(xt[2]), s1)
  xt[2] <- rnorm(1, mean21(xt[1]), s2)
  x[i,] <- xt
}

b <- burn + 1
x_final <- x[b:m,]

X <- x_final[,1]
Y <- x_final[,2]
lin.reg <- lm(Y ~ X)
lin.reg
res <- lin.reg$residuals
plot(x, cex = 0.5, main = "generated data", ylab = 'Y', xlab = 'X')

ks.test(res,"pnorm",mean=mean(res),sd=sqrt(var(res)))

qqnorm(res,main ="Q-Q图")
qqline(res)



plot(lin.reg$fitted.values, res, cex=0.25)
abline(h = 0)
hist(lin.reg$residuals, main = "residuals of linear model")

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
}

x <- matrix(0, 2, 1)
psi <- t(apply(x, 1, cumsum))
  for (i in 1:nrow(psi)){
    psi[i,] <- psi[i,] / (1:ncol(psi))
  }

i <- 1
r_hat <- c()
rho <- 0.9
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

mean12 <- function (x2) mu1 + rho*sigma1/sigma2*(x2 - mu2)
mean21 <- function (x1) mu2 + rho*sigma2/sigma1*(x1 - mu1)



while (i <= 1000 | Gelman.Rubin(psi) >= 1.2) {
  i <- i + 1
  xt <- x[,i-1]
  temp <- matrix(0, 2, 1)
  temp[1,1] <- rnorm(1, mean12(xt[2]), s1)
  temp[2,1] <- rnorm(1, mean21(xt[1]), s2)
  x <- cbind(x, temp)
  psi <- t(apply(x, 1, cumsum))
  for (j in 1:nrow(psi)){
    psi[j,] <- psi[j,] / (1:ncol(psi))
  }
  r_hat[i-1] <- Gelman.Rubin(psi)
}

plot(r_hat)

## -----------------------------------------------------------------------------
library(mediation)
data_gen <- function(N=10000, am=1.5, ay=3.5, alpha, beta, tau=1){
  x <- rnorm(N,4,5)
  em <- rnorm(N)
  m <- alpha*x + am + em
  ey <- rnorm(N)
  y <- beta*m + tau*x + ey + ay
  return(cbind(x,m,y))
}

per_f <- function(x,m,y,R=300, para){
  model_1 <- lm(m~x)
  model_2 <- lm(y~x+m)
  re0 <- mediate(model_1, model_2, sims = 200, treat = 'x', mediator = 'm')
  t0 <- re0$d0/sd(re0$d0.sims)
  
  t <- numeric(R)
  if(para == 1){
    for (i in 1:R) {
      mstar <- sample(m,replace = FALSE)
      xstar <- x
      ystar <- y
      model_1 <- lm(mstar~xstar)
      model_2 <- lm(ystar~mstar+xstar)
      re <- mediate(model_1, model_2, sims = 200, treat = 'xstar', mediator = 'mstar')
      t[i] <- re$d0/sd(re$d0.sims)
    }
    p_value <- mean(abs(c(t0,t))>=abs(t0))
    print(paste0('Case1 p_value: ', p_value))
  }
  if(para == 2){
    for (i in 1:R) {
      mstar <- m
      xstar <- sample(x,replace = FALSE)
      ystar <- y
      model_1 <- lm(mstar~xstar)
      model_2 <- lm(ystar~mstar+xstar)
      re <- mediate(model_1, model_2, sims = 200, treat = 'xstar', mediator = 'mstar')
      t[i] <- re$d0/sd(re$d0.sims)
    }
    p_value <- mean(abs(c(t0,t))>=abs(t0))
    print(paste0('Case2 p_value: ', p_value))
  }
  if(para == 3){
    for (i in 1:R) {
      mstar <- m
      xstar <- x
      ystar <- sample(y,replace = FALSE)
      model_1 <- lm(mstar~xstar)
      model_2 <- lm(ystar~mstar+xstar)
      re <- mediate(model_1, model_2, sims = 200, treat = 'xstar', mediator = 'mstar')
      t[i] <- re$d0/sd(re$d0.sims)
    }
    p_value <- mean(abs(c(t0,t))>=abs(t0))
    print(paste0('Case3 p_value: ', p_value))
  }
}

## -----------------------------------------------------------------------------
# data <- data_gen(alpha=0,beta=0)
# per_f(data[,1],data[,2],data[,3],para=1)
# data <- data_gen(alpha=2,beta=3)
# per_f(data[,1],data[,2],data[,3],para=1)

## -----------------------------------------------------------------------------
# data <- data_gen(alpha=0,beta=1)
# per_f(data[,1],data[,2],data[,3],para=2)
# data <- data_gen(alpha=2,beta=3)
# per_f(data[,1],data[,2],data[,3],para=1)

## -----------------------------------------------------------------------------
# data <- data_gen(alpha=1,beta=0)
# per_f(data[,1],data[,2],data[,3],para=3)
# data <- data_gen(alpha=2,beta=3)
# per_f(data[,1],data[,2],data[,3],para=3)

## -----------------------------------------------------------------------------

f <- function(N,b1,b2,b3,f0){
  
  x1 <- rpois(N,1)
  x2 <- rexp(N,1)
  x3 <- rbinom(N,1,0.5)
  
  g <- function(alpha){
    temp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+temp)
    mean(p)-f0
  }
  
  solution <- uniroot(g,c(-20,0))
  return(solution$root)
}


## -----------------------------------------------------------------------------
N <- 10^6
b1 <- 0
b2 <- 1
b3 <- -1
f0s <- c(0.1,0.01,0.001,0.0001)
alphas <- c()
for (f0 in f0s) {
  set.seed(123)
  alpha <- f(N,b1,b2,b3,f0)
  alphas <- c(alphas, alpha)
  print(alpha)
}


## -----------------------------------------------------------------------------
plot(f0s,alphas, xlab = 'f0', ylab = 'alpha')

## -----------------------------------------------------------------------------
u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3) 

g <- function(lambda){
  u <- c(11,8,27,13,16,0,23,10,24,2)
  v <- c(12,9,28,14,17,1,24,11,25,3)
  re <- 0
  n <- length(u)
  for (i in 1:n) {
    re <- re + (u[i]*exp(-lambda*u[i])-v[i]*exp(-lambda*v[i]))/(exp(-lambda*u[i])-exp(-lambda*v[i]))
  }
  return(re)
}

lambda_mle <- uniroot(g,c(0,10))$root
print(paste0('The numerical solution of MLE of lambda is: ',lambda_mle))

iter <- function(lambda){
  u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
  v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3) 
  re <- 0
  n <- length(u)
  for (i in 1:n) {
    re <- re + 1/lambda + (u[i]*exp(-lambda*u[i])-v[i]*exp(-lambda*v[i]))/(exp(-lambda*u[i])-exp(-lambda*v[i]))
  }
  return(n/re)
}
f <- function(lambda0){
  k <- 2
  n <- length(u)
  lambda_old <- iter(lambda0)
  lambda_new <- iter(lambda_old)
  while (abs(lambda_old-lambda_new)>=1e-6) {
    lambda_old <- lambda_new
    lambda_new <- iter(lambda_old)
    k <- k + 1
  }
  return(lambda_new)
}
lambda_em <- f(1)

print(paste0('The numerical solution of MLE of lambda obtained by EM algorithm is: ',lambda_em))

## -----------------------------------------------------------------------------
a <- list(list(1,2,3,4),'a')
print('result of as.vector():')
print(as.vector(a))
print('result of unlist():')
print(unlist(a))

## -----------------------------------------------------------------------------
dim(c(1,2,3,4))
dim(array(c(1,2,3,4),1))

## -----------------------------------------------------------------------------
a=matrix(c(1,2,3,'4',5,6),nrow=2)
is.matrix(a)
is.array(a)

## -----------------------------------------------------------------------------
a <- data.frame('num'=c(1,2,3),nrow='a',log=TRUE)
a
as.matrix(a)

## -----------------------------------------------------------------------------
a <- data.frame()

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}
head(cars)
data.frame(lapply(cars, function(x) if (is.numeric(x)) scale01(x) else x))
rm(list = ls())

## -----------------------------------------------------------------------------
head(cars)
vapply(cars, sd, numeric(1))
rm(list = ls())

## -----------------------------------------------------------------------------
head(iris)
vapply(iris[vapply(iris, is.numeric, logical(1))],sd, numeric(1))
rm(list = ls())

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

gibbsR <- function(m,burn,rho,mu1,mu2,sigma1,sigma2){
  x <- matrix(0, m, 2)

  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  
  mean12 <- function (x2) mu1 + rho*sigma1/sigma2*(x2 - mu2)
  mean21 <- function (x1) mu2 + rho*sigma2/sigma1*(x1 - mu1)
  
  x[1,] <- c(0, 0)
  
  for (i in 2:m) {
    xt <- x[i-1,]
    xt[1] <- rnorm(1, mean12(xt[2]), s1)
    xt[2] <- rnorm(1, mean21(xt[1]), s2)
    x[i,] <- xt
  }
  
  b <- burn + 1
  x_final <- x[b:m,]
}


data_dir <- getwd()
# sourceCpp(paste0("./gibbsC.cpp"))
set.seed(123)
gibbs_R <- gibbsR(10000,1000,0.9,0,0,1,1)
set.seed(123)
# gibbs_C <- gibbsC(10000,1000,0.9,0,0,1,1)


## -----------------------------------------------------------------------------
# qqplot(gibbs_C,gibbs_R)

## -----------------------------------------------------------------------------
# ts <- microbenchmark(gibbC=gibbsC(10000,1000,0.9,0,0,1,1),gibbR=gibbsR(10000,1000,0.9,0,0,1,1))
# summary(ts)[,c(1,3,5,6)]

