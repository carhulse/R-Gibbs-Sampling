setwd("~/STAT 462/Final Exam")
## Functions that draw from full conditionals

update.beta_0 = function(n, ybar, sig2, gamma_0, tau_0) {
  gamma0_0 = 1.0 / (n / sig2 + 1.0 / tau_0^2) 
  tau0_0 = sig2 * (n * ybar / sig2 + gamma_0 / tau_0^2) 
  rnorm(n=1, mean=gamma0_0, sd=sqrt(tau0_0))
}

update.beta_1 = function(n, ybar, sig2, gamma_1, tau_1) {
  gamma1_1 = 1.0 / (n / sig2 + 1.0 / tau_1 ^2) 
  tau1_1 = sig2 * (n * ybar / sig2 + gamma_1 / tau_1 ^2) 
  rnorm(n=1, mean=gamma1_1, sd=sqrt(tau1_1))
}
#gamma1_1 and tau1_1 are the updated gamma and tau for beta_1


update.sig2 = function(y, n, delta, nu, x, beta_0, beta_1) {
  shape1 = (delta + n)/2 
  sumsq = sum( (y-(beta_0 + beta_1*x))^2 ) 
  scale1 = (delta*nu^2 + sumsq) / 2 
  out_gamma = rgamma(n=1, shape = shape1, rate = scale1) 
  1.0 / out_gamma
}

## Perform Gibbs Sampling

gibbs = function(y, x, n_iter, init, priors) {
  ybar = mean(y) 
  n = length(y)
  
  beta0_out = numeric(n_iter)
  beta1_out = numeric(n_iter)
  sig2_out = numeric(n_iter)
  
  
  beta0_now = init$beta0
  beta1_now = init$beta1
  
  ## Gibbs sampler
  for (i in 1:n_iter){
    sig2_now = update.sig2(y=y, n=n, delta=priors$delta, nu=priors$nu, x=x, beta_0 = beta0_now, beta_1 = beta1_now) 
    beta1_now = update.beta_1(n=n, ybar=mean(y), sig2 = sig2_out, gamma_1 = priors$gamma_1, tau_1 = priors$tau_1)
    beta0_now = update.beta_0(n=n, ybar=mean(y), sig2=sig2_out, gamma_0 = priors$gamma_0, tau_0 = priors$tau_0)
    
    
    beta1_out[i] = beta1_now
    beta0_out[i] = beta0_now
    sig2_out[i] = sig2_now
  }
  cbind(beta0=beta0_out, beta1=beta1_out, sig2=sig2_out) #3 columns, n_iter rows
}

#update.beta_0(n, ybar, sig2, gamma_0, tau_0)
#update.beta_1(n, ybar, sig2, gamma_1, tau_1)

#hyperparameters gamma_0, tau_0, gamma_1, tau_1, delta, nu


data = read.table("data_q3.txt", header = TRUE, sep = "", dec = ".")
y=data[,1]
x=data[,2]

ybar=mean(y)
n=length(y)

priors = list()

#gamma_0=0, tau_0^2=5, gamma_1=0, tau_1^2=5, delta=3, nu=1

priors$gamma_0 = 0
priors$tau_0 = sqrt(5)
priors$gamma_1 = 0 
priors$tau_1 = sqrt(5)
priors$delta = 3
priors$nu = 1

init=list()
init$beta0=1
init$beta1=1

set.seed(18)
results = gibbs(y=y, x=x, n_iter=1000, init=init, priors=priors)

results.beta0 = results[,1]
results.beta1 = results[,2]
results.sig2 = results[,3]

mean(na.omit(results.beta0)) #=0.09888593
mean(na.omit(results.beta1)) #=0.1039395
mean(na.omit(results.sig2))  #=15.61165

par(mfrow=c(1,1))
hist(results.beta0, xlab="beta0", main="Histogram of beta0 estimate")
hist(results.beta1, xlab="beta1", main="Histogram of beta1 estimate")
hist(results.sig2, xlab="sig2", main="Histogram of sig2 estimate")

library(coda)
plot(as.mcmc(na.omit(results)))
summary(as.mcmc(na.omit(results)))

#BEFORE:
hist(y, freq = FALSE, xlim = c(-1.0, 3.0))
curve(dnorm(x=x, mean = priors$gamma_0, sd = priors$tau_0), lty=2, add=TRUE)

######################################################

## compare with the true standard normal 
par(mfrow=c(2,2))
X0 = seq(-100,100,by = 0.1)
hist(X[1000:N],breaks = 100, freq = FALSE,xlab="X")
dX0 = dnorm(X0)
lines(X0,dX0)
hist(Y[1000:N],breaks = 100, freq = FALSE,xlab="Y")
lines(X0,dX0)

## correlation 
cor(cbind(X,Y))

## compare with the samples generated using mvrnorm
library(MASS)
mu = c(0,0)
Sigma = rbind(c(1,rho),c(rho,1))
U = mvrnorm(N,mu,Sigma)
hist(U[,1],breaks = 100, freq = FALSE,xlab="X")
dX0 = dnorm(X0)
lines(X0,dX0)
hist(U[,2],breaks = 100, freq = FALSE,xlab="Y")
lines(X0,dX0)
## correlation 
cor(U)

