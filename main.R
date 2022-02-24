rm(list = ls())      #remove all vars
library(pracma)
library(mvtnorm)
library(coda)

set.seed(123456)
# sample parameters
n = 100
A = matrix(runif(n^2), ncol = n)
Sigma = A%*%t(A)
mu = runif(n,-5,5)
a = runif(n,-5,5)

#sample data
Z_true = mu + A%*%rnorm(n)
Y = Z_true > 0

#get the conditional standard deviations and vectors to calculate the conditional means (regr_mat)
sigma_cond = c()
regr_mat = c()
for( i in seq(n)){
  regr_vec = Sigma[i,-i]%*%inv(Sigma[-i,-i])
  regr_mat = rbind(regr_mat, regr_vec) 
  sigma_cond = c(sigma_cond, Sigma[i,i] - regr_vec%*%Sigma[-i,i])
}


#gibbs sampling
tic()
M = 1e5
burnin = M/2
Z_cur = mu + A%*%rnorm(n) #initial guess
# Z_vals = zeros(M, n)
estimate = c(0)
for( i in seq(M+burnin)){
  print(i)
  for(j in seq(n)){
    mu_cond = mu[j] + regr_mat[j,]%*%(Z_cur[-j] - mu[-j])
    
    p_0 = pnorm(0, mu_cond, sqrt(sigma_cond[j]))
    if(Y[j]==1){
      P = p_0 + runif(1)*(1-p_0)
    }else{
      P = p_0*( 1- runif(1))
    }
    Z_cur[j] = qnorm(min(1-1e-7, max(P,1e-7)), mean = mu_cond, sd =sqrt(sigma_cond[j]))
  }
  
  if(i>burnin){
    estimate = c(estimate, t(a)%*%Z_cur)
    # Z_vals[i-burnin,] = Z_cur
  }
}
t1 = toc()

#plot convergence of the estimate
plot(cumsum(estimate)/seq(M+1), type = 'l', xlab = '# iteration', ylab = 'estimate', lwd = 2)

#print the estimate
range = qnorm(0.975)*std(estimate)/effectiveSize(estimate)
c(mean(estimate)-range, mean(estimate) +range)

#print the probability of the region induced by Y
lb = sapply(Y, function(y){if(y == 1){0}else{-Inf}})
ub = sapply(Y, function(y){if(y == 1){Inf}else{0}})
pmvnorm(lb, ub, mean = mu, sigma = Sigma)


