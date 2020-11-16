#Author: Zishu Zhan, Xiangjie Li
#Update: 2020.11.12
#Contact: 13023988263@163.com (Zishu Zhan), ele717@163.com (Xiangjie Li)
#Author: Zishu Zhan
#Update: 2020.11.12
#function includes: kern_fun, indx_comp, kern_ilse, kern_prime, est_fun

#the reproduce script for case 1 in PRIME paper

R <- 0.7 # R^2=0.2
n <- 100 # sample size equals to 100
mr <- '60%' # missing rate equals to 60%NS
NS <- 100
mis.para <- c(-1,-2,0.1)
Set <- paste('MR','=',mr,', ','n','=',n,', ','R-square','=',R)
p <- 12
B <- 100
beta0 <- c(1, -2, 1.5, 1, 1.2, 2, -1.8, -3.1, 1.3, 0.5, 1.1, -1.4)
Sigma0 <- matrix(0.5,p,p)
diag(Sigma0) <- 1
sig <- 1
distr <- 'rchange'
mu0 <- rep(0,p)
pattern_mis <- list(c(1,2,3),c(4,5,6),c(7,8,9),c(1:6),c(1,2,3,7,8,9),c(4:9),
                    c(1,4),c(1,7),c(4,7),c(3),c(6),c(9))

beta_all <- list()
pre_all <- list()
beta_ind <- list()
beta_rank <- list()
num_cc <- vector()
beta_mse <- matrix(0, NS, 5)
colnames(beta_mse) <- c('full', 'prime','ilse','cc','ml')
colnames(beta_dis) <- c('full', 'prime','ilse','cc','ml')

for(k in 1:NS){
  print(k)
  X <- mvrnorm(n, mu=mu0, Sigma = Sigma0)
  X_mis <- X
  if (distr=='normal'){
    eps <- rnorm(n, 0, sig)  
    R <- var(X%*%beta0)/(var(X%*%beta0)+sig^2)
  }
  
  if (distr=='hete'){
    eps <- rnorm(n, 0, sig)*X[,6]
    R <- var(X%*%beta0)/(var(X%*%beta0)+sig^2)
  }
  
  if (distr=='chi'){
    eps <- rchisq(n, 3)
    R <- var(X%*%beta0)/(var(X%*%beta0)+sig^2)
  }
  
  if (distr=='rchange'){
    sig <- sqrt(var(X%*%beta0)*(1-R)/R)
    eps <- rnorm(n, 0, sig)  
  }
  #rate_mis <- 1/(1+exp(-1.7*X[,6]-3)) 
  rate_mis <- 1/(1+exp(mis.para[1]*eps+mis.para[2])) 
  delta1 <- rbinom(n, 1, mis.para[3])
  delta2 <- rbinom(n, 1, rate_mis)
  index_mis <- rep(0, n)
  ll <- length(pattern_mis)
  for (i in 1:n) {
    if (delta1[i] == 1) {
      #randomly assign missing patterns
      index_mis[i] <- sample((1:round(ll/2,0)), 1)
      pattern <- pattern_mis[[index_mis[i]]]
      X_mis[i, pattern] <- NA
    }
    
    if (delta2[i] == 1) {
      #randomly assign missing patterns
      index_mis[i] <- sample(((round(ll/2,0)+1):ll), 1)
      pattern <- pattern_mis[[index_mis[i]]]
      X_mis[i, pattern] <- NA
    }
  }
  
  Y <- X%*%beta0 + eps
  res.ilse <- est_fun(Y, X_mis, func='ilse', method='Par.cond', intercept = F)
  beta.ilse <- res.ilse$beta
  res.rds <- est_fun(Y, X_mis, func='rds', method='Par.cond', intercept = F)
  beta.rds <- res.rds$beta
  X_cc <- na.omit(X_mis)
  data.cc <- data.frame(na.omit(cbind(Y,X_mis)))
  num_cc[k] <- nrow(data.cc)
  print(paste(Set,(n-num_cc[k])/n)) 
  if(nrow(data.cc) >= p) beta.cc <- lm(X1~.+0, data=data.cc)$coef
  if(nrow(data.cc) < p) beta.cc <- as.vector(ginv(t(X_cc)%*%X_cc)%*%t(X_cc)%*%data.cc[,1])
  data.mis <- data.frame(Y,X_mis)
  beta.ml <- miss.lm(Y~.+0, data=data.mis)$coef
  beta.full <- as.vector(solve(t(X)%*%X)%*%t(X)%*%Y)
  beta_mse[k, 1] <- sum((beta.full - beta0)^2)
  beta_mse[k, 2] <- sum((beta.rds - beta0)^2)
  beta_mse[k, 3] <- sum((beta.ilse - beta0)^2)
  beta_mse[k, 4] <- sum((beta.cc - beta0)^2)
  beta_mse[k, 5] <- sum((beta.ml - beta0)^2)
  
  beta_ind[[k]] <- rbind(abs((beta.full - beta0)/beta0), abs((beta.rds - beta0)/beta0), 
                         abs((beta.ilse - beta0)/beta0), abs((beta.cc - beta0)/beta0), 
                         abs((beta.ml - beta0)/beta0))
  beta_rank[[k]] <- apply(beta_ind[[k]],2,rank)
  
  beta_all[[k]] <- matrix(c(beta.full, beta.rds, beta.ilse, beta.cc, beta.ml), nrow = 5, byrow = T)
}

