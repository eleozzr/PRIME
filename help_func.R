#Author: Zishu Zhan, Xiangjie Li
#Update: 2020.11.12
#function includes: kern_fun, indx_comp, kern_ilse, kern_rds, est_fun
#Contact: 13023988263@163.com (Zishu Zhan), ele717@163.com (Xiangjie Li)
#Author: Zishu Zhan
#Update: 2020.11.12
#function includes: kern_fun, indx_comp, kern_ilse, kern_prime, est_fun

library(MASS)
library(misaem)
library(ggplot2)
library(gridExtra)

# kern_fun function has already been defined in the package 'ILSE' and is used to return kernel result 
# e.g.
# data.kernel <- data.fun(u, type = "gaussian")
kern_fun <- function (u, type = "gaussian") 
{
  if (!is.element(type, c("epk", "biweight", "triangle", 
                          "gaussian", "triweight", "tricube", 
                          "cosine", "uniform",'eucle'))) 
    stop("type must belong to 'epk', 'biweight', 'triangle',\n         'guassian', 'triweight', 'tricube', 'cosine', 'uniform'!")
  if (type == "epk") 
    f = 0.75 * (1 - u^2) * (u <= 1 & u >= -1)
  else if (type == "biweight") 
    f = (15/16) * (1 - u^2)^2 * (u <= 1 & u >= -1)
  else if (type == "triangle") 
    f = (1 - abs(u)) * (u <= 1 & u >= -1)
  else if (type == "gaussian") 
    f = dnorm(u)
  else if (type == "triweight") 
    f = 35/32 * (1 - u^2)^3 * (u <= 1 & u >= -1)
  else if (type == "tricube") 
    f = 70/81 * (1 - abs(u)^3)^3 * (u <= 1 & u >= -1)
  else if (type == "cosine") 
    f = pi/4 * cos(pi * u/2)
  else if (type == "uniform") 
    f = 1/2 * (u <= 1 & u >= -1)
  return(f)
}

# indx_comp function has already been defined in the package 'ILSE'
# and is used to the get the all the indexes of units i' with missing data pattern Pi' including Pi and j
indx_comp <- function (Xmat) 
{
  p <- ncol(Xmat)
  IDX <- NULL
  for (j in 1:p) {
    idxj <- list(which(!is.na(Xmat[, j])))
    IDX <- c(IDX, idxj)
  }
  IDX
}

# kern_ilse function has already been defined in the package 'ILSE' and is used to get the imputation results by ILSE
kern_ilse <- function (ind, beta, Xmat, Y, IDX, bw, k.type = NULL, K, bw.type) 
{
  p <- ncol(Xmat)
  if (is.null(k.type)) 
    k.type = "gaussian"
  Zi <- Xmat[ind, ]
  Pi <- which(!is.na(Xmat[ind, ]))
  w <- sum(Xmat[ind, Pi] * beta[Pi])
  i.peer <- Reduce(intersect, IDX[Pi])
  Pi.peer <- which(is.na(Xmat[ind, ]))
  for (j in 1:p) {
    ij.peer <- intersect(i.peer, IDX[[j]])
    n_ijp <- length(ij.peer)
    if (bw.type == "var.bw") {
      Xbeta <- matrix(Xmat[ij.peer, Pi], nrow = n_ijp, 
                      ncol = length(Pi)) %*% matrix(beta[Pi], nrow = length(Pi), 
                                                    ncol = 1)
      Xbeta <- Xbeta[, 1]
      Ind.knn <- order(abs(Xbeta - w))[1:K]
      h <- (max(Xbeta[Ind.knn], na.rm = T) - min(Xbeta[Ind.knn], 
                                                 na.rm = T))/2
      h <- max(h, min(abs(Xmat[ij.peer, Pi] %*% beta[Pi] - 
                            w)) + 0.01)
    }
    if (bw.type == "fix.bw") {
      Xbeta <- matrix(Xmat[ij.peer, Pi], nrow = length(ij.peer), 
                      ncol = length(Pi)) %*% matrix(beta[Pi], nrow = length(Pi), 
                                                    ncol = 1)
      Xbeta <- Xbeta[, 1]
      h <- sd(Xbeta + rnorm(n_ijp, sd = 1e-05)) * bw
    }
    u <- (Xbeta - w)/(h)
    if (sum(u <= 1 & u >= -1) == 0 & k.type != "gaussian") {
      warning("Warning: the bandwidth is too small.. \n")
      Zi <- NULL
      return(Zi)
    }
    u.kern <- kern_fun(u, type = 'gaussian')
    Zi[j] <- sum(Xmat[ij.peer, j] * u.kern)/sum(u.kern)
  }
  return(Zi)
}

# geom_mean function is used to get the geometric mean
geo_mean <- function(data){
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

# kern_prime function is used to get the imputation results by PRIME
kern_prime <- function (ind, Xmat, Y, IDX, bw, k.type = NULL, K, bw.type) 
{
  Sigma <- diag(rep(1,p))
  t <- mvrnorm(B, rep(0,p), Sigma)
  #t <- t(replicate(B,runif(p,-1,1)))
  #tssq <- apply(t0^2, 1, sum)
  #t <- t0/matrix(rep(sqrt(tssq), p), B)
  #t <- t0
  p <- ncol(Xmat)
  if (is.null(k.type)) 
    k.type = "gaussian"
  ZI <- Xmat[ind, ]
  Pi <- which(!is.na(Xmat[ind, ]))
  #w <- sum(Xmat[ind, Pi] * beta[Pi])
  i.peer <- Reduce(intersect, IDX[Pi])
  Pi.peer <- which(is.na(Xmat[ind, ]))
  
  for(j in 1:p){
    ij.peer <- intersect(i.peer, IDX[[j]])
    Wb <- matrix(0, B, length(ij.peer))
    for (ii in 1:B) {
      w <- sum(Xmat[ind, Pi] * t[ii,Pi])
      n_ijp <- length(ij.peer)
      if (bw.type == "var.bw") {
        Xt <- matrix(Xmat[ij.peer, Pi], nrow = n_ijp, 
                     ncol = length(Pi)) %*% matrix(t[ii,Pi], nrow = length(Pi), 
                                                   ncol = 1)
        Xt <- Xt[, 1]
        Ind.knn <- order(abs(Xt - w))[1:K]
        h <- (max(Xt[Ind.knn], na.rm = T) - min(Xt[Ind.knn], 
                                                na.rm = T))/2
        h <- max(h, min(abs(Xmat[ij.peer, Pi] %*% t[ii,Pi] - 
                              w)) + 0.01)
      }
      if (bw.type == "fix.bw") {
        Xt <- matrix(Xmat[ij.peer, Pi], nrow = length(ij.peer), 
                     ncol = length(Pi)) %*% matrix(t[ii,Pi], nrow = length(Pi), 
                                                   ncol = 1)
        Xt <- Xt[, 1]
        h <- sd(Xt + rnorm(n_ijp, sd = 1e-05)) * bw
      }
      u <- (Xt - w)/h
      if (sum(u <= 1 & u >= -1) == 0 & k.type != "gaussian") {
        warning("Warning: the bandwidth is too small.. \n")
        Zi <- NULL
        return(Zi)
      }
      u.kern <- kern_fun(u, type = 'gaussian')
      Wb[ii,] <- u.kern
    } 
    
    Wij <- apply(Wb, 2, geo_mean)
    ZI[j] <- sum(Xmat[ij.peer, j] * Wij)/sum(Wij)
  }
  return(ZI)
}

# est_fun function is used to get the estimation results of ILSE and PRIME
est_fun <- function (Y, X, func, bw = NULL, intercept = F, k.type = NULL, K=NULL, bw.type = "fix.bw", 
                      method = "Par.cond", max.iter = 100, peps = 1e-05, feps = 1e-07, 
                      infor_output = F) 
{
  if (is.null(Y) || is.null(X)) 
    stop("\"X\",\"Y\" must be given\n simutaneously!")
  if (is.null(n <- nrow(X))) 
    stop("'X' must be a matrix")
  if (n == 0L) 
    stop("0 (non-NA) cases")
  if (!is.null(bw) && !is.numeric(bw)) 
    stop("\"bw\" must be NULL or a positive scalar!")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  if (intercept) {
    X <- cbind(1, X)
    colnames(X)[1] <- "intercept"
  }
  p <- ncol(X)
  i.com <- which(complete.cases(X))
  Xmat0 <- X
  Xmat0[is.na(X)] <- 0
  Z2 <- Xmat0
  IDX <- indx_comp(X)
  cc.coef <- as.numeric(lm(Y ~ . + 0, data.frame(X), subset = i.com)$coef)
  if (any(is.na(cc.coef))) 
    cc.coef[is.na(cc.coef)] <- 1
  if (is.null(bw)) 
    bw <- n^(-1/3) * sd(Xmat0 %*% matrix(cc.coef, p, 1))
  bw <- ifelse(bw > 1, bw, 1)
  beta <- cc.coef
  Bmat <- beta
  rss <- rss.old <- sqrt(sum((Y - Xmat0 %*% matrix(beta, p, 1))^2))
  k <- 0
  if (infor_output == T) {
    cat("iter:", k, "d.fn:", NA, "d.par:", 
        NA, "\n")
    cat("par:", as.vector(round(beta, 2)), "\n")
  }
  k <- k + 1
  NA_ind <- which(apply(is.na(X), 1, sum) != 0)
  if(func == 'ilse'){
    while (k < max.iter) {
      W <- Xmat0 %*% beta
      Z1 <- X
      for (i in NA_ind) {
        temp_Z1 <- kern_ilse(ind = i, beta = beta, Xmat = X, 
                             Y = Y, IDX = IDX, bw = bw, k.type = k.type, K = K, 
                             bw.type = "fix.bw")
        if (is.null(temp_Z1)) {
          return(NULL)
        }
        else {
          Z1[i, ] <- temp_Z1
        }
      }
      Z2[is.na(X)] <- Z1[is.na(X)]
      if (method == "Full.cond") {
        xtx <- t(Z1) %*% Z2
        xty <- t(Z1) %*% Y
        beta.new <- solve(xtx) %*% xty
      }
      else if (method == "Par.cond") {
        beta.new <- lm(Y ~ . + 0, data.frame(Z2))$coef
      }
      rss.new <- sqrt(sum((Y - Xmat0 %*% beta.new)^2))
      d.fn <- abs(rss.new - rss.old)/rss.old
      d.par <- max(abs(beta.new - beta))/max(abs(beta))
      if (infor_output == T) {
        cat("iter:", k, "d.fn:", d.fn, "d.par:", 
            d.par, "\n")
        cat("par:", as.vector(round(beta.new, 2)), 
            "\n")
      }
      beta <- beta.new
      Bmat <- rbind(Bmat, matrix(beta.new, nrow = 1))
      rss.old <- rss.new
      rss <- c(rss, rss.new)
      if ((method != "SGD" & d.par < peps) | (method != 
                                              "SGD" & d.fn < feps)) {
        break
      }
      else if (method == "SGD" & max(abs(2 * t(Z2) %*% 
                                         (Y - Z2 %*% beta.new))) < 0.5) {
        break
      }
      k <- k + 1
      if (k == max.iter) {
        warning("algorithm may not converge!")
      }
    }
    residuals <- lm(Y ~ . + 0, data.frame(Z2))$residuals
    #row.names(Bmat) <- paste0("iter", 0:k)
  }
  
  if(func == 'rds'){
    #W <- Xmat0 %*% beta
    Z1 <- X
    for (i in NA_ind) {
      temp_Z1 <- kern_prime(ind = i, Xmat = X, 
                          Y = Y, IDX = IDX, bw = bw, k.type = k.type, K = K, 
                          bw.type = "fix.bw")
      if (is.null(temp_Z1)) {
        return(NULL)
      }
      else {
        Z1[i, ] <- temp_Z1
      }
    }
    Z2[is.na(X)] <- Z1[is.na(X)]
    if (method == "Full.cond") {
      xtx <- t(Z1) %*% Z2
      xty <- t(Z1) %*% Y
      beta.new <- solve(xtx) %*% xty
    }
    if (method=="Full.repl") {
      xtx <- t(Z1) %*% Z1
      xty <- t(Z1) %*% Y
      beta.new <- solve(xtx) %*% xty
    }
    else if (method == "Par.cond") {
      beta.new <- lm(Y ~ . + 0, data.frame(Z2))$coef
    }
    rss.new <- sqrt(sum((Y - Xmat0 %*% beta.new)^2))
    d.fn <- NULL
    d.par <- NULL
    
    if (infor_output == T) {
      cat("iter:", k, "d.fn:", d.fn, "d.par:", 
          d.par, "\n")
      cat("par:", as.vector(round(beta.new, 2)), 
          "\n")
    }
    Bmat <- NULL
    residuals <- NULL
  }
  
  res <- list(beta = as.vector(beta.new), Bmat = Bmat, d.fn = d.fn, 
              d.par = d.par, iterations = k, residuals = residuals, 
              inargs = list(bw = bw, intercept = intercept, k.type = k.type, bw.type = bw.type, 
                            K = K, method = method, max.iter = max.iter, peps = peps, 
                            feps = feps, infor_output = infor_output))
  return(res)
}
