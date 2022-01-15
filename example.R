rm(list=ls())

library(mvtnorm)
source("mixturesplmm.R")
library(penalized)
library(emulator)
library(miscTools)
library(Rcpp)
library(tidyverse)
library(glmnet)
library(ggplot2)
library(longclust)
library(MASS)

source("helpers.R")
source("splmmControl.R")
sourceCpp("splmm_rcpp.cpp")


### function sim() generateds simulated data with parameters N1, N2, p, q, ni, rho. Simulated data have 3 true fixed effects and 2 true random effects 
sim = function(N1,N2,p,q,ni,rho){
  
  # N = number of groups
  # p = number of covariates (including intercept)
  # q = number of random effect covariates
  
  ni1 <- rep(ni,N1)    # observations per group
  ntot1 <- sum(ni1)   # total number of observations
  
  ni2 <- rep(ni,N2)    # observations per group
  ntot2 <- sum(ni2)   # total number of observations
  
  
  
  ni.tot <- c(ni1,ni2)
  time <- rep(1:ni,N1+N2)
  
  grp <- factor(rep(1:(N1+N2),each=ni.tot)) # grouping variable
  
  
  beta1 <- c(2,c(1.5,0,1.8,2,0),rep(0,p-5)) # fixed-effects coefficients, true nonzero effects = 1.5, 1.8, 2.0 at position 2, 4, 5
  beta2 <- c(-2,c(-1.5,0,-1.8,-2,0),rep(0,p-5)) # fixed-effects coefficients, true nonzero effects = -1.5, -1.8, -2.0 at position 2, 4, 5
  
  covmat <- matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      covmat[i,j] <- rho^(abs(i-j))
    }
  }
  
  ntot <- ntot1+ntot2
  x <- cbind(1,rmvnorm(ntot, rep(1,p), covmat))
  
  bi <- rep(rnorm((N1+N2),0,0.5),each=ni) 
  bi12 <- rep(rnorm(N1+N2,0,1),each=ni)
  bi13 <- rep(rnorm(N1+N2,0,0.8),each=ni)
  bi22 <- rep(rnorm(N1+N2,0,0.8),each=ni)
  bi23 <- rep(rnorm(N1+N2,0,1),each=ni)
  
  bi.rest <- matrix(0,nrow = q-3,ncol = ntot)
  bi1 <- rbind(bi,bi.rest,bi12,bi13) # random effects variance, true nonzero effects = position q, q-1
  bi2 <- rbind(bi,bi.rest,bi22,bi23) # random effects variance, true nonzero effects = position q, q-1
  
  z <- x[,1:q,drop=FALSE]
  y <- numeric(ntot)
  for (k in 1:ntot1) y[k] <- x[k,]%*%beta1 + t(z[k,])%*%bi1[,grp[k]] + rnorm(1, 0, 0.5)
  
  for (k in (ntot1+1):(ntot1+ntot2)) {
    y[k] <- x[k,]%*%beta2 + t(z[k,])%*%bi2[,grp[k]] + rnorm(1, 0, 0.5)
  }
  
  
  return(list(y=y,x=x,z=z,grp=grp,time=time))
}

N1 <- 50           # number of groups in cluster 1
N2 <- 50           # number of groups in cluster 2
p <- 50            # number of covariates (excluding intercept)
q <- 10         # number of random effect covariates
ni <- 3         # number of repeated measurements per subject
rho <- 0.0      # correlation strengh between covariates


## simulated data

i=1609
set.seed(i)
data1 <- sim(N1,N2,p,q,ni,rho)


## trajectories plot by clusters

cluster <- c(rep("one",ni*N1),rep("two",ni*N2))
plotdata <- cbind.data.frame(data1$grp, data1$y, data1$time,cluster)
names(plotdata) <- c("id","y","time","cluster")
pii <- ggplot(data = plotdata, aes(x = time, y = y, group = id, color=cluster))
pii + geom_line()

## fit mixture-splmm
set.seed(i)
fit <- mixturesplmm(x=data1$x,y=data1$y,z=data1$z,grp=data1$grp,time=data1$time, lam1=c(0.15,0.1),lam2=0.5,nCluster=2,penalty.b="scad", penalty.L="lasso")

## estimated fixed-effects coefficients
fit$coefficients

## estimated random-effects variance
diag(fit$D[[1]])
diag(fit$D[[2]])

## cluster size
sum(round(fit$membership)[1,1:N1])
sum(round(fit$membership)[2,1:N1])
sum(round(fit$membership)[1,(N1+1):(N1+N2)])
sum(round(fit$membership)[2,(N1+1):(N1+N2)])

