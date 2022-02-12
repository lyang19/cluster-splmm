rm(list=ls())

source("clustersplmm.R")
library(penalized)
library(emulator)
library(miscTools)
library(Rcpp)
library(tidyverse)
library(glmnet)
library(longclust)
library(MASS)

source("helpers.R")
source("splmmControl.R")
sourceCpp("splmm_rcpp.cpp")
load("data1.RData")

## fit cluster-splmm
i=1609
set.seed(i)
## nCluster is the number of cluster input to the function
## penalty.b is the penalty for the fixed effects
## penalty.L is the penalty for the random effects
## lam1 and lam2 are tuning parameters for fixed effects and random effects respectively, and they can be input as either vector or scalar 
## in this example, lam1= c(0.15,0.1) is a vector where 0.15 is for the 1st cluster and 0.1 is for the 2nd cluster
## lam2 = 0.5 is a scalar, and it will be populated to a nCluster-length vector of the input value. In this example, the lam2 actually used in the function is c(0.5,0.5) 

fit <- clustersplmm(x=data1$x,y=data1$y,z=data1$z,grp=data1$grp,time=data1$time, lam1=c(0.15,0.1),lam2=0.5,nCluster=2,penalty.b="scad", penalty.L="lasso")

## estimated fixed-effects coefficients
fit$coefficients

## estimated random-effects variance
diag(fit$D[[1]])
diag(fit$D[[2]])

## cluster size
N1 <- 50 
N2 <- 50
sum(round(fit$membership)[1,1:N1])
sum(round(fit$membership)[2,1:N1])
sum(round(fit$membership)[1,(N1+1):(N1+N2)])
sum(round(fit$membership)[2,(N1+1):(N1+N2)])
