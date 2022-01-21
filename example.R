rm(list=ls())

source("mixturesplmm.R")
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
