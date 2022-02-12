# cluster-splmm

This repository contains the implemented algorithm of a model-based clustering method for longitudinal data, cluster-splmm[1]. The names and descriptions are shown below:

code/clustersplmm.R: The R script contains the main function "clustersplmm" for model fitting. "clustersplmm" was developed based on the R package "splmm"[2]. For function usage please see example.R and the "splmm" mannual: https://cran.r-project.org/web/packages/splmm/splmm.pdf.

code/splmm_rcpp.cpp: The cpp file contains several functions written using Rcpp that used in the main function clustersplmm.

code/helpers.R: The R scirpt contains functions that used in the main function clustersplmm.

code/splmmControl.R: The Rscirpt defines various kinds of parameters used in the main function clustersplmm.

code/data.R: The Rscript contains a function that can generate a 2-cluster data, each cluster has 3 true fixed effects and 2 true random effects. See script for more details. Executing this file requires the following R packages:

mvtnorm

ggplot2



data1.RData: A 2-cluster data generated from data.R. See data.R for more details.

example.R: The R script demonstrates the usage of the implemented algorithm using the simulated data data1.RData. See the script for details. Executing the example requires the following R packages:

penalized

emulator

miscTools

Rcpp

tidyverse

glmnet

longclust

MASS

In addition, because the function contains C++ chunks, Rtools is required for the function to run in Windows/Mac environment. To run the example, please have the required packages installed, and download the "data1.RData", "example.R" and all files in the "code" folder under the same folder. 

[1] Yang, Luoying. Model-Based Clustering of Longitudinal Data in High Dimensions. Diss. University of Rochester, 2021.

[2] Luoying Yang and Tong Tong Wu (2021). splmm: Simultaneous Penalized Linear Mixed Effects Models. R package version 1.1.3.
  https://CRAN.R-project.org/package=splmm
