# mixture-splmm

This repository contains the implemented algorithm of a model-based clustering method for longitudinal data, mixture-splmm[1]. The names and descriptions are shown below:

code/mixturesplmm.R: The R script contains the main function "mixturesplmm" for model fitting. For function usage please see example.R.

code/splmm_rcpp.cpp: The cpp file contains several functions written using Rcpp that used in the main function mixturesplmm.

code/helpers.R: The R scirpt contains functions that used in the main function mixturesplmm.

code/splmmControl: The Rscirpt defines various kinds of parameters used in the main function mixturesplmm.

example.R: The R script simulates a simple dataset to demonstrate the usage of the implemented algorithm. See the script for details. Executing the example requires the following R packages:

library(mvtnorm)

library(penalized)

library(emulator)

library(miscTools)

library(Rcpp)

library(tidyverse)

library(glmnet)

library(ggplot2)

library(longclust)

library(MASS)

To run the example, please have the required packages installed, and download all code scripts under the same folder. 

[1] Yang, Luoying. Model-Based Clustering of Longitudinal Data in High Dimensions. Diss. University of Rochester, 2021.
