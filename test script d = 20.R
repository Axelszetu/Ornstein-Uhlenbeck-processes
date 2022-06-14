#Test script for 20X20 drift matrix
d <- 20
t <- 100
N <- 10000
sparsity <- 0.3

A <- make_drift_matrix(d = d, sparsity = sparsity)
