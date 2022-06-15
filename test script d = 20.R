#Test script for 20X20 drift matrix
d <- 20
t <- 100
N <- 10000
dt <- t/N
sparsity <- 0.3

A <- make_drift_matrix(d = d, sparsity = sparsity)

C_infty <- C_infty_calculator(A, ds = 0.1)

C_h <- C_hat(A, dt = dt)

process <- SimOU(A0 = A, t = t, N = N, burn = 1000)

A_mle <- OU_MLE_analytical(X = process, dt = dt)

A_Lasso <- OU_Lasso(X = process, dt = dt, lambda = 0.005)
A_Lasso

A_Dantzig <- OU_Dantzig(X = process, dt, lambda = 0.005)
A_Dantzig
