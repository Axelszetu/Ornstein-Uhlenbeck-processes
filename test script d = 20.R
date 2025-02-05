#Test script for 20X20 drift matrix
d <- 20
t <- 100
N <- 10000
dt <- t/N
sparsity <- 0.3

A <- make_drift_matrix(d = d, sparsity = sparsity)

#C_infty <- C_infty_calculator(A, ds = 0.1)

C_h <- C_hat(A, dt = dt)

process <- SimOU(A0 = A, t = t, N = N, burn = 1000)

A_hat <- OU_MLE_analytical(X = process, dt = dt)

A_Lasso <- OU_Lasso(X = process, dt = dt, lambda = 0.005)
A_Lasso

A_Dantzig <- OU_Dantzig(X = process, dt, lambda = 0.005)
A_Dantzig

#Below we test the cross validator
#We note that it takes about 11 seconds to run OU_Lasso with current parameters.
#This gives us an idea of how long we should let the cross validator run before we start to worry.


bases <- (6:-5)
lambdas <- 10^(bases)
lambda_L <- cross_validator(process, dt = dt, f = OU_Lasso, pars = lambdas, split = 0.8)
A_Lasso_CV <- OU_Lasso(X = process, dt = dt, lambda = lambda_L$lambda)
A_Lasso_CV #Looks good to me

lambda_D <- cross_validator(process, dt = dt, f = OU_Dantzig, pars = lambdas, split = 0.8)
A_Dantzig_CV <- OU_Dantzig(X = process, dt = dt, lambda = lambda_D)
A_Dantzig

#Below are tests for OU_Lasso with F-norm
A_Lasso_F <- OU_Lasso(X = process, dt = dt, lambda = lambda_L, penalization = "F")
A_Lasso_F #This looks good too