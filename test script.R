#Test scipt
 
N <- 10000
t <- 100
dt <- t/N
A <- matrix(data = c(1,1/2,-1/2,1), nrow = 2)
A

test_array <- SimOU(A0 = A, t = t, N = N)
A_mle <- OU_MLE_analytical(X = test_array, dt = dt)
A_mle

A_mle_numeric <- OU_MLE_numeric(X = test_array, dt = dt)
A_mle_numeric

lambda = 0.5
A_lasso <- OU_Lasso(X = test_array, dt = dt, lambda = lambda)
A_lasso

