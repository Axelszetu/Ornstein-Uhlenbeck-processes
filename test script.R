#Test scipt
 
N <- 10000
t <- 100
dt <- t/N
d <- 2
A <- matrix(data = c(1,1/2,0,1), nrow = 2)
A

test_array <- SimOU(d = d, A0 = A, t = t, N = N)
A_mle <- OU_MLE_analytical(test_array, dt = dt)
A_mle

A_mle_numeric <- OU_MLE_numeric(X = test_array, dt = dt)
A_mle_numeric
