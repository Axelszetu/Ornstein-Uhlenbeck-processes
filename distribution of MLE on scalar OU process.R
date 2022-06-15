#This script simulates the distribution of MLE of an OU process in 1 dimension.
SimScalarOU <- function(W_0 = 0, t = 10, N = 1000, theta = 0.1, sigma = 1){
  dW <- rnorm(n = N, mean = 0, sd = sqrt(t/N)*sigma)
  W <- numeric(N+1)
  W[1] <- W_0
  for (i in (1:N)){
    W[i+1] <- W[i] -theta*W[i]*(t/N) + sigma*dW[i]
  }
  W <- t(W)
  W
}

n <- 1000
theta <- 0.1
t <- 10
N <- 1000
dt <- t/N

simulated_values <- numeric(length = n)
for (i in(1:n)){
  path <- SimScalarOU(theta = theta, t = t, N = N)
  theta_hat <- OU_MLE_analytical(path, dt = dt)
  simulated_values[i] <- theta_hat
}

hist(simulated_values)