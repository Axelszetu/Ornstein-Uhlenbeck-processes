# Core functions
# Initial commit to test if the project is connected to GitHub

SimOU <- function(d = 1, A0 = diag(d), t = 10, N = 1000){
  #We need a d-dmiensional brownian motion.
  samples = rnorm(n = d*N, mean = 0, sd = sqrt(t/N))
  increment_matrix <- matrix(data = samples, nrow = d, ncol = N, byrow = T)
  #Create matrix to contain OU process
  OU <- matrix(nrow = d, ncol = N+1)
  #Initializing
  OU[,1] <- numeric(d)
  #Computation of reusable constant
  reversion_factor <- diag(d)-A0*(t/N)
  #Computation of OU paths
  for (i in 2:(N+1)){
    OU[,i] <- reversion_factor%*%OU[,i-1] + increment_matrix[,i-1]
  }
  OU
}

OU_MLE_analytical <- function(X, dt){
  dX <- apply(X, MARGIN = 1, FUN = diff)
}