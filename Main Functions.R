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
  d <- nrow(X)
  N <- ncol(X)
  
  #Create matrix of increments
  dX <- apply(X, MARGIN = 1, FUN = diff)
  #Transpose it, because apply returns the result of each application as a column.
  dX <- t(dX)
  
  #Create matrix to iteratively contain terms of integral sum
  integral_matrix_1 <- matrix(data = 0, nrow = d, ncol = d)
  
  #Compute each step of integral sum
  for (i in 1:(N-1)){
    M <- dX[,i]%*%t(X[,i])
    integral_matrix_1 <- integral_matrix_1 + M
  }
  
  integral_matrix_2 <- matrix(data = 0, nrow = d, ncol = d)
  for (i in 1:(N-1)){
    M <- X[,i]%*%t(X[,i])*dt
    integral_matrix_2 <- integral_matrix_2 + M
  }
  integral_matrix_2_inverted <- solve(integral_matrix_2)
  
  A_hat <- -1 * integral_matrix_1%*%integral_matrix_2_inverted
  A_hat
}

likelihood_evaluator <- function(A, X, dt){
  N <- ncol(X)
  AX <- A%*%X
  dX <- apply(X, MARGIN = 1, FUN = diff)
  dX <- t(dX)
  
  terms_in_first_sum <- numeric(length = N-1)
  for (i in (1:(N-1))){
    terms_in_first_sum[i] <- crossprod(AX[,i], dt[,i])
  }
  first_sum <- sum(terms_in_first_sum)/(N-1)
  
  terms_in_second_sum <- numeric(length = N)
  for (i in (1:N)){
    terms_in_second_sum[i] <- 
  }
}

OU_MLE_numeric <- function(){
  
}