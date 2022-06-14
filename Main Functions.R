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

likelihood_evaluator_trace <- function(A, B, C){
  d <- nrow(B)
  A <- matrix(data = A, nrow = d, ncol = d)
  
  likelihood <- sum(diag(B%*%t(A))) + sum(diag(A%*%C%*%t(A)))
  likelihood
}

OU_MLE_numeric <- function(X, dt){
  d <- nrow(X)
  N <- ncol(X)
  dX <- apply(X, MARGIN = 1, FUN = diff)
  dX <- t(dX)
  
  B <- matrix(data = 0, nrow = d, ncol = d)
  for (i in (1:(N-1))){
    M <- dX[,i]%*%t(X[,i])
    B <- B + M
  }
  B <- B/(N-1)
  
  C <- matrix(data = 0, nrow = d, ncol = d)
  for (i in (1:N)){
    M <- X[,i]%*%t(X[,i])
    C <- C + M
  }
  C <- C/(2*N)
  C <- C*dt
  
  par <- diag(d)
  par <- as.vector(par)
  optimal_pars_vector <- optim(par = par, fn = likelihood_evaluator_trace, B = B, C = C)$par
  optimal_pars_matrix <- matrix(data = optimal_pars_vector, nrow = d, ncol = d)
  optimal_pars_matrix
}

lasso_score_trace <- function(A, B, C, lambda){
  d <- nrow(B)
  A <- matrix(data = A, nrow = d, ncol = d)
  
  one_norm <- sum(abs(A))/(d^2)
  penalty <- lambda*one_norm
  
  likelihood <- sum(diag(B%*%t(A))) + sum(diag(A%*%C%*%t(A))) + penalty
  likelihood
}

OU_Lasso <- function(X, dt, lambda){
  d <- nrow(X)
  N <- ncol(X)
  dX <- apply(X, MARGIN = 1, FUN = diff)
  dX <- t(dX)
  
  B <- matrix(data = 0, nrow = d, ncol = d)
  for (i in (1:(N-1))){
    M <- dX[,i]%*%t(X[,i])
    B <- B + M
  }
  B <- B/(N-1)
  
  C <- matrix(data = 0, nrow = d, ncol = d)
  for (i in (1:N)){
    M <- X[,i]%*%t(X[,i])
    C <- C + M
  }
  C <- C/(2*N)
  C <- C*dt
  
  par <- diag(d)
  par <- as.vector(par)
  optimal_pars_vector <- optim(par = par, fn = lasso_score_trace, B = B, C = C, lambda = lambda)$par
  optimal_pars_matrix <- matrix(data = optimal_pars_vector, nrow = d, ncol = d)
  optimal_pars_matrix
}

