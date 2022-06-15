# Core functions
library(lpSolve)
library(Matrix)
library(extraDistr) #for bernoulli distribution
library(MASS)


make_drift_matrix <- function(d, sparsity){
  A <- diag(d)
  n <- floor((d^2)*sparsity) - d
  values <- rbern(n = n, prob = 0.5)
  values <- (values - 0.5)*2
  
  #browser()
  values_inserted <- 0
  while(values_inserted < n){
    coord_found <- 0
    attempts <- 0
    while(coord_found == 0){
      coords <- ceiling(runif(n = 2, min = 0, max = d))
      if (A[coords[1],coords[2]] == 0){
        A[coords[1],coords[2]] <- values[values_inserted+1]
        coord_found <- 1
      }
    }
    
    if (sum(Re(eigen(A)$values) > 0) == d){
      values_inserted <- values_inserted + 1
    }
    else{
      A[coords[1],coords[2]] <- 0
    }
    attemps <- attempts + 1
    if(attempts %% 100 == 0){
      cat(sum(abs(A)), "of", d^2, "found", "\n")
    }
  }
  A
}

C_infty_calculator <- function(A, ds = 0.01, maxiter = 10000, epsilon = 0.00001){
  #browser()
  d <- ncol(A)
  integral_matrix <- matrix(data = 0, nrow = d, ncol = d)
  for (i in (1:maxiter)){
    s <- i*0.01 - 0.005
    M1 <- expm(-s*A)
    M2 <- expm(-s*t(A))
    term_in_sum <- (M1%*%M2)*ds
    integral_matrix <- integral_matrix + term_in_sum
    if ((sum((abs(term_in_sum)) > epsilon)) == 0){
      print("Convergence = 1")
      break
    }
  }

  integral_matrix
}

C_hat_calculator <- function(X){
  d <- nrow(X)
  N <- ncol(X)
  C <- matrix(data = 0, nrow = d, ncol = d)
  for (i in (1:N)){
    M <- X[,i]%*%t(X[,i])
    C <- C + M
  }
  C <- C/N
  C
}

SimOU <- function(A0, t = 10, N = 1000, burn = 1000){
  d <- ncol(A0)
  #We need a d-dmiensional brownian motion.
  samples = rnorm(n = d*(N+burn), mean = 0, sd = sqrt(t/N))
  increment_matrix <- matrix(data = samples, nrow = d, ncol = N+burn, byrow = T)
  #Create matrix to contain OU process
  OU <- matrix(nrow = d, ncol = N+1+burn)
  
  #Initializing
  #Computation of covariance of stationary distribution
  
  #init <- C_infty(A_0)
  #OU[,1] <- init
  
  OU[,1] <- numeric(d)
  
  #Computation of reusable constant
  reversion_factor <- diag(d)-A0*(t/N)
  #Computation of OU paths
  for (i in 2:(N+1+burn)){
    OU[,i] <- reversion_factor%*%OU[,i-1] + increment_matrix[,i-1]
  }
  OU <- OU[,burn:(N+burn)]
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

lasso_score_trace <- function(A, B, C, lambda, penalization = "L1"){
  d <- nrow(B)
  A <- matrix(data = A, nrow = d, ncol = d)
  
  if(penalization == "L1"){
    one_norm <- sum(abs(A))/(d^2)
    penalty <- lambda*one_norm
  }
  if(penalization == "F"){
    F_norm <- sum(A^2)/d^2
    penalty <- lambda*F_norm
  }
  
  likelihood <- sum(diag(B%*%t(A))) + sum(diag(A%*%C%*%t(A))) + penalty
  likelihood
}

OU_Lasso <- function(X, dt, lambda, penalization = "L1"){
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
  optimal_pars_vector <- optim(par = par, fn = lasso_score_trace, B = B, C = C, lambda = lambda, penalization = penalization, method = "BFGS")$par
  optimal_pars_matrix <- matrix(data = optimal_pars_vector, nrow = d, ncol = d)
  optimal_pars_matrix
}

OU_Dantzig <- function(X, dt, lambda){
  d <- nrow(X)
  N <- ncol(X)
  t <- dt*N
  dX <- apply(X, MARGIN = 1, FUN = diff)
  dX <- t(dX)
  
  M_1 <- matrix(data = 0, nrow = d, ncol = d)
  for (i in 1:(N-1)){
    M <- dX[,i]%*%t(X[,i])
    M_1 <- M_1 + M
  }
  
  M_2 <- matrix(data = 0, nrow = d, ncol = d)
  for (i in 1:(N-1)){
    M <- X[,i]%*%t(X[,i])*dt
    M_2 <- M_2 + M
  }
  
  A_mle <- -1 * M_1%*%solve(M_2)
  
  c <- rep(1, 2*d^2)
  
  diag_block <- bdiag(replicate(d, t(M_2), simplify = FALSE))
  B_left <- rbind(diag_block, -diag_block)
  B_right <- rbind(-diag_block, diag_block)
  B <- cbind(B_left, B_right)
  B <- as.matrix(B)
  
  dir <- rep("<=", 2*d^2)
  
  M_ij_vector <- as.vector(t(M_1))
  rhs_top <- -M_ij_vector + t*lambda
  rhs_bot <- M_ij_vector + t*lambda
  rhs <- c(rhs_top, rhs_bot)
  
  lp_problem <- lp(direction = "min", objective.in = c, const.mat = B, const.dir = dir, const.rhs = rhs)
  
  solution <- lp_problem$solution
  A_pos_values <- solution[1:d^2]
  A_pos <- matrix(data = A_pos_values, nrow = d, ncol = d, byrow = T)
  A_neg_values <- solution[(d^2 + 1):(2*d^2)]
  A_neg <- matrix(data = A_neg_values, nrow = d, ncol = d, byrow = T)
  
  A <- A_pos - A_neg
  
  A
}

cross_validator <- function(X, dt, f, pars, split = 0.8){
  #Setting constant values
  no_pars <- length(pars)
  scores <- numeric(length = no_pars)
  #estimator <- f
  d <- nrow(X)
  N <- ncol(X)
  
  #Splitting the sample
  split_value <- ceiling(N*split)
  training_set <- X[,(1:split_value)]
  validation_set <- X[,((split_value+1):N)]
  
  #The following will be used in the evaluation likelihood of validation set
  dtraining_set <- apply(training_set, MARGIN = 1, FUN = diff)
  dtraining_set <- t(dtraining_set)
  
  B <- matrix(data = 0, nrow = d, ncol = d)
  for (i in (1:(split_value-1))){
    M <- dtraining_set[,i]%*%t(training_set[,i])
    B <- B + M
  }
  B <- B/(split_value-1)
  
  C <- matrix(data = 0, nrow = d, ncol = d)
  for (i in (1:split_value)){
    M <- training_set[,i]%*%t(training_set[,i])
    C <- C + M
  }
  C <- C/(2*split_value)
  C <- C*dt
  
  #For each lambda, estimate A_lambda and compute (negative log-)likelihood or validation set under A_lambda
  for (i in (1:no_pars)){
    A_lambda <- f(training_set, dt, pars[i])
    scores[i] <- likelihood_evaluator_trace(A_lambda, B, C)
  }
  
  #Extracting optimal pars value
  index_of_minimum <- which.min(scores)
  optimal_parameter <- pars[index_of_minimum]
  optimal_parameter
}
