#Working example for Dantzig estimator
rm(list = ls())
library(lpSolve)
library(Matrix)
library(extraDistr) #for bernoulli distribution


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

SimOU <- function(A0, t = 10, N = 1000){
  d <- ncol(A0)
  #We need a d-dmiensional brownian motion.
  samples = rnorm(n = d*N, mean = 0, sd = sqrt(t/N))
  increment_matrix <- matrix(data = samples, nrow = d, ncol = N, byrow = T)
  #Create matrix to contain OU process
  OU <- matrix(nrow = d, ncol = N+1)
  
  #Initializing
  #Computation of covariance of stationary distribution
  
  #init <- C_infty(A_0)
  #OU[,1] <- init
  
  OU[,1] <- numeric(d)
  
  #Computation of reusable constant
  reversion_factor <- diag(d)-A0*(t/N)
  #Computation of OU paths
  for (i in 2:(N+1)){
    OU[,i] <- reversion_factor%*%OU[,i-1] + increment_matrix[,i-1]
  }
  OU
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


N <- 10000
t <- 100
dt <- t/N
lambda <- 0.005
d <- 6
sparsity <- 0.35


A2 <- make_drift_matrix(d = d, sparsity = sparsity)
test2 <- SimOU(A0 = A2, t = t, N = N)

A2_Dantzig <- OU_Dantzig(X = test2, dt = dt, lambda = lambda)
A2_Dantzig
A2