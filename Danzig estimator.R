#Script for testing lpSolve
library(lpSolve)
library(Matrix)


Dantzig_estimator <- function(X, dt, lambda){
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
  #browser()
  
  lp_problem <- lp(direction = "min", objective.in = c, const.mat = B, const.dir = dir, const.rhs = rhs)
  
  solution <- lp_problem$solution
  A_pos_values <- solution[1:d^2]
  A_pos <- matrix(data = A_pos_values, nrow = d, ncol = d, byrow = T)
  A_neg_values <- solution[(d^2 + 1):(2*d^2)]
  A_neg <- matrix(data = A_neg_values, nrow = d, ncol = d, byrow = T)
  
  A <- A_pos - A_neg
  
  A
}

lambda = 0.0005
A_Dantzig <- Dantzig_estimator(X = test_array, dt = dt, lambda = lambda)
A_Dantzig
