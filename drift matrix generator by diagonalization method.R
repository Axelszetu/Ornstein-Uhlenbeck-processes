# Making sparse drift matrix
library(extraDistr) #For discrete uniform
library(gtools) #For permute function

make_drift_matrix <- function(d=4, s=3){
  valid_matrix_found <- 0
  while(valid_matrix_found == 0){
    #browser()
    #Making diagonal matrix
    D <- matrix(data = 0, nrow = d, ncol = d)
    D_diagonal_terms <- rchisq(n = d, df = 3)
    diag(D) <- D_diagonal_terms
    #Making invertible matrix P
    P <- matrix(data = 0, nrow = d, ncol = d)
    P_diagonal_terms <- rnorm(n = d)
    diag(P) <- P_diagonal_terms
    #Making P s-sparse
    for (i in (1:d)){ #In each row of the matrix, do:
      row_nonzeros <- rdunif(1, 0, s-1) #Decide on the number of nonzero elements off the diagonal.
      nonzero_values <- rnorm(n = row_nonzeros)
      #We now need to make a permutation to decide the position of the nonzero entries, that preserves the diagonal.
      # Idea: Make a permutation of length d-1, slice at position i and insert the diagonal value.
      zeros <- numeric(length = d-row_nonzeros-1)
      nondiagonal_values <- c(nonzero_values, zeros)
      nondiagonal_values <- permute(nondiagonal_values)
      if (i == 1){
        P[1,] <- c(P_diagonal_terms[1], nondiagonal_values)
      } else if (i == d){
        P[d,] <- c(nondiagonal_values, P_diagonal_terms[d])
      } else{
        row_values <- c(nondiagonal_values[1:(i-1)], P_diagonal_terms[i], nondiagonal_values[i:(d-1)])
        P[i,] <- row_values
      }
      # The above should work for values 2 to d-1. I don't think it works for i=1. Unsure about i=d.
      # We can maybe fix this with if-statements for i = 1 and i = d.
      # As you can see, I've tried doing that. Let's see if it works.
      # It works!
      
      # Next problem:
      # When computing the product PDP^-1, we get some entries that seem to be numerically zero.
      # Maybe we can fix this using the Sparse class in the Matrix package.
      #Another idea: if the value is close enough to 0, we simply manually make it zero and see if we retain positive eigenvalues.
    }
    A <- P%*%D%*%solve(P)
    A[abs(A) < 0.5] <- 0
    no_negatives <- sum(A < 0)
    eigen_values <- eigen(A)$values
    if((sum(Re(eigen_values) <= 0)) == 0)
      valid_matrix_found <- 1
  }
  A
}

drift <- make_drift_matrix(d = 20, s = 6)
drift
eigen(drift)
