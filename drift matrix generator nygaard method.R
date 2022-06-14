make_drift_matrix <- function(d, sparsity){
  A <- diag(d)
  n <- floor((d^2)*sparsity) - d
  values <- rbern(n = n, prob = 0.5)
  values <- (values - 0.5)*2
  
  #browser()
  values_inserted <- 0
  while(values_inserted < n){
    coord_found <- 0
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
  }
  A
}

drift_matrix <- make_drift_matrix(d = 5, sparsity = 0.3)
eigen(drift_matrix)$values
View(drift_matrix)
