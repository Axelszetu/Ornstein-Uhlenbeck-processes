#Test script for 20X20 drift matrix
d <- 15
t <- 100
N <- 10000
dt <- t/N
sparsity <- 0.3

A <- make_drift_matrix(d = d, sparsity = sparsity)

#C_infty <- C_infty_calculator(A, ds = 0.1)

C_h <- C_hat(A, dt = dt)

process <- SimOU(A0 = A, t = t, N = N, burn = 10000)

A_mle <- OU_MLE_analytical(X = process, dt = dt)

A_Lasso <- OU_Lasso(X = process, dt = dt, lambda = 0.05)
A_Lasso

A_Dantzig <- OU_Dantzig(X = process, dt, lambda = 0.005)
A_Dantzig

#Below we test the cross validator
#We note that it takes about 11 seconds to run OU_Lasso with current parameters.
#This gives us an idea of how long we should let the cross validator run before we start to worry.

bases <- (5:-8)
lambdas <- exp(bases)
lambda_L <- cross_validator(process, dt = dt, f = OU_Lasso, pars = lambdas, split = 0.8)
A_Lasso_CV <- OU_Lasso(X = process, dt = dt, lambda = lambda_L)
A_Lasso_CV #Looks good to me

lambda_D <- cross_validator(process, dt = dt, f = OU_Dantzig, pars = lambdas, split = 0.8)
A_Dantzig_CV <- OU_Dantzig(X = process, dt = dt, lambda = lambda_D)
A_Dantzig

#Below are tests for OU_Lasso with F-norm
A_Lasso_F <- OU_Lasso(X = process, dt = dt, lambda = lambda_L, penalization = "F")
A_Lasso_F #This looks good too

ggplot(data = A_hat_melt, mapping = aes(x = Column, y = Row, fill = value)) + geom_tile(color = "white") +
  scale_fill_gradient2(low = rgb(0,0,219,maxColorValue = 255), high = rgb(219,0,0,maxColorValue = 255), mid = "white",
                       midpoint = 0, space = "Lab", 
                       name="Drift") +
  scale_y_reverse() + coord_fixed() + ggtitle(A_hat)

A_D <- OU_Dantzig(X = process, dt = dt, lambda = 0.05)
A_D_melt <- melt(t(A_D), varnames = c("Column", "Row"))
ggplot(data = A_D_melt, mapping = aes(x = Column, y = Row, fill = value)) + geom_tile(color = "white") +
  scale_fill_gradient2(low = rgb(0,0,219,maxColorValue = 255), high = rgb(219,0,0,maxColorValue = 255), mid = "white",
                       midpoint = 0, limit = limit, space = "Lab", 
                       name="Drift") +
  scale_y_reverse() + coord_fixed() + ggtitle("A_D")

normalization <- sum(abs(A))
MLE_error <- sum(abs(A-A_mle))
Dantzig_error <- sum(abs(A-A_D))
Lasso_CV_error <- sum(abs(A-A_L_CV))
Lasso_error <- sum(abs(A-A_Lasso))
A_Lasso
