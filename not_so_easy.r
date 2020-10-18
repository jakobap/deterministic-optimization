library(CVXR)
library(MASS)


exercise_2 <- function(S, pre = 1e-1, max_iter = 100){
  n <- dim(S)[1]
  
  K  <- Variable(n, n, PSD = TRUE)     ## Variable constrained to positive semidefinite cone
  
  obj <- Maximize(log_det(K) - matrix_trace(S %*% K))

  ix <- 1
  constr <- list()
  
  for (i in 1:n){
    for (j in 1:n){
      if (i == j){
        constr[ix] <- K[i, j] >= 0
      }else{
        constr[ix] <- K[i, j] <= 0
      }
      ix <- ix + 1
    }
  }
  
  prob <- Problem(obj,constr)
  result <- solve(prob)
  round(result$getValue(K),2)
}