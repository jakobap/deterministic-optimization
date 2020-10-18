source('new_way.r')

optimize_f_ij <- function(S, K, i, j){
  f <- poli_fun(S, K, i, j)
  
  if (i == j){
    x0 <- (K[i, j])^.5
  }else{
    x0 <- (-K[i, j])^.5
  }
  
  
  result <- optim(x0, f, control=list(fnscale=-1, reltol=1e-8))
  
  if(i == j){
    K[i,j] <- (result$par)^2
  }else{
    K[i, j] <- -(result$par)^2
    K[j, i] <- -(result$par)^2
  }
  K
}

exercise_3 <- function(S, pre = 1e-1, max_iter = 100){
  
  n <- dim(S)[1]
  K <- diag(n)
  K <- K - 0.1*matrix(1, nrow = n, ncol = n) 
  # K <- S
  # 
  # for (i in 1:n){
  #   for (j in 1:n){
  #     if (i!= j){
  #       if (K[i, j]>0){
  #         K[i, j] <- -K[i, j]
  #       }
  #     }
  #   }
  # }
  
  last_value <- Inf
  
  for (i in 1:n){
    for (j in i:n){
      aux <- optimize_f_ij(S, K, i, j)
      K <- aux
    }
  }
  
  for (mi in 1:max_iter){
    ij = find_ij(S, K)
    
    i <- ij[1]
    j <- ij[2]
    
    aux <- optimize_f_ij(S, K, i, j)
    K <- aux
    
    # if (abs(last_value - aux[2]$value)< pre){
    #   break
    # }else{
    #   last_value <- aux[2]$value
    # }
  }
  round(K, 2)
}