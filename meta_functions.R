
my_sign <- function(n){
  if((n %% 2) == 0) {
    return(-1)
  } else {
    return(1)
  }
}

my_det_1 <- function(K){
  
  det_ <- 0
  
  n <- dim(K)[1]
  for(i in 1:n){
      idx_i <- c(1:n)
      idx_j <- c(2:n)
      a <- my_sign(i)*K[idx_i[i], 1]
      idx_i <- idx_i[-i]
      if (n == 2){
        c <- K[idx_i, idx_j]
      }else{
        c <- det(K[idx_i, idx_j])
      }
      det_ <- det_+a*c
  }
  det_
}

my_det_2 <- function(K){
  
  det_ <- 0
  
  n <- dim(K)[1]
  for(i in 1:n){
    for(j in 1:(n-1)){
      idx_i <- c(1:n)
      idx_j <- c(3:n)
      a <- my_sign(i)*K[idx_i[i], 1]
      b <- my_sign(j)*K[idx_i[-i][j], 2]
      idx_i <- idx_i[-i][-j]
      c <- det(K[idx_i, idx_j])
      det_ <- det_+a*b*c
    }
  }
  det_
}

permute<- function(v, i, j){
  aux <- v[i]
  v[i] <- v[j]
  v[j] <- aux
  v
}

my_det_ij <- function(K, i, j){
  n <- dim(K)[1]
  idx_i <- c(1:n)
  idx_j <- c(1:n)
  idx_i <- permute(idx_i, i, 1)
  idx_j <- permute(idx_j, j, 1)
  idx_i <- permute(idx_i, j, 2)
  idx_j <- permute(idx_j, i, 2)
  
  my_det_12(K[idx_i, idx_j])
}