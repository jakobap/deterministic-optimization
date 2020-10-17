library(CVXR)
library(psych) 

# =====================
my_sign <- function(n){
  if((n %% 2) == 0) {
    return(-1)
  } else {
    return(1)
  }
}

# =====================
my_sign_perm <- function(i, j){
  if(i==j){return(1)}
  else(return(-1))
}

# =====================
permute<- function(v, i, j){
  aux <- v[i]
  v[i] <- v[j]
  v[j] <- aux
  v
}

# =====================
my_det_0 <- function(K){
  if (is.matrix(K)){
    return(det(K))
  }else{
    return(K)
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
        c <- my_det_0(K[idx_i, idx_j])
      }
      det_ <- det_+a*c
  }
  det_
}

my_det_1_par <- function(K){
  n <- dim(K)[1]
  if (is.array(K[c(2:n), c(2:n)])){
    return(det(K[c(2:n), c(2:n)]))
  }else{
    return(K[c(2:n), c(2:n)])
  }
}

my_det_1_fun <- function(K, x){
  n <- dim(K)[1]

  det_ <- x*my_det_0(K[c(2:n), c(2:n)])
  
  for(i in 2:n){
    idx_i <- c(1:n)
    idx_j <- c(2:n)
    a <- my_sign(i)*K[idx_i[i], 1]
    idx_i <- idx_i[-i]
    c <- my_det_0(K[idx_i, idx_j])
    det_ <- det_+a*c
  }
  det_
}


# =====================
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
      c <- my_det_0(K[idx_i, idx_j])
      det_ <- det_+a*b*c
    }
  }
  det_
}

my_det_2_par <- function(K){
  det_ <- 0
  
  n <- dim(K)[1]
  for(i in 1:n){
    for(j in 1:(n-1)){
      idx_i <- c(1:n)
      idx_j <- c(3:n)
        if (i == 1 && j == 1){
        det_ <- det_ + 2*K[1, 1]*my_det_0(K[c(3:n), c(3:n)])
      }else{
        a <- my_sign(i)*K[idx_i[i], 1]
        b <- my_sign(j)*K[idx_i[-i][j], 2]
        idx_i <- idx_i[-i][-j]
        if (1 %in% idx_i || 2 %in% idx_i){
          c <- my_det_0(K[idx_i, idx_j])
          det_ <- det_+a*b*c
        }
      }
    }
  }
  det_
}

my_det_2_fun <- function(K, x){
  det_ <- 0
  
  n <- dim(K)[1]
  for(i in 1:n){
    for(j in 1:(n-1)){
      idx_i <- c(1:n)
      idx_j <- c(3:n)
      if (i == 1 && j == 1){
        det_ <- det_ + x^2*my_det_0(K[c(3:n), c(3:n)])
      }else{
        a <- my_sign(i)*K[idx_i[i], 1]
        b <- my_sign(j)*K[idx_i[-i][j], 2]
        idx_i <- idx_i[-i][-j]
        c <- my_det_0(K[idx_i, idx_j])
        if (1 %in% idx_i){ 
          det_ <- det_+my_sign(i)*x*b*c
        }else if (2 %in% idx_i){
          det_ <- det_+a*my_sign(j)*x*c
        }else
          det_ <- det_+a*b*c
      }
    }
  }
  det_
}

# =====================
my_det_ij <- function(K, i, j, x=NULL){
  if (!is.null(x)){
    K[i, j] <- x
    K[j, i] <- x
  }
  n <- dim(K)[1]
  idx_i <- c(1:n)
  idx_j <- c(1:n)
  if (i !=j){
    idx_i <- permute(idx_i, i, 1)
    idx_j <- permute(idx_j, j, 1)
    j2 <- match(j, idx_i)
    i2 <- match(i, idx_j)
    idx_i <- permute(idx_i, j2, 2)
    idx_j <- permute(idx_j, i2, 2)
    return(my_det_2(K[idx_i, idx_j])*my_sign_perm(i, 1)*my_sign_perm(j, 1)*my_sign_perm(i2, 2)*my_sign_perm(j2, 2))
  }else{
    idx_i <- permute(idx_i, i, 1)
    idx_j <- permute(idx_j, j, 1)
    return(my_det_1(K[idx_i, idx_j])*my_sign_perm(i, 1)*my_sign_perm(j, 1))
  }
}

my_det_ij_par <- function(K, i, j){
  n <- dim(K)[1]
  idx_i <- c(1:n)
  idx_j <- c(1:n)
  if (i !=j){
    idx_i <- permute(idx_i, i, 1)
    idx_j <- permute(idx_j, j, 1)
    j2 <- match(j, idx_i)
    i2 <- match(i, idx_j)
    idx_i <- permute(idx_i, j2, 2)
    idx_j <- permute(idx_j, i2, 2)
    return(my_det_2_par(K[idx_i, idx_j])*my_sign_perm(i, 1)*my_sign_perm(j, 1)*my_sign_perm(i2, 2)*my_sign_perm(j2, 2))
  }else{
    idx_i <- permute(idx_i, i, 1)
    idx_j <- permute(idx_j, j, 1)
    return(my_det_1_par(K[idx_i, idx_j])*my_sign_perm(i, 1)*my_sign_perm(j, 1))
  }
}

my_det_ij_fun <- function(K, i, j, x){
  n <- dim(K)[1]
  idx_i <- c(1:n)
  idx_j <- c(1:n)
  if (i !=j){
    if((i==1 &j==2)|(i==2&j==1)){
      idx_i <- permute(idx_i, 1, 2)
      return(-1*my_det_2_fun(K[idx_i, idx_j], x))
    }else{
      idx_i <- permute(idx_i, i, 1)
      idx_j <- permute(idx_j, j, 1)
      j2 <- match(j, idx_i)
      i2 <- match(i, idx_j)
      idx_i <- permute(idx_i, j2, 2)
      idx_j <- permute(idx_j, i2, 2)
      return(my_det_2_fun(K[idx_i, idx_j], x)*my_sign_perm(i, 1)*my_sign_perm(j, 1)*my_sign_perm(i2, 2)*my_sign_perm(j2, 2))
    }
  }else{
    idx_i <- permute(idx_i, i, 1)
    idx_j <- permute(idx_j, j, 1)
    return(my_det_1_fun(K[idx_i, idx_j], x)*my_sign_perm(i, 1)*my_sign_perm(j, 1))
  }
}

# =====================
my_trace_ij_par <- function(S, K, i, j){
  if (i == j){
    return(S[i, j])
  }else{
    return(S[i, j]+S[j, i])
  }
}

my_trace_ij_fun <- function(S, K, i, j, x){
  t_ <- tr(S %*% K)
  if (i == j){
    t_ <- t_ - S[i, j]*K[i, j] + S[i, j]*x
  }else{
    t_ <- t_ - S[j, i]*K[i, j] + S[j, i]*x - S[i, j]*K[j, i] + S[i, j]*x
  }
  t_
}

my_trace_ij <- function(S, K, i, j, x=NULL){
  if (!is.null(x)){
    K[i, j] <- x
    K[j, i] <- x
  }
  tr(S%*%K)
}

# =====================
my_fun_fun <- function(S, K, i, j, x){
  log(my_det_ij_fun(K, i, j, x)) - my_trace_ij_fun(S, K, i, j, x)
}

my_fun <- function(S, K, i, j, x= NULL){
  if (!is.null(x)){
    K[i, j] <- x
    K[j, i] <- x
  }
  log(det(K)) - tr(S%*%K)
}
  
my_fun_par <- function(S, K, i, j, det_=NULL){
  if (is.null(det_)){
    det_ = det(K)
  }
  my_det_ij_par(K, i, j)/det_ - my_trace_ij_par(S, K, i, j)
}

find_ij <- function(S, K){
  n <- dim(K)[1]
  
  the_i <- 0
  the_j <- 0
  par  <- -1
  det_  <- det(K)
  
  for (i in 1:n){
    for (j in i:n){
      p <- abs(my_fun_par(S, K, i, j, det_))
      if (p > par){
        the_i <- i
        the_j <- j
        par <- p
      }
    }
  }
  c(the_i, the_j)
}

#==========
# Testing
# =========
print('Testing Determinant')
m <- matrix(rexp(36, rate=.5), ncol=6)
t <- det(m)

for (i in 1:6){
  for (j in 1:6){
    t_ <- my_det_ij(m, i, j)
    
    if ((t-t_)^2>.01) {
      print(i)
      print(j)
      print(t)
      print(t_)
    }
  }
}

print('Testing fun_fun')
S <- matrix(rnorm(16), ncol=4)
set.seed(1)
nums <- c()
for (i in 1:9){
  nums[i] <- sample(0:3, 1)
}

K <- matrix(nums, ncol=3)
K <- K+t(K)

for (i in 1:3){
  for (j in 1:3){
    d1 <- my_det_ij(K, i, j)
    d2 <- my_det_ij_fun(K, i, j, K[i, j])
    if (is.nan(d2)){
      print('NAN')
      print(i)
      print(j)
    }
    if (abs(d1 - d2) > 1){
      print(i)
      print(j)
      print(my_det_ij(K, i, j))
      print(my_det_ij_fun(K, i, j, K[i, j]))
      print('que')
    }
  }
}

print('done ;-)')
#rm(m, t, t_, i, j)
