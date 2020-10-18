library(psych) 


poli_det <- function(K, i, j){
  if (i==j){
    # Calculate
    K[i, j] <- 0
    V0 <- det(K)
    K[i, j] <- 1
    V1 <- det(K)
    
    # Solve
    c <- V0
    a <- V1 - c
    
    f <- function(x){
      a*x + c
    }
  }else{
    # Calculate
    K[i, j] <- 0
    K[j, i] <- 0
    V0 <- det(K)
    K[i, j] <- 1
    K[j, i] <- 1
    V1 <- det(K)
    K[i, j] <- -1
    K[j, i] <- -1
    V_1 <- det(K)
    
    # Solve
    c <- V0
    a <- (V1 + V_1 - 2*c)/2
    b <- V1 - c - a
    
    f <- function(x){
      a*x^2+b*x++c
    }
  }
  f
}

poli_trace <- function(S, K, i, j){
  if(i==j){
    # Calculate
    K[i, j]<- 0
    V0 <- tr(S%*%K)
    K[i, j]<- 1
    V1 <- tr(S%*%K)
    
    # Solve
    b <- V0
    a <- V1 - b
    
    f <- function(x){
      a*x+b
    }
  }else{
    # Calculate
    K[i, j] <- 0
    K[j, i] <- 0
    V0 <- tr(S%*%K)
    K[i, j] <- 1
    K[j, i] <- 1
    V1 <- tr(S%*%K)

    # Solve
    b <- V0
    a <- V1 - b
    
    f <- function(x){
      a*x+b
    }
  }
}

poli_fun <- function(S, K, i, j){
  pd <- poli_det(K, i, j)
  pt <- poli_trace(S, K, i, j)
  if (i==j){
    function(x){
      log(pd(x^2))-pt(x^2)
    }
  }else{
    function(x){
      log(pd(-x^2))-pt(-x^2)
    }
  }
}

fun <- function(S, K){
  log(det(K)) - tr(S%*%K)
}

poli_par <- function(f, x0){
  eps = 0.0001
  (f(x0+eps)-f(x0-eps))/(2*eps)
}

poli_par_ij <- function(S, K, i, j){
  f <- poli_fun(S, K, i, j)
  if (i==j){
    par <- poli_par(f, (K[i, j])^.5)
  }else{
    par <- poli_par(f, (-(K[i, j]))^.5)
  }
  par
}

find_ij <- function(S, K){
  n <- dim(K)[1]
  
  the_i <- 0
  the_j <- 0
  the_v <- -Inf
  
  for (i in 1:n){
    for (j in i:n){
      if (i == j){
        new_v <- abs(poli_par(poli_fun(S, K, i, j), (K[i, j])^.5))
      }else{
        new_v <- abs(poli_par(poli_fun(S, K, i, j), (-K[i, j])^.5))
      }
      
      if (new_v > the_v){
        the_v <- new_v
        the_i <- i
        the_j <- j
      }
    }
  }
  c(the_i, the_j)
}
#===============================
# testing
n <- 5     ## Dimension of matrix
S <- matrix(rnorm(n*n), ncol=n)
S <- S*t(S)

K <- matrix(rnorm(n*n), ncol=n)
K <- K*t(K)

print('poli 1, 1')
f <- poli_det(K, 1, 1)

K1 <- K
K1[1,1] <- 0
print(f(0))
print(det(K1))

K1[1,1] <- 3.14
print(f(3.14))
print(det(K1))


print('poli 2, 3')
f <- poli_det(K, 2, 3)

K1 <- K
K1[2,3] <- 0
K1[3,2] <- 0
print(f(0))
print(det(K1))

K1[2,3] <- 3.14
K1[3,2] <- 3.14
print(f(3.14))
print(det(K1))

print('poli_trace 1, 1')
g <- poli_trace(S, K, 1, 1)

K1 <- K
K1[1,1] <- 0
print(g(0))
print(tr(S%*%K1))

K1[1,1] <- 3.14
print(g(3.14))
print(tr(S%*%K1))


print('poli_trace 2, 3')
g <- poli_trace(S, K, 2, 3)

K1 <- K
K1[2,3] <- 0
K1[3,2] <- 0
print(g(0))
print(tr(S%*%K1))

K1[2,3] <- 3.14
K1[3,2] <- 3.14
print(g(3.14))
print(tr(S%*%K1))
