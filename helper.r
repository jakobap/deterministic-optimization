source('new_life_nice_easy.r')
source('not_so_easy.r')
library(MASS)

#set.seed(1)

n <- 9     ## Dimension of matrix
S <- matrix(rnorm(n*n), ncol=n)
S <- S%*%t(S)

n <- dim(S)[1]
K <- diag(n)

print('exercise_2')
Kr_2 <- exercise_2(S, pre=1e-2, max_iter = 1000)
print(Kr_2)
print('exercise_3')
Kr_3 <- exercise_3(S, pre=1e-2, max_iter = 100)
print(Kr_3)