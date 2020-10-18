source('new_life_nice_easy.r')
source('not_so_easy.r')
library(MASS)

#set.seed(1)

performance_test <- function(lower, upper){

  test_results_df <- data.frame(matrix(ncol = 3, nrow = 0)) ## Dataframe for test results
  colnames(test_results_df) <- c("dimensions", "solving_time_e2", "solving_time_e3")

  for (n in lower:upper) {
    S <- matrix(rnorm(n * n), ncol = n)
    S <- S %*% t(S)

    start_time <- Sys.time()
    res_e2 <- exercise_2(S, pre = 1e-2, max_iter = 1000)
    end_time <- Sys.time()
    time_e2 <- end_time - start_time

    start_time <- Sys.time()
    res_e3 <- exercise_3(S, pre = 1e-2, max_iter = 100)
    end_time <- Sys.time()
    time_e3 <- end_time - start_time


    df <- data.frame(n, time_e2, time_e3) ## result index 1 expected as time
    colnames(df) <- c("dimensions", "solving_time_e2", "solving_time_e3")
    test_results_df <- rbind(test_results_df, df)

  }
  return(test_results_df)
}

test_results_df <- performance_test(2, 15)

ggplot(test_results_df, aes(x=dimensions)) +
  geom_line(aes(y = as.numeric(solving_time_e2)), color = "red") +
  geom_line(aes(y = as.numeric(solving_time_e3)), color = "green") +
  labs(x="Dimensions of S", y="Solving Time in Seconds")