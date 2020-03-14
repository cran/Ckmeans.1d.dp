# test_MDW_Ckmeans.1d.dp.R
#
# Hua Zhong
# Created: Aug 11, 2019
# Modified:
#   March 13, 2020. MS. Changed the random data from uniform integer
#     to standard normal to avoid examples with multile optimal
#     solutions occuring to integers.

library(testthat)
library(Ckmeans.1d.dp)
context("MDW_Ckmeans")

test_that("Test 1D weights comparing Ckmeans and MDW_Ckmeans", {
  total <- 10000

  i <- 0
  for(j in c(1:total)){

    # X <- sort(sample(x = c(1:100), size = 20, replace = TRUE)) # MS March 13, 2020
    X <- sort(sample(x = rnorm(100), size = 20, replace = TRUE)) # MS March 13, 2020

    Y <- matrix(sample(x = c(1:100), size = 20, replace = TRUE),
                ncol=1, nrow=length(X))

    cluster.num <- c(5:10)

    res.1 <- MultiChannel.WUC(x = X, y = Y, k = cluster.num)

    res.2 <- suppressWarnings(
      Ckmeans.1d.dp(x = X, y = as.numeric(Y), k = cluster.num,
      method = "linear", estimate.k = "BIC"))

    # expect_equal(res.1$cluster, res.2$cluster)
    # expect_equal(res.1$centers, res.2$centers)
    # expect_equal(res.1$withinss, res.2$withinss)
    # expect_equal(res.1$size, res.2$size)
    # #expect_equal(res.1$BIC, res.2$BIC)

    if(all(c(all(res.1$cluster == res.2$cluster),
             all(res.1$centers - res.2$centers < 1e-13),
             all(res.1$withinss - res.2$withinss < 1e-13),
             all(res.1$size - res.2$size < 1e-13),
             all(res.1$BIC - res.2$BIC < 1e-13)
      ))) {
      i <- i + 1
    }else{
      # stop()
    }
  }
  # print(i)
  expect_equal((total - i) / total * 1000 , 0, tolerance = 1)

})

test_that("Test to compare linear and quadratic algorithms", {
  i <- 0
  total <- 10000
  for(j in c(1:total)){
    X <- sort(sample(x = c(1:100), size = 20, replace = TRUE))
    Y <- matrix(sample(x = c(1:100), size = 20, replace = TRUE), ncol=1, nrow=length(X))

    cluster.num <- c(5:10)

    res.1 <- MultiChannel.WUC(x = X, y = Y, k = cluster.num)

    res.2 <- MultiChannel.WUC(x = X, y = Y, k = cluster.num)

    test <- all(all(res.1$cluster == res.2$cluster),
                all(res.1$centers == res.2$centers),
                all(res.1$withinss == res.2$withinss),
                all(res.1$size == res.2$size),
                all(res.1$BIC == res.2$BIC))
    i <- i + 1
  }
  # print(i)
  expect_equal((total - i) / total * 1000 , 0, tolerance = 1)

})
