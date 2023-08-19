# test_MC_WUC.R
#
# Hua Zhong
# Created: Aug 11, 2019
# Modified:
#   March 13, 2020. MS. Changed the random data from uniform integer
#     to standard normal to avoid examples with multiple optimal
#     solutions occurring to integers.
#
#   April 7, 2022. MS. Added "Test number-of-clusters selection".

library(testthat)
library(Ckmeans.1d.dp)
context("MDW_Ckmeans")

test_that("Test number-of-clusters selection", {
  # Example 1.
  n <- c(20, 20, 20)
  x <- c(rnorm(n[1], mean=-6), rnorm(n[2], mean=0), rnorm(n[3], mean=6))
  ks <- 3
  res <- Ckmeans.1d.dp(x, k=ks)
  if(0) {
    plot(x, res$cluster)
    expect_equal(length(res$size), 3)
  }

  Y <- matrix(c(
    rep(c(1,0,0), times=n[1]),
    rep(c(0,1,0), times=n[2]),
    rep(c(0,0,1), times=n[3])
  ), byrow=TRUE, nrow=length(x))

  if(0) {
    o <- order(x)
    x <- x[o]
    Y <- Y[o, ]
  }

  res <- MultiChannel.WUC(x = x, y = Y, k = ks)

  plot(x, res$cluster)
  expect_equal(length(res$size), 3)

  # Example 2.
  n <- 100
  x <- rnorm(n)
  res <- Ckmeans.1d.dp(x, k=1:10)
  expect_equal(length(res$size), 1)

  Y <- t(sapply(x, function(val) {
    if(val > 0) {c(1,0)}
    else {c(0,1)}
  }))

  res <- MultiChannel.WUC(x = x, y = Y, k = 2) # 1:10)
  expect_equal(length(res$size), 2)
})

test_that("Test 1D weights comparing Ckmeans and MDW_Ckmeans", {
  total <- 10000

  i <- 0
  for(j in c(1:total)){

    # X <- sort(sample(x = c(1:100), size = 20, replace = TRUE)) # MS March 13, 2020
    # X <- sort(sample(x = rnorm(100), size = 20, replace = TRUE)) # MS March 13, 2020

    # MS April 9, 2022: It is no longer necessary to have X sorted:
    X <- sample(x = rnorm(100), size = 20, replace = TRUE)

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

    # MS April 9, 2022: It is no longer necessary to have X sorted:
    # Replay the following line
    # X <- sort(sample(x = c(1:100), size = 20, replace = TRUE))
    # by
    X <- sample(x = c(1:100), size = 20, replace = TRUE)

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
