# test_Ckmeans.1d.dp.R
#
# Joe Song
# Created: May 3, 2016
# Updated: May 19, 2016

library(testthat)
library(Ckmeans.1d.dp)
context("Checking on several examples")

test_that("Given the number of clusters", {

  x <- c(-1, 2, -1, 2, 4, 5, 6, -1, 2, -1)

  result <- Ckmeans.1d.dp(x, 3)
  expect_equal(result$size, c(4,3,3))

  cluster.truth <- c(1,2,1,2,3,3,3,1,2,1)
  expect_equal(result$cluster, cluster.truth)

  centers.truth <- c(c(-1, 2, 5))
  expect_equal(result$centers, centers.truth)
  withinss.truth <- c(0,0,2)
  expect_equal(result$withinss, withinss.truth)

  totss.truth <- sum(scale(x, scale=FALSE)^2)
  expect_equal(result$totss, totss.truth)
  expect_equal(result$tot.withinss, 2)
  expect_equal(result$betweenss, totss.truth - sum(withinss.truth))

})

test_that("n<=k", {

  x <- c(3, 2, -5.4, 0.1);
  res <- Ckmeans.1d.dp(x, 4)

  cluster.truth <- c(4, 3, 1, 2)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(-5.4, 0.1, 2, 3)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(0, 0, 0, 0)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(1, 1, 1, 1)
  expect_equal(res$size, size.truth)

})

test_that("k==2", {

  x <- 1:10
  res <- Ckmeans.1d.dp(x, 2)

  cluster.truth <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(3, 8)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(10, 10)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(5, 5)
  expect_equal(res$size, size.truth)

})

test_that("k==1", {

  x <- c(-2.5, -2.5, -2.5, -2.5)
  res <- Ckmeans.1d.dp(x, 1)

  cluster.truth <- c(1, 1, 1, 1)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(-2.5)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(0)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(4)
  expect_equal(res$size, size.truth)

  x <- rep(1, 100)
  result <- Ckmeans.1d.dp(x, 1)
  expect_equal(result$size, 100)

})

test_that("n==10, k==3", {

  x <- c(3, 3, 3, 3, 1, 1, 1, 2, 2, 2)
  res <- Ckmeans.1d.dp(x, 3)

  cluster.truth <- c(3, 3, 3, 3, 1, 1, 1, 2, 2, 2)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(1, 2, 3)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(0, 0, 0)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(3, 3, 4)
  expect_equal(res$size, size.truth)

})


test_that("n==14, k==8", {

  x <- c(-3, 2.2, -6, 7, 9, 11, -6.3, 75, 82.6, 32.3, -9.5, 62.5, 7, 95.2)
  res <- Ckmeans.1d.dp(x, k=8)

  cluster.truth <- c(2, 2, 1, 3, 3, 3, 1, 6, 7, 4, 1, 5, 3, 8)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(-7.266666667, -0.400000000, 8.500000000, 32.300000000,
                     62.500000000, 75.000000000, 82.600000000, 95.200000000)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(7.526666667, 13.520000000, 11.000000000, 0.000000000,
                      0.000000000, 0.000000000, 0.000000000, 0.000000000)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(3, 2, 4, 1, 1, 1, 1, 1)
  expect_equal(res$size, size.truth)

})

test_that("Estimating k example set 1", {

  x <- c(.9, 1, 1.1, 1.9, 2, 2.1)
  result <- Ckmeans.1d.dp(x, c(1,6))
  expect_equal(result$size, c(3,3))

  x <- rev(x)
  result <- Ckmeans.1d.dp(x, c(1,6))
  expect_equal(result$size, c(3,3))

  x <- 1:10
  result <- Ckmeans.1d.dp(x, k=c(1,10))
  expect_equal(result$size, 10)

})

test_that("Estimating k example set 2", {

  x <- c(3.5, 3.6, 3.7, 3.1, 1.1, 0.9, 0.8, 2.2, 1.9, 2.1)
  res <- Ckmeans.1d.dp(x, k=c(2,5))

  cluster.truth <- c(3, 3, 3, 3, 1, 1, 1, 2, 2, 2)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(0.933333333333, 2.066666666667, 3.475000000000)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(0.0466666666667, 0.0466666666667, 0.2075000000000)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(3, 3, 4)
  expect_equal(res$size, size.truth)

  x <- cos((-10:10))
  res <- Ckmeans.1d.dp(x)
  # format(res, digits=10)

  cluster.truth <- c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(-0.6592474631, 0.6751193405)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(1.0564793100, 0.6232976959)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(12, 9)
  expect_equal(res$size, size.truth)

  x <- dgamma(seq(1,10, by=0.5), shape=2, rate=1)
  res <- Ckmeans.1d.dp(x)
  # format(res, digits=10)

  cluster.truth <- c(3, 3, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(0.01702193495, 0.15342151455, 0.32441508262)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(0.006126754998, 0.004977009034, 0.004883305120)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(13, 3, 3)
  expect_equal(res$size, size.truth)

})
