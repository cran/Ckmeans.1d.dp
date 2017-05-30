library(testthat)
library(Ckmeans.1d.dp)
context("Ckmedian.1d.dp")

test_methods <- c("linear", "loglinear", "quadratic")

test_that("Unweighted k-median", {
  for(method in test_methods) {

    x <- c(-1, 2, 4, 5, 6)

    result <- Ckmedian.1d.dp(x, 3, method)
    expect_equal(result$size, c(1,1,3))
    expect_equal(result$cluster, c(1,2,3,3,3))
    expect_equal(result$centers, c(-1, 2, 5))
    expect_equal(result$withinss, c(0,0,2))

    x <- c(-.9, 1, 1.1, 1.9, 2, 2.05)
    result <- Ckmedian.1d.dp(x, c(1,6), method)
    expect_equal(result$size, c(1,2,3))
    expect_equal(result$centers, c(-0.9, 1, 2))
  }
})

