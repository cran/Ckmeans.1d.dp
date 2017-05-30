library(testthat)
library(Ckmeans.1d.dp)
context("Cksegs.1d.dp")

test_methods <- c("quadratic")

test_that("k-segs", {
  for(method in test_methods) {

    y <- c(-1, 2, 4, 5, 6)

    result <- Cksegs.1d.dp(y, 3, method=method)
    expect_equal(result$size, c(1,1,3))
    expect_equal(result$cluster, c(1,2,3,3,3))
    expect_equal(result$centers, c(1, 2, 4))
    expect_equal(result$withinss, c(0,0,2))

    y <- c(-.9, 1, 1.1, 1.9, 2, 2.05)
    result <- Cksegs.1d.dp(y, c(1,6), method=method)
    expect_equal(result$size, c(1,2,3))
    expect_equal(result$centers, c(1, 2.5, 5))
  }
})

