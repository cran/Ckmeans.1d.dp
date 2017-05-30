# print.R -- define S3 print functions for Ckmeans.1d.dp,
#            Ckmedian.1d.dp, Cksegs.1d.dp
#
# Joe Song
#
# Created: May 29, 2017.

## print method for Ckmeans.1d.dp object
print.Ckmeans.1d.dp <- function(x, ...)
{
  with(x, {
    cat("\nCluster centers:\n")
    print(centers, ...)
    cat("\nCluster index associated with each element:\n")
    print(cluster, ...)
    cat("\nWithin-cluster sum of squares:\n")
    print(withinss, ...)
    if(length(size) > 1) {
      cat("Ckmeans.1d.dp returns", length(size), "optimal clusters of sizes",
          paste(size, collapse=", "), "\n")
      cat("  minimum total within-cluster sum of squares:", tot.withinss, "\n")
      cat("  maximum between-cluster sum of squares:", betweenss, "\n")
      cat("  total sum of squares of input vector:", totss, "\n")
      if(totss > 0) {
        cat("  maximum (between-SS / total-SS):",
            round(1000 * betweenss / totss)/10, "%\n")
      }
    } else {
      cat("Ckmeans.1d.dp returns one cluster containing all elements.\n")
    }
  })

  cat("\nAvailable components:\n")
  print(names(x))
  invisible(x)
}

print.Ckmedian.1d.dp <- function(x, ...)
{
  with(x, {
    cat("\nCluster centers:\n")
    print(centers, ...)
    cat("\nCluster index associated with each element:\n")
    print(cluster, ...)
    cat("\nWithin-cluster sum of L1 distances:\n")
    print(withinss, ...)
    if(length(size) > 1) {
      cat("Ckmedian.1d.dp returns", length(size), "optimal clusters of sizes",
          paste(size, collapse=", "), "\n")
      #cat("  minimum total within-cluster sum of squares:", tot.withinss, "\n")
      #cat("  maximum between-cluster sum of squares:", betweenss, "\n")
      #cat("  total sum of squares of input vector:", totss, "\n")
      #if(totss > 0) {
      #  cat("  maximum (between-SS / total-SS):",
      #      round(1000 * betweenss / totss)/10, "%\n")
      #}
    } else {
      cat("Ckmedian.1d.dp returns one cluster containing all elements.\n")
    }
  })

  cat("\nAvailable components:\n")
  print(names(x))
  invisible(x)
}

print.Cksegs.1d.dp <- function(x, ...)
{
  with(x, {
    cat("\nSegment centers:\n")
    print(centers, ...)
    cat("\nSegment index associated with each element:\n")
    print(cluster, ...)
    cat("\nWithin-segment sum of L2 distances in Y:\n")
    print(withinss, ...)
    if(length(size) > 1) {
      cat("Cksegs.1d.dp returns", length(size), "optimal segments of sizes",
          paste(size, collapse=", "), "\n")
      #cat("  minimum total within-cluster sum of squares:", tot.withinss, "\n")
      #cat("  maximum between-cluster sum of squares:", betweenss, "\n")
      #cat("  total sum of squares of input vector:", totss, "\n")
      #if(totss > 0) {
      #  cat("  maximum (between-SS / total-SS):",
      #      round(1000 * betweenss / totss)/10, "%\n")
      #}
    } else {
      cat("Cksegs.1d.dp returns one segment containing all elements.\n")
    }
  })

  cat("\nAvailable components:\n")
  print(names(x))
  invisible(x)
}
