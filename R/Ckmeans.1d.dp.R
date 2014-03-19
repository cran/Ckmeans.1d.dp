##   Ckmeans.1d.dp()-----Interface funcition to call C++ version kmeans.1d.dp()
##               
##   Haizhou Wang
##   Computer Science Department
##   New Mexico State University
##   hwang@cs.nmsu.edu
##
##   Created: Oct 20, 2009
##

## modelled on print methods in the cluster package
print.Ckmeans.1d.dp <- function(x, ...)
{
  if(length(x$size) > 1) {
    cat("Ckmeans.1d.dp returns ", length(x$size), " clusters, each containing ", 
        paste(x$size, collapse=", "), " elements, respectively.\n", sep="")    
  } else {
    cat("Ckmeans.1d.dp returns one cluster containing all elements.\n")        
  }
  
  cat("\nCluster means:\n")
  print(x$centers, ...)
  cat("\nCluster id of each element:\n")
  print(x$cluster, ...)
  cat("\nWithin-cluster sum of squares:\n")
  print(x$withinss, ...)
  cat("\nAvailable components:\n")
  print(names(x))
  invisible(x)
}

##Ckmeans.1d.dp : function which implement optimal one-dimensional clustering
## x is one-dimensional input vector
## k indicates cluster level
Ckmeans.1d.dp <- function( x, k=c(1,9) )   
{
  if(is.null(k)) {
    k <- 1: min( 9, length(x) )
  }
  
  if(length(k) > 1) {
    k.min <- min(k)
    k.max <- max(k)  
  } else {
    k.min <- k
    k.max <- k
  }
  
	##Check to see if k is less than 0.
	##If k is less than 0, stop and exit.
  if( k.max <= 0 ) {
		stop ("Can NOT classify vector into 0 or less cluster\n")
  }

	##Check to see if cluster level bigger than the unique number of the input vector
	##If k is bigger than the unique number of the input vector, 
	##force k set to the number of unique number of input.
  if(length(unique(x)) < k.min) {
    
		print ("Number of clusters is bigger than the unique number of the input vector\n k is set to the number of unique number of input.")
		k <- 1:length(unique(x))
  }
	
	##Form data which will be passed to external C++ function.
  clusters <- vector("integer", length(x))
  center <- vector("double", k.max)
  withinss <- vector("double", k.max)
  size <- vector("integer", k.max)

	#Call external C++ function
  result <- .C("Ckmeans_1d_dp", PACKAGE="Ckmeans.1d.dp", 
               data=as.double(x), length=as.integer(length(x)), 
               levels=as.integer(k), n.levels=as.integer(length(k)),  
               cluster=as.integer(clusters), centers=as.double(center), 
               withinss=as.double(withinss), size=as.integer(size))

  k.opt <- length(unique(result$cluster))
  
  if(k.min < k.max) {
    if (k.opt == k.min) {
      cat("WARNING: Min number of clusters used. Consider decreasing it!\n")
    } else if(k.opt == k.max) {
      cat("WARNING: Max number of clusters used. Consider increasing it!\n")
    }
  }
  
  r <- structure(list(cluster = result$cluster, centers = result$centers[1:k.opt],
                      withinss = result$withinss[1:k.opt], size = result$size[1:k.opt]),
                 class = "Ckmeans.1d.dp")
			
	return( r )
} ##end of Ckmeans.1d.dp()