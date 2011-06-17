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
    cat("Ckmeans.1d.dp clustering with ", length(x$size), " clusters of sizes ",paste(x$size, collapse=", "), "\n", sep="")
    cat("\nCluster means:\n")
    print(x$centers, ...)
    cat("\nClustering vector:\n")
    print(x$cluster, ...)
    cat("\nWithin cluster sum of squares by cluster:\n")
    print(x$withinss, ...)
	cat("\nAvailable components:\n", sep="\n")
    print(names(x))
    invisible(x)
}

##Ckmeans.1d.dp : function which implement optimal one-dimensional clustering
## x is one-dimensional input vector
## k indicates cluster level
Ckmeans.1d.dp <- function( x, k )   
{
	##Check to see if k is less than 0.
	##If k is less than 0, stop and exit.
    if( k <= 0 )
    {
		stop ("Can NOT classify vector into 0 or less cluster\n")
    }
	 
	##Check to see if cluster level bigger than the unique number of the input vector
	##If k is bigger than the unique number of the input vector, 
	##force k set to the number of unique number of input.
    if(length(unique(x)) < k)
    {
		print ("Number of clusters is bigger than the unique number of the input vector\n k is set to the number of unique number of input.")
		k <- length(unique(x))
    }
	
	##Form data which will be passed to external C++ function.
    clusters <- c(1:length(x))
    center <- c(1:k)
    withinss <- c(1:k)
    size <- c(1:k)

	#Call external C++ function
    result <- .C("Ckmeans_1d_dp", PACKAGE="Ckmeans.1d.dp", n=as.double(x), length=as.integer(length(x)), level=as.integer(k),cluster=as.integer(clusters), centers=as.double(center), withinss=as.double(withinss), size=as.integer(size))
    
	r = structure(list(cluster = result$cluster, centers = result$centers,
                   withinss = result$withinss,size = result$size),
			class = "Ckmeans.1d.dp")
			
	return (r)
}##end of Ckmeans.1d.dp()