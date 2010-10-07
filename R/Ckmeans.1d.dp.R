##   Ckmeans.1d.dp()-----Interface funcition to call C++ version kmeans.1d.dp()
##               
##   Haizhou Wang
##   Computer Science Department
##   New Mexico State University
##   hwang@cs.nmsu.edu
##
##   Created: Oct 20, 2009
##


##load dynamic library
##.dll under Windows
##.so  under UNIX/Linux
dyn.load( paste("src/Ckmeans.1d.dp", .Platform$dynlib.ext,sep="") )

##Ckmeans.1d.dp : function which implement optimal one-dimensional clustering
## x is one-dimensional input vector
## K indicates cluster level
Ckmeans.1d.dp <- function( x, K)   
{
	##Check to see if K is less than 0.
	##If K is less than 0, stop and exit.
    if( K <= 0 )
    {
		print ("Error! Can NOT classify vector into 0 or less cluster\n")
		return (1)
    }
	 
	##Check to see if cluster level bigger than the unique number of the input vector
	##If K is bigger than the unique number of the input vector, 
	##force K set to the number of unique number of input.
    if(length(unique(x)) < K)
    {
		K <- length(unique(x))
    }
	
	##Form data which will be passed to external C++ function.
    cluster <- c(1:length(x))
    centers <- c(1:K)
    withinss <- c(1:K)
    size <- c(1:K)

	#Call external C++ function
    result <- .C("Ckmeans_1d_dp", PACKAGE="Ckmeans.1d.dp", n=as.double(x), length=as.integer(length(x)), level=as.integer(K),clusters=as.integer(cluster), center=as.double(centers), withinss=as.double(withinss), size=as.integer(size))
    
	##Give back clustering result
	##Result will have structre like:
	##############################################
	##$n                                        ##
	##[1]  47 21 65 52 97 66 61 81  3 53        ##  Show the actural input vector
    ##                                          ##
	##$length                                   ##
	##[1] 10                                    ##  Show length of the input vector
	##                                          ##
	##$level                                    ##
	##[1] 2                                     ##  Show cluster level
	##                                          ##
	##$clusters                                 ##  Show which cluster each element belongs to
	##[1] 2 1 2 2 2 2 2 2 1 2                   ##
	##                                          ##
	##$center                                   ##  Show center value(mean) of each cluster
	##[1] 12.00 65.25                           ##
	##                                          ##
	##$withinss                                 ##  within-cluster sum of squares for each cluster
	##[1]  162.0 1933.5                         ##
	##                                          ##
	##$size                                     ##  Size of each cluster
	##[1] 2 8                                   ##
	##############################################
    return (result)
}##end of Ckmeans.1d.dp()
