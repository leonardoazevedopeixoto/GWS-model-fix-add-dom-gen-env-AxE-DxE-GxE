## Given the vectors of sire and dam, return additive relationship matrix 'A'.

## Arguments
## s: a vector of sire
## d: a vector of dam 

## Note: Unknown parents should be coded as zero.

## Literature: Henderson, C. R. 1976. Simple Method for Computing the Inverse of a Numerator Relationship Matrix Used in Prediction of Breeding Values. Biometrics 32:69-83.

## Author: Gota Morota <morota at wisc dot edu>
## Create: 16-Apr-2009
## Last-Modified: 1-Apr-2010
## License: GPLv3 or later

`createA` <-
function(s, d){
	if (nargs()==1){
		stop("sire vector and dam vector are required")
	}
	
	if (length(s) != length(d)){
		stop("size of the sire vector and dam vector are different!")
	}
	 
	n <- length(s)
	N <- n + 1
	A <- matrix(0, ncol=N, nrow=N)
		
	s <- (s == 0)*(N) + s
	d <- (d == 0)*N + d
				
	for(i in 1:n){
		
		A[i,i] <- 1 + A[s[i], d[i]]/2
			
		for(j in (i+1):n){
			if (j > n) break
			A[i,j] <- ( A[i, s[j]] + A[i,d[j]] )/2
		  	A[j,i] <- A[i,j] 	
		}			
	}
		
		return(A[1:n, 1:n])
	
}

