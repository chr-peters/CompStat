# Name: Christian Peters

library(microbenchmark)
library(ggplot2)

# this library is used to execute the variance simulation in parallel
library(parallel)

# No 1)
#======
# i)
print(1.01 + 1.02 == 1.03)
# ii)
print(0.1*0.05/0.05 == 0.1)
# iii)
print(sqrt(2)**2 == 2)
# iv)
print(2 + 1e32 - 1e32 == 2)
# v)
print(exp(log(3)) == 3)
# vi)
print(choose(23, 2) == factorial(23) / (factorial(2)*factorial(21)))
# example 1: roots of a polynomial
# the polynomial (x-2)^2 = x^2 - 4x + 4 has 2 as the only root
print(polyroot(c(4, -4, 1))==2)
# example 2: sin and cos
# sin(x)**2 + cos(x)**2 = 1
print(sin(0.08)**2 + cos(0.08)**2 == 1)


# No 3)
#======

#' Computes the variance of a vector using the two pass algorithm.
#' 
#' @param x A numeric vector.
#' 
#' @return The variance of x.
varTwoPass <- function(x) {
  stopifnot(is.numeric(x))
  return(sum((x-mean(x))**2)/(length(x)-1))
}

#' Computes the variance of a vector using the textbook algorithm.
#' 
#' @param x A numeric vector.
#' 
#' @return The variance of x.
varTextBook <- function(x) {
  stopifnot(is.numeric(x))
  return((sum(x**2) - sum(x)**2/length(x))/(length(x)-1))
}

#' Computes the variance of a vector 1...n or an arbitrary permutation
#' thereof analytically.
#' 
#' @param x A numeric vector containing the first n natural numbers.
#' 
#' @return The variance of x.
varAnalytical <- function(x) {
  stopifnot(is.numeric(x))
  n <- length(x)
  return((n*(n+1))/12.)
}

varEval <- function() {
  # simulation parameters
  minExponent <- 1
  maxExponent <- 26
  simCount <- 10
  
  # the medians of the errors from ten simulations are stored here
  errors_twoPass <- numeric(maxExponent-minExponent+1)
  errors_textBook <- numeric(maxExponent-minExponent+1)
  errors_R <- numeric(maxExponent-minExponent+1)
  
  # get number of cores used for parallel computations
  cores <- detectCores()
  
  # initialize the cluster
  cluster <- makeCluster(cores)
  clusterExport(cluster, c("varAnalytical", "varTextBook", "varTwoPass"))
  
  # simulation loop
  for (exponent in minExponent:maxExponent) {
    # compute a matrix consisting of the intermediate errors (do this in parallel)
    tmp_errors <- parSapply(cluster, 1:simCount, function(x){
      # create a random permutation of the integers from 1:2^maxExponent
      permutation <- sample.int(2**exponent)
      
      # calculate the real variance and the errors of the algorithms
      realVar <- varAnalytical(permutation)
      return(c(abs(varTwoPass(permutation) - realVar),
               abs(varTextBook(permutation) - realVar),
               abs(var(permutation) - realVar)))
    })
    
    # calculate the medians and store them
    errors_twoPass[exponent-minExponent+1] <- median(tmp_errors[1,])
    errors_textBook[exponent-minExponent+1] <- median(tmp_errors[2,])
    errors_R[exponent-minExponent+1] <- median(tmp_errors[3,])
  }
  
  # stop the cluster
  stopCluster(cluster)
  
  # create the dataframe used for plotting
  res <- data.frame(exponent=rep(minExponent:maxExponent, 3),
                    type=rep(c('twoPass', 'textBook', 'R'),
                             each=(maxExponent-minExponent+1)),
                    error=c(errors_twoPass, errors_textBook, errors_R))

  # plot the results
  print(
    ggplot(res, aes(x=exponent, y=error, color=type)) +
      geom_line(size=2, alpha=0.5) + 
      ggtitle("Median Errors of Variance Algorithms") +
      scale_y_continuous(name="Median Absolute Error of Ten Simulations") + 
      scale_x_continuous(name="Binary Logarithm of Vector Size",
                         breaks=minExponent:maxExponent) +
      theme(plot.title = element_text(hjust = 0.5))
  )
}

#print(microbenchmark(varEval(), times=1))

# The algorithms produce exact results when the exponent is in the interval [1, 21].
#
# As seen in the plot, R most likely uses the two pass algorithm.
#
# The order of the numbers affects the accuracy of the computations. If a vector
# is sorted in ascending order, the variance calculations will be more accurate
# than for vectors that are sorted in descending order.
# This is due to the fact that adding small values to a big sum
# (like it is done when computing the variance)
# does not affect the big sum, but adding small values successively gradually
# increases the sum without disturbing the overall accuracy.