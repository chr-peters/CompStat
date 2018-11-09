# Name: Christian Peters

library(testthat)
library(ggplot2)

# No. 1)
#=======

#' Calculates the sample variance of a vector x using the Youngs and Cramer
#' updating algorithm.
#' 
#' @param x The vector to calculate the variance of. It must be numeric.
#' 
#' @return The variance of x.
youngsCramerVar <- function(x) {
  stopifnot(is.numeric(x), length(x) > 0)
  
  n <- length(x)
  # The variance of a vector with one element is zero.
  if (n==1){
    return(0);
  }
  
  # Implementation of the Youngs Cramer formulas.
  t <- x[1]
  s <- 0
  for(j in 2:n) {
    t <- t + x[j]
    s <- s + (j*x[j]-t)**2/(j*(j-1))
  }
  
  # Return the sample variance.
  return(s/(n-1))
}

#' This function contains various test cases for the youngsCramerVar function.
test_youngsCramerVar <- function() {
  test_that("Test of youngsCramerVar function", {
    # edge cases
    expect_error(youngsCramerVar(numeric()))
    expect_equal(youngsCramerVar(rnorm(1)), 0)
    # minimal example
    expect_true(all.equal(youngsCramerVar(c(0, 1)), 0.5))
    # some bigger examples
    a <- rnorm(100)
    # use all.equal to test "near equality"
    expect_true(all.equal(youngsCramerVar(a), var(a)))
    a <- sample.int(100)
    expect_true(all.equal(youngsCramerVar(a), var(a)))
    a <- rnorm(10000)
    expect_true(all.equal(youngsCramerVar(a), var(a)))
  })
}

test_youngsCramerVar()

# No. 2)
#=======

#' Computes the variance of a vector using the two pass algorithm.
#' 
#' @param x A numeric vector.
#' @param m An optional parameter indicating the sample size based on which
#'          the shift value is calculated (shift = mean of sample of size m).
#'          It must be greater than or equal to zero.
#' 
#' @return The variance of x.
varTwoPass <- function(x, m=0) {
  if (m > 0) {
    x <- x - mean(sample(x, m))
  }
  return(sum((x-mean(x))**2)/(length(x)-1))
}

#' Computes the variance of a vector using the textbook algorithm.
#' 
#' @param x A numeric vector.
#' @param m An optional parameter indicating the sample size based on which
#'          the shift value is calculated (shift = mean of sample of size m).
#'          It must be greater than or equal to zero.
#' 
#' @return The variance of x.
varTextBook <- function(x, m=0) {
  stopifnot(is.numeric(x), m>=0)
  if (m > 0) {
    x <- x - mean(sample(x, m))
  }
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

#' This function visualizes the results of several shifting parameters.
startSimulation <- function(){
  # first determine the real variance
  realVariance <- varAnalytical(1:1000)
  
  # now get the results of twoPass and textBook for different shift sample sizes
  resultsTwoPass <- sapply(0:10, function(m){
    return(varTwoPass(1:1000 + 1e16, m))
  })
  resultsTextBook <- sapply(0:10, function(m){
    return(varTextBook(1:1000 + 1e16, m))
  })
  
  # get the absolute errors
  absTwoPass <- abs(realVariance-resultsTwoPass)
  absTextBook <- abs(realVariance-resultsTextBook)
  
  # create data.frame used for plotting
  df <- data.frame(m=rep(0:10, 2), type=rep(c("TwoPass", "TextBook"), each=11),
                   error=c(absTwoPass, absTextBook))
  
  # plot the results
  print (
    ggplot(df, aes(x=m, y=error, color=type)) + 
      geom_line(size=1.5) +
      ggtitle("The Effect of Shifting on Variance Calculation") +
      scale_x_continuous("m (sample size used for shifting)", breaks=0:10) +
      scale_y_continuous("absolute error") +
      theme(plot.title = element_text(hjust = 0.5))
  )
}

startSimulation()

# In this case no improvements in precision are gained by using a shift parameter
# greater than one. There is however a striking difference in accuracy between
# m=1 and m=0 (not shifting at all).
#
# Based on these results, I recommend a value of m=1 because any further gain in
# precision when estimating the mean does not necessarily lead to further
# improvements when calculating the variance as such.
# It is therefore sufficient to end up at a rough estimate of the numerical
# scale we are dealing with and m=1 should usually be an adequate solution to
# accomplish this.