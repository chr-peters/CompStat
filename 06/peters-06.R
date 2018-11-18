# Name: Christian Peters

library(testthat)

# No. 1)
#=======

#' Implementation of greville's algorithm for matrix inversion.
#' 
#' @param x The matrix to be inverted.
#' 
#' @return The inverted matrix.
greville <- function(x) {
  stopifnot(is.numeric(x), is.matrix(x))
  
  # the algorithm only works if x has a full row rank
  transposed <- FALSE
  if (nrow(x) > ncol(x)) {
    # it is ok to transpose x first and transpose it back after inversion,
    # because transposing and inverting a matrix is a commutative operation
    x <- t(x)
    transposed <- TRUE
  }
  
  # initialize the result
  x_inv <- matrix(0, nrow=ncol(x), ncol=nrow(x))
  
  # first iteration
  x_inv[,1] <- 1/sum(x[1,]*x[1,])*x[1,]
  
  # no further need for inversion
  if (nrow(x) == 1) {
    return(x_inv)
  }
  
  # greville's algorithm
  for(j in 2:nrow(x)) {
    d <- x[j,] %*% x_inv[,1:(j - 1)]
    c <- x[j,] - d %*% x[1:(j - 1),]
    b <- 1 / sum(c * c) * t(c)
    x_inv[,1:(j - 1)] <- x_inv[,1:(j - 1)] - b %*% d
    x_inv[,j] <- b
  }
  
  # transpose back if necessary
  if (transposed) {
    x_inv <- t(x_inv)
  }
  
  return(x_inv)
}

#' This function can be used to test a matrix inversion algorithm.
test_matrixInversion <- function(algo) {
  test_that("Test of matrix inversion algorithm", {
    # edge cases
    expect_equal(algo(matrix(1)), matrix(1))
    expect_equal(algo(matrix(2)), matrix(0.5))
    
    # square matrix
    a <- matrix(sample.int(10*10), nrow=10)
    expect_equal(algo(a) %*% a, diag(10))
    expect_equal(algo(algo(a)), a)
    expect_equal(a %*% algo(a) %*% a, a)
    
    # more rows than columns
    a <- matrix(sample.int(100*3), nrow=100)
    expect_equal(algo(a) %*% a, diag(3))
    expect_equal(algo(algo(a)), a)
    expect_equal(a %*% algo(a) %*% a, a)
    
    # more columns than rows
    a <- matrix(sample.int(100*3), nrow=3)
    expect_equal(a %*% algo(a), diag(3))
    expect_equal(algo(algo(a)), a)
    expect_equal(a %*% algo(a) %*% a, a)
  })
}

test_matrixInversion(greville)