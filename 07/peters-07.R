# Name: Christian Peters

library(testthat)

givens <- function(x) {
  q <- diag(nrow(x))
  
  n <- ncol(x)
  m <- nrow(x)
  
  for (j in 1:n) {
    if (j+1 > m) {
      break
    }
    for (i in (j+1):m) {
      J <- diag(2)
      if (isTRUE(all.equal(x[j, j], 0)) && isTRUE(all.equal(x[i, j], 0))) {
        # do nothing
      } else {
        a <- numeric(1)
        b <- numeric(1)
        if (abs(x[i, j]) > abs(x[j, j])) {
          a <- 1/sqrt(1 + (x[j, j]/x[i, j])**2)
          b <- (x[j, j]/x[i, j]) * a
        } else {
          b <- 1/sqrt(1+(x[i, j]/x[j, j])**2)
          a <- (x[i, j]/x[j, j]) * b
        }
        J <- matrix(c(b, -a, a, b), nrow=2)
      }
      
      x[c(j, i),] <- J %*% x[c(j, i),]
      q[c(j, i),] <- J %*% q[c(j, i),]
    }
  }
  
  if (n+1 < m) {
    x <- x[-((n+1):m),]
    q <- q[-((n+1):m),]
  }
  
  return(backsolve(x, q))
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
    
  })
}

test_matrixInversion(givens)