# Name: Christian Peters

library(testthat)

# No. 1)
#=======

greville <- function(x) {
  x_inv <- matrix(0, nrow=ncol(x), ncol=nrow(x))
  
  x_inv[,1] <- 1/sum(x[1,]*x[1,])*x[1,]
  
  if (ncol(x) == 1 && nrow(x) == 1) {
    return(x_inv)
  }
  
  for(j in 2:nrow(x)) {
    d <- x[j,] %*% x_inv[,1:(j-1)]
    c <- x[j,] - d %*% x[1:(j-1),]
    b <- 1/sum(c*c)*t(c)
    x_inv[,1:(j-1)] <- x_inv[,1:(j-1)] - b %*% d
    x_inv[,j] <- b
  }
  
  return(x_inv)
}

test_matrixInversion <- function(algo) {
  expect_equal(algo(matrix(1)), matrix(1))
  expect_equal(algo(matrix(2)), matrix(0.5))
  
  a <- matrix(sample.int(10*10), nrow=10)
  expect_equal(a %*% algo(a), diag(10))
  expect_equal(a %*% algo(a) %*% a, a)
  
  a <- matrix(sample.int(100*3), nrow=3)
  expect_equal(a %*% algo(a), diag(3))
  expect_equal(a %*% algo(a) %*% a, a)
}

test_matrixInversion(greville)