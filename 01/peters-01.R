library(testthat)
library(microbenchmark)

#' This function determines the greatest common divisor of two natural numbers
#' using euclid's algorithm.
#' 
#' In each iteration of the algorithm, the remainder of the two numbers is
#' calculated. If it is zero, the smaller of the two numbers is the correct
#' result. If it is not zero, the remainder is the next candidate.
#' 
#' @param a The first natural number. It has to be greater than zero.
#' @param b The second natural number. It has to be greater than zero.
#' 
#' @return The greatest common divisor of the natural numbers a and b.
ggT <- function(a, b) {
  # test if a and b are both integers and greater than zero
  stopifnot(is.numeric(a), is.numeric(b),
            a %% 1 == 0, b %% 1 == 0, a > 0, b > 0)
  
  # euclid's algorithm
  remainder <- a %% b
  while (remainder != 0) {
    a <- b
    b <- remainder
    remainder <- a %% b
  }
  return(b)
}

#' This function tests the results of the ggT-function using a few examples.
test_ggT <- function() {
  test_that("Euclidean ggT test", {
    # this leads to an error because both numbers have to be greater than zero
    expect_error(ggT(0, 1))
    expect_equal(ggT(1, 2), 1)
    expect_equal(ggT(26, 34), 2)
    expect_equal(ggT(58, 145), 29)
    expect_equal(ggT(1000001, 1048576), 1)
  })
}

test_ggT()

#' This function computes the n'th fibonacci number using a naive recursive
#' approach.
#' 
#' @param n An integer describing which element of the fibonacci sequence to
#' return.
#' 
#' @return Returns the n'th fibonacci number or n, if n is less than 2.
fiboSimple <- function(n) {
  stopifnot(is.numeric(n), n %% 1 == 0)
  if (n < 2) {
    return(n)
  }
  # definition of the fibonacci sequence
  return(fiboSimple(n-1) + fiboSimple(n-2))
}

#' This function tests the simple fibonacci implementation on some examples.
test_fiboSimple <- function() {
  test_that("Simple fibonacci test", {
    expect_equal(fiboSimple(1), 1)
    expect_equal(fiboSimple(2), 1)
    expect_equal(fiboSimple(3), 2)
    expect_equal(fiboSimple(12), 144)
    expect_equal(fiboSimple(15), 610)
  })
}

test_fiboSimple()

#' This function computes the n'th fibonacci number using a more efficient
#' approach.
#' 
#' In every iteration of the algorithm, only the two most recent fibonacci
#' numbers are kept in memory. The next number of the sequence is calculated
#' based on these stored numbers in each step. This way, no fibonacci number has
#' to be computed twice as it is the case with the simple approach.
#' 
#' @param n An integer describing which element of the fibonacci sequence to
#' return.
#' 
#' @return Returns the n'th fibonacci number or n, if n is less than 2.
fiboEfficient <- function(n) {
  stopifnot(is.numeric(n), n %% 1 == 0)
  if (n < 2) {
    return(n)
  }
  # these are the two most recently calculated fibonacci numbers
  a <- 1
  b <- 1
  # now iterate until the desired fibonacci number is reached
  while (n > 2) {
    # calculate the next fibonacci number
    c <- a + b
    # update the intermediate results
    a <- b
    b <- c
    # decrement the counter
    n <- n - 1
  }
  # return the most recent fibonacci number
  return(b)
}

#' This function tests the efficient fibonacci implementation on some examples.
test_fiboEfficient <- function() {
  test_that("Efficient fibonacci test", {
    expect_equal(fiboEfficient(1), 1)
    expect_equal(fiboEfficient(2), 1)
    expect_equal(fiboEfficient(3), 2)
    expect_equal(fiboEfficient(12), 144)
    expect_equal(fiboEfficient(15), 610)
    expect_equal(fiboEfficient(30), 832040)
  })
}

test_fiboEfficient()

# performance benchmark
n <- 15 # fibonacci numbers to generate
times <- 20 # used to average runtime
runtimesSimple <- vector(mode="numeric", length=n+1)
runtimesEfficient <- vector(mode="numeric", length=n+1)
for (i in 0:n) {
  runtimesSimple[i+1] <- sum(microbenchmark(fiboSimple(i), times=times)$time)/times
  runtimesEfficient[i+1] <- sum(microbenchmark(fiboEfficient(i), times=times)$time)/times
}

plot(0:n, runtimesSimple, main="Runtimes of naive Fibonacci implementation",
     xlab="Fibonacci number", ylab="Runtime in nanoseconds")
plot(0:n, runtimesEfficient, main="Runtimes of efficient Fibonacci implementation",
     xlab="Fibonacci number", ylab="Runtime in nanoseconds")