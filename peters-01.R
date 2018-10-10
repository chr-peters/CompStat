library(testthat)

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
  stopifnot(a %% 1 == 0, b %% 1 == 0, a > 0, b > 0)
  
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