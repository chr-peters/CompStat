# Name: Christian Peters

library(testthat)

#' Vectorized implementation of the simple quick sort algorithm.
#' 
#' @param x A vector to sort. It must be numeric and only contain finite elements.
#' 
#' @return A list containing the sorted vector as the first element and the number
#'         of comparisons as the second element.
#'         Both elements are of the 'numeric' data type.
simpleQSort <- function(x) {
  # Check if the input is of the correct type
  stopifnot(
    is.vector(x),
    is.numeric(x),
    all(is.finite(x))
  )
  
  # No sorting necessary
  if (length(x) <= 1L)
    return(list(res = x, ncmp = 0))
  
  # Select a random pivot element
  pivot <- sample(x, 1)
  
  # Determine the elements that are less than, greater than or equal to the pivot
  # element
  less_than_pivot <- x < pivot
  equal_to_pivot <- x == pivot
  greater_than_pivot <- x > pivot
  
  # Execute quick sort recursively on the left as well as on the right partition
  # of the list
  left <- Recall(x[less_than_pivot])
  center <- x[equal_to_pivot]
  right <- Recall(x[greater_than_pivot])
  
  # Concatenate the results
  list(res = c(left$res, center, right$res),
       ncmp = left$ncmp + right$ncmp + length(x))
}

#' This function is used to test a generic sort function.
test_sortGeneric <- function(sortFunc) {
  test_that("Test sorting function", {
    # edge cases
    expect_error(sortFunc(c()))
    expect_equal(sortFunc(numeric())[[1]], numeric())
    expect_equal(sortFunc(c(0))[[1]], c(0))
    # minimal examples
    expect_equal(sortFunc(c(0, 1))[[1]], c(0, 1))
    expect_equal(sortFunc(c(1, 0))[[1]], c(0, 1))
    expect_equal(sortFunc(c(-1, -2))[[1]], c(-2, -1))
    # repetitions
    expect_false(is.unsorted(sortFunc(sample(c(-1, -2, 1, 2, 3, 4, 5), 100, replace=TRUE))[[1]]))
    # sorting some bigger samples
    expect_false(is.unsorted(sortFunc(sample.int(100))[[1]]))
    expect_false(is.unsorted(sortFunc(rnorm(100))[[1]]))
  })
}

# test the simple quicksort
test_sortGeneric(simpleQSort)

#' This function is used to sort a three element numerical vector.
sortThree <- function(x) {
  # Check if the input is of the correct type
  stopifnot(
    is.vector(x),
    is.numeric(x),
    all(is.finite(x)),
    length(x)==3
  )
  # sort the three element list by conducting the necessary comparisons
  if (x[1] > x[2]) {
    tmp <- x[1]
    x[1] <- x[2]
    x[2] <- tmp
  }
  if (x[2] > x[3]) {
    tmp <- x[2]
    x[2] <- x[3]
    x[3] <- tmp
    if (x[1] > x[2]) {
      tmp <- x[1]
      x[1] <- x[2]
      x[2] <- tmp
    }
  }
  return(x)
}

test_sortThree <- function() {
  test_that('Test three element sort', {
    expect_equal(sortThree(c(0, 0, 0)), c(0, 0, 0))
    expect_false(is.unsorted(sortThree(c(1, 2, 3))))
    expect_false(is.unsorted(sortThree(c(1, 3, 2))))
    expect_false(is.unsorted(sortThree(c(2, 1, 3))))
    expect_false(is.unsorted(sortThree(c(2, 3, 1))))
    expect_false(is.unsorted(sortThree(c(3, 1, 2))))
    expect_false(is.unsorted(sortThree(c(3, 2, 1))))
    expect_false(is.unsorted(sortThree(rnorm(3))))
    expect_false(is.unsorted(sortThree(sample.int(3))))
  })
}

# test three element sort function
test_sortThree()

#' Vectorized implementation of the clever quick sort algorithm.
#' 
#' @param x A vector to sort. It must be numeric and only contain finite elements.
#' 
#' @return A list containing the sorted vector as the first element and the number
#'         of comparisons as the second element.
#'         Both elements are of the 'numeric' data type.
cleverQSort <- function(x) {
  # Check if the input is of the correct type
  stopifnot(
    is.vector(x),
    is.numeric(x),
    all(is.finite(x))
  )
  
  # No sorting necessary
  if (length(x) <= 1L)
    return(list(res = x, ncmp = 0))
  
  # select the pivot as the median of three (first, middle, last element)
  medianCandidates <- c(x[1], x[length(x) %/% 2], x[length(x)])
  pivot <- sortThree(medianCandidates)[2]
  
  # Determine the elements that are less than, greater than or equal to the pivot
  # element
  less_than_pivot <- x < pivot
  equal_to_pivot <- x == pivot
  greater_than_pivot <- x > pivot
  
  # Execute quick sort recursively on the left as well as on the right partition
  # of the list
  left <- Recall(x[less_than_pivot])
  center <- x[equal_to_pivot]
  right <- Recall(x[greater_than_pivot])
  
  # Concatenate the results
  list(res = c(left$res, center, right$res),
       ncmp = left$ncmp + right$ncmp + length(x))
}

test_sortGeneric(cleverQSort)

# No. 3:
#
# a) (9.8)_10 = (1001.11001100...)_2 -> (0.10100e(10000-1000))_2 = (0.10100e100)_2
#
# b) 1. Compare exponents -> e_w = (100)_2
#    2. Add mantissas:
#                        ( 0 . 1 0 1 0 0 )_2
#                      + ( 0 . 1 0 1 0 0 )_2
#                      ---------------------
#                      = ( 1 . 0 1 0 0 0 )_2
#    3. Normalize: (1.01000e100)_2 -> (0.101000e101)_2
#    4. Round: (0.101000e101)_2 -> (0.10100e101)_2
#
# c) (19.6)_10 = (10011.10011001...)_2 -> (0.10100e101)_2
#
# d) The rounding errors occured when converting the number (9.8)_10
#    to the inexact binary representation. When performing the addition,
#    no rounding errors happened, but due to the fact that the inputs were no
#    exact representations of the number (9.8)_10, the result was inexact as well.