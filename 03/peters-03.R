# Name: Christian Peters

#' Vectorized implementation of the simple quick sort algorithm.
#' 
#' @param x A vector to sort. It must be numeric and only contain finite elements.
#' 
#' @return A list containing the sorted vector as the first element and the number
#'         of comparisons as the second element.
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
    expect_false(is.unsorted(sortFunc(sample(c(-1, -2, 1, 2, 3, 4, 5), 100, replace=TRUE))))
    # sorting some samples ascending as well as descending
    expect_false(is.unsorted(sortFunc(sample.int(100))[[1]]))
    expect_false(is.unsorted(sortFunc(rnorm(100))[[1]]))
  })
}

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