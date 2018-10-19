# Name: Christian Peters

library(testthat)

#' Implementation of the bubble sort algorithm.
#' 
#' @param a The vector to sort. It has to be numeric.
#' @param decreasing Optional parameter describing whether to sort a in
#' descending order or not. The default value is FALSE which means that a will
#' be sorted in ascending order.
#' 
#' @return A list containing the sorted version of a as the first element and
#' the number of comparisons as the second element.
bubbleSort <- function(a, decreasing=FALSE) {
  stopifnot(is.numeric(a))
  count <- 0
  # no need to sort if there are 0 or 1 elements
  if (length(a) <= 1) {
    return(list(a, count))
  }
  # j is the last index to be compared in each run
  # the indices greater than j+1 are already sorted and don't need comparison
  for (j in (length(a)-1):1) {
    # swapped is used to stop early once the vector is sorted
    swapped <- FALSE
    # i is the index of the element that is currently compared
    for (i in 1:j) {
      count <- count + 1
      # carry out the swap if necessary
      if (a[i] > a[i+1] && !decreasing || a[i] < a[i+1] && decreasing) {
        tmp <- a[i+1]
        a[i+1] <- a[i]
        a[i] <- tmp
        swapped <- TRUE
      }
    }
    # stop if the vector is sorted
    if (!swapped) {
      break
    }
  }
  return(list(a, count))
}

#' This function is used to test the bubbleSort function.
test_bubbleSort <- function() {
  test_that("Test function bubbleSort", {
    # edge cases
    expect_error(bubbleSort(c()))
    expect_equal(bubbleSort(rnorm(0))[[1]], numeric())
    expect_equal(bubbleSort(c(0))[[1]], c(0))
    expect_equal(bubbleSort(c(0))[[2]], 0)
    # minimal examples
    expect_equal(bubbleSort(c(0, 1))[[1]], c(0, 1))
    expect_equal(bubbleSort(c(0, 1))[[2]], 1)
    expect_equal(bubbleSort(c(1, 0))[[1]], c(0, 1))
    expect_equal(bubbleSort(c(1, 0))[[2]], 1)
    expect_equal(bubbleSort(c(-1, -2))[[1]], c(-2, -1))
    # sorting some samples ascending as well as descending
    expect_false(is.unsorted(bubbleSort(sample.int(100))[[1]]))
    expect_false(is.unsorted(bubbleSort(rnorm(100))[[1]]))
    expect_false(is.unsorted(rev(bubbleSort(sample.int(100), TRUE)[[1]])))
    expect_false(is.unsorted(rev(bubbleSort(rnorm(100), TRUE)[[1]])))
  })
}

test_bubbleSort()

bucketSort <- function(a, n) {
  stopifnot(is.numeric(a), is.numeric(n), n > 0)
  maxValue <- max(a)
  minValue <- min(a) - 5 * .Machine$double.eps
  buckets <- numeric(length=length(a))
  count <- 0
  for (i in 1:length(a)) {
    buckets[i] <- ceiling(n * (a[i] - minValue) / (maxValue - minValue))
  }
  res <- numeric()
  for (j in 1:n) {
    tmp <- bubbleSort(a[buckets==j])
    res <- c(res, tmp[[1]])
    count <- count + tmp[[2]]
  }
  return(list(res, count))
}

test_bucketSort <- function() {
  test_that("Test function bucketSort", {
    # edge cases
    expect_error(bucketSort(c()))
    # minimal examples
    expect_equal(bucketSort(c(0, 1), 1)[[1]], c(0, 1))
    expect_equal(bucketSort(c(0, 1), 1)[[2]], 1)
    expect_equal(bucketSort(c(1, 0), 1)[[1]], c(0, 1))
    expect_equal(bucketSort(c(1, 0), 1)[[2]], 1)
    expect_equal(bucketSort(c(-1, -2), 1)[[1]], c(-2, -1))
    # sorting some samples
    expect_false(is.unsorted(bucketSort(sample.int(100), 10)[[1]]))
    expect_false(is.unsorted(bucketSort(rnorm(100), 10)[[1]]))
  })
}

test_bucketSort()
