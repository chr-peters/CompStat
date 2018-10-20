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
    # sorting some samples ascending as well as descending
    expect_false(is.unsorted(sortFunc(sample.int(100))[[1]]))
    expect_false(is.unsorted(sortFunc(rnorm(100))[[1]]))
    expect_false(is.unsorted(rev(sortFunc(sample.int(100), TRUE)[[1]])))
    expect_false(is.unsorted(rev(sortFunc(rnorm(100), TRUE)[[1]])))
  })
}

test_sortGeneric(bubbleSort)

#' Implementation of the bucket sort algorithm.
#' 
#' Complexity: O(n + buckets*(n/buckets)^2)
#' 
#' @param a The vector to sort. It has to be numeric.
#' @param n The number of buckets. It has to be greater than zero.
#' 
#' @return A list containing the sorted version of a as the first element and
#' the number of comparisons as the second element.
bucketSort <- function(a, n) {
  stopifnot(is.numeric(a), is.numeric(n), n > 0)
  
  # no need to sort
  if (length(a) <= 1) {
    return(list(a, 0))
  }
  
  # calculate maximum and minimum values
  maxValue <- max(a)
  # subtract a very small value to avoid the 0 bucket
  minValue <- min(a) - 5 * .Machine$double.eps
  buckets <- numeric(length=length(a))
  count <- 0
  
  # determine the bucket for each element
  for (i in 1:length(a)) { # complexity: O(n)
    buckets[i] <- ceiling(n * (a[i] - minValue) / (maxValue - minValue))
  }
  
  # apply bubbleSort to each bucket
  res <- numeric()
  for (j in 1:n) {
    tmp <- bubbleSort(a[buckets==j]) # complexity: O(n^2)
    res <- c(res, tmp[[1]])
    count <- count + tmp[[2]]
  }
  return(list(res, count))
}

# now test bucket sort
for (buckets in 1:20) {
  test_sortGeneric(function (a, descending=FALSE) {
    if (descending) {
      # reverse the result if the testing function expects a descending list
      res <- bucketSort(a, buckets)
      return(list(rev(res[[1]]), res[[2]]))
    } else {
      return(bucketSort(a, buckets))
    }
  })
}

# simulation function
startSimulation <- function() {
  # simulation parameters
  numBuckets <- 100
  numRuns <- 100
  
  buckets <- 1:numBuckets
  # for every bucket count, the mean of the comparisons is stored here
  meanComparisons <- numeric(length = numBuckets)
  
  for (n in buckets) {
    # the results of each run are stored here
    comparisons <- numeric(length = numRuns)
    for (i in 1:numRuns) {
      comparisons[i] <- bucketSort(sample.int(1000), n)[[2]]
    }
    # save the mean
    meanComparisons[n] <- mean(comparisons)
  }
  plot(buckets, meanComparisons, xlab="Number of Buckets",
       ylab="Mean Number of Comparisons", main="Analysis of Bucket Sort")
}

#startSimulation()

# If the elements of the list are not equally distributed (like samples from
# a normal distribution for example), some buckets contain more elements than
# others. This means that bucketsort is unfair to those buckets that contain
# more elements.
