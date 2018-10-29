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
    expect_false(is.unsorted(sortFunc(sample(c(-1, -2, 1, 2, 3, 4, 5), 100,
                                             replace=TRUE))[[1]]))
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

## runTuringMachine - runs the simulator for a generic turing machine
##
## Input
##   I       - the instruction set. Must be a list, that contains a 
##             list for every instruction. Every single instruction is
##             a list with 5 entries (in this order)
##             * The cell number
##             * The character from the tape
##             * The character, which will be printed on the tape 
##             * The instruction for the cell selector, must be "L",
##               "R" or "N"
##             * The next cell number
##   tape    - the initial tape of the TM
##   L0      - the initial instruction of the TM
##   F       - the final states of the TM
##   ini.pos - the initiel position of the RW
##
## Output
##   The finale tape is returned invisible. As a side effect after every step
##   the tape is  printed on the R console
runTuringMachine <- function(I, tape, L0, F, ini.pos) {
  ## Save the current of the TM in a list. This list contains:
  ## cell.selector - the current value of the cell selector
  ## tape          - the current tape
  ## pos           - the current position of the RW
  tm <- list(cell.selector = L0, tape = tape, pos = ini.pos)
  repeat {
    ## some nice output on the console
    tape.str <- paste(" ", paste(tm$tape, collapse = " "), sep = "")
    i <- tm$pos * 2 
    substr(tape.str, i - 1, i - 1) = "*"
    substr(tape.str, i + 1, i + 1) = "*"
    cat("tape:", tape.str , " cell.selector:", tm$cell.selector)
    ## error, if the machine runs of the tape
    if (tm$pos < 1 || tm$pos > length(tape))
      stop("The TM ran of the tape")
    ## stop if final state is reaches
    if(tm$cell.selector %in% F)
      break  
    ## execute one Step
    tm <- oneStep(I, tm$tape, tm$cell.selector, tm$pos)
  }
  cat("\n")
  return(invisible(tm$tape))
}

## oneStep - executes one Step of a given Turing Machine.
##
## Input
##   I             - the instruction set, as described in the 
##                   documentation of turingMachine
##   tape          - the current tape of the TM
##   cell.selector - (cell selector) the current state of the TM
##   pos           - the current position of the RW
##
## Output
##   A list with the named elements
##   tape  - the tape of the TM after the execution of the step
##   state - the new state of the TM
##   pos   - the new position of the RW
oneStep <- function(I, tape, cell.selector, pos) {
  ## 1st cycle
  command.register <- Filter(function(p) p[[1]] == cell.selector, I)
  ## 2nd cycle
  input.register <- tape[pos]
  ## 3rd cycle
  tmp <- Find(function(p) p[[2]] == input.register, command.register)
  operation.register <- tmp[3:4]
  cell.selector <- tmp[[5]]
  ## some nice output
  cat("  next:", unlist(tmp), "\n")
  ## 4th cycle
  tape[pos] <- operation.register[[1]]
  ## 5th cycle
  pos <- if(operation.register[[2]] == "L") pos - 1L else pos
  pos <- if(operation.register[[2]] == "R") pos + 1L else pos
  ## return the new TM state
  return(list(cell.selector = cell.selector, tape = tape, pos = pos))
}

#' Thre program of the turing machine.
#' 
#' In this program, the additional character 'd' is introduced, which indicates
#' that either an x or y was present at this place in the beginning.
#' This character is only temporary though and removed during the program
#' execution.
prog <- list(
  # starting state
  list(1, "b", "b", "R", 2),
  # search for a y
  list(2, "y", "d", "L", 3),
  list(2, "x", "x", "R", 2),
  list(2, "d", "d", "R", 2),
  list(2, "b", "b", "L", 7),
  # back to starting position (look for x afterwards)
  list(3, "x", "x", "L", 3),
  list(3, "d", "d", "L", 3),
  list(3, "b", "b", "R", 4),
  # search for first x
  list(4, "y", "y", "R", 4),
  list(4, "d", "d", "R", 4),
  list(4, "x", "d", "R", 5),
  # search for second x
  list(5, "y", "y", "R", 5),
  list(5, "d", "d", "R", 5),
  list(5, "x", "d", "L", 6),
  list(5, "b", "b", "L", 7), # no x left, y is the result
  # back to starting position (look for y afterwards)
  list(6, "y", "y", "L", 6),
  list(6, "d", "d", "L", 6),
  list(6, "b", "b", "R", 2),
  # check if any x are left, otherwise y is the result
  list(7, "d", "b", "L", 7),
  list(7, "b", "y", "N", 0),
  list(7, "x", "b", "L", 8),
  list(7, "y", "b", "L", 7),
  # x is the result
  list(8, "x", "b", "L", 8),
  list(8, "d", "b", "L", 8),
  list(8, "b", "x", "N", 0)
)
# test cases
ini.tape.1 <- c("b", "b")
ini.tape.2 <- c("b", "x", "y", "x", "b")
ini.tape.3 <- c("b", "x", "y", "x", "x", "x", "x", "y", "b")
ini.tape.4 <- c("b", "y", "x", "b")
ini.tape.5 <- c("b", "y", "x", "y", "y", "x", "x", "b")

# execute test cases
print("Example 1, expected result: *y*b")
runTuringMachine(prog, tape = ini.tape.1, L0 = 1L, F = 0L,
                 ini.pos = 1L)
cat("\n")
print("Example 2, expected result: *y*b b b b")
runTuringMachine(prog, tape = ini.tape.2, L0 = 1L, F = 0L,
                 ini.pos = 1L)
cat("\n")
print("Example 3, expected result: *x*b b b b b b b b")
runTuringMachine(prog, tape = ini.tape.3, L0 = 1L, F = 0L,
                 ini.pos = 1L)
cat("\n")
print("Example 4, expected result: *y*b b b")
runTuringMachine(prog, tape = ini.tape.4, L0 = 1L, F = 0L,
                 ini.pos = 1L)
cat("\n")
print("Example 5, expected result: *y*b b b b b b b")
runTuringMachine(prog, tape = ini.tape.5, L0 = 1L, F = 0L,
                 ini.pos = 1L)

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