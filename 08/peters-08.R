# Name: Christian Peters

library(testthat)

# No. 1)
# ======

#' Implementation of the quadratic interpolation optimization algorithm.
#' 
#' @param f     A function that is to be minimized.
#' @param left  The left endpoint of the search interval. It must be numeric.
#' @param right The right endpoint of the search interval. It must be numeric.
#' @param eps   Tolerance threshold. Used to adjust the precision of the optimizer.
quadratic_interpolation <- function(f, left, right, eps=1e-4) {
  # type checks
  stopifnot(is.function(f), is.numeric(c(left, right)), left<=right, is.numeric(eps))
  
  # calculate initial values
  best <- (left+right)/2
  f_right <- f(right)
  f_left <- f(left)
  f_best <- f(best)
  new <- 1/2 * (f_right * (left**2 - best**2) + f_best*(right**2-left**2) +
                  f_left * (best**2 - right**2)) / (f_right*(left-best) + 
                                                    f_best*(right-left) + 
                                                    f_left * (best-right))
  # iterate until the desired precision is reached
  while (abs(right - left) > eps) {
    f_new <- f(new)
    if (f_new < f_best) {
      tmp <- new
      new <- best
      best <- tmp
      tmp <- f_new
      f_new <- f_best
      f_best <- tmp
    }
    if (new < best) {
      left <- new
      f_left <- f_new
    } else {
      right <- new
      f_right <- f_new
    }
    new_old <- new
    new <- 1/2 * (f_right * (left**2 - best**2) + f_best*(right**2 - left**2) +
                    f_left * (best**2 - right**2)) / (f_right*(left-best) + 
                                                      f_best*(right-left) + 
                                                      f_left * (best-right))
    
    # Also check if new is NA. This occurs when dividing by zero when computing
    # 'new' and can happen if 'left' or 'right' is the minimum.
    if (is.na(new) || isTRUE(all.equal(new, new_old)) || new <= left || new >= right) {
      return(best)
    }
  }
  return(best)
}

#' This function is used to test the quadratic interpolation algorithm.
test_quadratic_interpolation <- function() {
  test_that("Test of the quadratic interpolation algorithm.", {
    # simple quadratic function
    expect_equal(quadratic_interpolation(function(x) x**2, -1, 1), 0, tolerance=1e-4)
    # shifted 4'th degree polynomial
    expect_equal(quadratic_interpolation(function(x) (x-3)**4 + 1, 1, 4), 3, tolerance=1e-4)
    # find pi/2
    expect_equal(quadratic_interpolation(function(x) -sin(x), 1, 2), pi/2, tolerance=1e-4)
    # additional logarithm
    expect_equal(quadratic_interpolation(function(x) log((x+2)**2)+3, -10, 10), -2, tolerance=1e-4)
  })
}

test_quadratic_interpolation()